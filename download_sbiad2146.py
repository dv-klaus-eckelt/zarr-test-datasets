#!/usr/bin/env python3
"""
Reusable download pipeline for S-BIAD2146 Xenium FFPE spatial transcriptomics slides.

Downloads .zarr.zip files from BioStudies, validates against expected sizes from the API,
and supports resume of interrupted downloads.

Usage:
    python download_sbiad2146.py --list
    python download_sbiad2146.py --tissue skin --output /tmp/slides
    python download_sbiad2146.py --slide skin_s1,kidney_s0 --output /tmp/slides
    python download_sbiad2146.py --all --output /tmp/all
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import requests
except ImportError:
    print(
        "Error: 'requests' is required. Install with: pip install requests",
        file=sys.stderr,
    )
    sys.exit(1)

try:
    from tqdm import tqdm
except ImportError:
    print(
        "Error: 'tqdm' is required. Install with: pip install tqdm",
        file=sys.stderr,
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
METADATA_PATH = SCRIPT_DIR / "sbiad2146_metadata.json"

BIOSTUDIES_API_URL = (
    "https://www.ebi.ac.uk/biostudies/api/v1/files/S-BIAD2146?start=0&length=500"
)
DOWNLOAD_BASE_URL = (
    "https://www.ebi.ac.uk/biostudies/files/S-BIAD2146"
    "/STHELAR/sdata_slides/sdata_{slide_id}.zarr.zip"
)

CHUNK_SIZE = 8192


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

def load_metadata() -> Dict[str, Dict[str, str]]:
    """Load slide metadata from the JSON config file."""
    with open(METADATA_PATH) as f:
        data = json.load(f)
    return data["slides"]


def get_tissues(metadata: Dict[str, Dict[str, str]]) -> Dict[str, List[str]]:
    """Group slide IDs by tissue type."""
    tissues: Dict[str, List[str]] = {}
    for slide_id, info in metadata.items():
        tissues.setdefault(info["tissue"], []).append(slide_id)
    return tissues


# ---------------------------------------------------------------------------
# BioStudies API
# ---------------------------------------------------------------------------

def fetch_file_sizes() -> Dict[str, int]:
    """
    Query the BioStudies file API for S-BIAD2146 and return expected file sizes.

    Returns:
        Dict mapping slide_id (e.g. ``"skin_s1"``) to file size in bytes.
    """
    resp = requests.get(BIOSTUDIES_API_URL, timeout=30)
    resp.raise_for_status()
    payload = resp.json()

    sizes: Dict[str, int] = {}
    for entry in payload["data"]:
        path: str = entry["path"]
        # e.g. "STHELAR/sdata_slides/sdata_skin_s1.zarr.zip"
        filename = os.path.basename(path)
        if filename.startswith("sdata_") and filename.endswith(".zarr.zip"):
            slide_id = filename[len("sdata_") : -len(".zarr.zip")]
            sizes[slide_id] = int(entry["Size"])
    return sizes


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def download_slide(
    slide_id: str,
    output_dir: Path,
    expected_size: Optional[int] = None,
) -> Path:
    """
    Download a single slide's ``.zarr.zip`` file with resume support.

    * If the target file already exists **and** its size matches
      *expected_size* the download is skipped.
    * If the target file exists but is **smaller** than expected (interrupted
      download) the script issues a ``Range`` request to resume from where it
      left off.
    * After completion the file size is validated against *expected_size*
      (when available).

    Args:
        slide_id: Slide identifier, e.g. ``"skin_s1"``.
        output_dir: Directory to write the file into (created if needed).
        expected_size: Expected file size in bytes.  When ``None`` size
            validation is skipped.

    Returns:
        Absolute path to the downloaded (or already present) file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = f"sdata_{slide_id}.zarr.zip"
    filepath = output_dir / filename
    url = DOWNLOAD_BASE_URL.format(slide_id=slide_id)

    # ------------------------------------------------------------------
    # Check for existing file (skip / resume)
    # ------------------------------------------------------------------
    existing_size = 0
    if filepath.exists():
        existing_size = filepath.stat().st_size
        if expected_size is not None:
            if existing_size == expected_size:
                print(f"[SKIP] {slide_id}: already downloaded ({_format_size(existing_size)})")
                return filepath
            if existing_size > expected_size:
                print(
                    f"[ERROR] {slide_id}: existing file ({existing_size} B) is "
                    f"larger than expected ({expected_size} B) -- remove it manually",
                    file=sys.stderr,
                )
                sys.exit(1)
            print(
                f"[RESUME] {slide_id}: partial download "
                f"({_format_size(existing_size)} / {_format_size(expected_size)}) "
                f"-- resuming..."
            )
        else:
            # No expected size to validate against -- trust the existing file.
            print(f"[SKIP] {slide_id}: file exists ({_format_size(existing_size)})")
            return filepath

    # ------------------------------------------------------------------
    # Perform the HTTP request
    # ------------------------------------------------------------------
    headers: Dict[str, str] = {}
    if existing_size > 0:
        headers["Range"] = f"bytes={existing_size}-"

    resp = requests.get(url, stream=True, allow_redirects=True, headers=headers, timeout=30)

    # If server doesn't support Range (200 instead of 206), start fresh
    if existing_size > 0 and resp.status_code == 200:
        print(f"[WARN] {slide_id}: server does not support Range; restarting from scratch")
        filepath.unlink()
        existing_size = 0
        mode = "wb"
        resp = requests.get(url, stream=True, allow_redirects=True, timeout=30)
    else:
        mode = "ab" if existing_size > 0 else "wb"

    resp.raise_for_status()

    # ------------------------------------------------------------------
    # Total size for progress bar
    # ------------------------------------------------------------------
    if expected_size is not None:
        total_size = expected_size
    else:
        # Content-Length in a 206 response is the remaining bytes
        content_length = resp.headers.get("Content-Length")
        if content_length is not None and existing_size > 0:
            total_size = existing_size + int(content_length)
        elif content_length is not None:
            total_size = int(content_length)
        else:
            total_size = None  # indeterminate progress

    # ------------------------------------------------------------------
    # Stream to disk
    # ------------------------------------------------------------------
    with open(filepath, mode) as f, tqdm(
        desc=slide_id,
        initial=existing_size,
        total=total_size,
        unit="B",
        unit_scale=True,
        unit_divisor=1024,
        miniters=1,
    ) as pbar:
        for chunk in resp.iter_content(chunk_size=CHUNK_SIZE):
            if chunk:
                f.write(chunk)
                pbar.update(len(chunk))

    # ------------------------------------------------------------------
    # Validate final size
    # ------------------------------------------------------------------
    final_size = filepath.stat().st_size
    if expected_size is not None and final_size != expected_size:
        print(
            f"[ERROR] {slide_id}: downloaded file size ({final_size} B) "
            f"does not match expected size ({expected_size} B)",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"[DONE] {slide_id}: {_format_size(final_size)}")
    return filepath


def download_slides(slide_ids: List[str], output_dir: Path, sizes: Dict[str, int]) -> List[Path]:
    """Download a list of slide IDs."""
    paths: List[Path] = []
    for slide_id in slide_ids:
        p = download_slide(slide_id, output_dir, expected_size=sizes.get(slide_id))
        paths.append(p)
    return paths


def download_tissue(tissue: str, output_dir: Path, sizes: Dict[str, int]) -> List[Path]:
    """Download all slides belonging to a given tissue type."""
    metadata = load_metadata()
    slide_ids = [sid for sid, info in metadata.items() if info["tissue"] == tissue]
    if not slide_ids:
        available = sorted({info["tissue"] for info in metadata.values()})
        print(
            f"[ERROR] No slides found for tissue '{tissue}'. "
            f"Available tissues: {', '.join(available)}",
            file=sys.stderr,
        )
        sys.exit(1)
    return download_slides(slide_ids, output_dir, sizes)


# ---------------------------------------------------------------------------
# Listing
# ---------------------------------------------------------------------------

def list_slides() -> None:
    """Print all 31 slides grouped by tissue with metadata."""
    metadata = load_metadata()
    tissues = get_tissues(metadata)

    total = len(metadata)
    cancerous = sum(1 for v in metadata.values() if v["patient_status"] == "cancerous")
    non_cancerous = total - cancerous

    print(f"S-BIAD2146 Slides: {total} total ({cancerous} cancerous, {non_cancerous} non-cancerous)")
    print(f"Tissue types: {len(tissues)}")
    print()

    for tissue in sorted(tissues):
        slides = sorted(tissues[tissue])
        print(f"  {tissue} ({len(slides)}):")
        for sid in slides:
            info = metadata[sid]
            print(f"    {sid:<25} {info['disease']:<55} [{info['patient_status']}]")
        print()

    # Summary line
    tissue_counts = ", ".join(f"{t}:{len(slides)}" for t, slides in sorted(tissues.items()))
    print(f"Summary: {tissue_counts}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _format_size(size_bytes: int) -> str:
    """Return a human-readable string for a byte count."""
    for unit in ("B", "KB", "MB", "GB"):
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Download slides from BioStudies S-BIAD2146 (Xenium FFPE spatial "
            "transcriptomics data)."
        ),
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List all available slides with tissue, disease, and patient status",
    )
    parser.add_argument(
        "--tissue",
        type=str,
        metavar="TISSUE",
        help="Download all slides for a tissue type (e.g., 'skin')",
    )
    parser.add_argument(
        "--slide",
        type=str,
        metavar="SLIDES",
        help="Comma-separated list of slide IDs (e.g., 'skin_s1,kidney_s0')",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Download all 31 slides",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=".",
        help="Output directory for downloaded files (default: current directory)",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # --list mode: just print and exit
    if args.list:
        list_slides()
        return

    # Any download mode requires --output (or default '.')
    if args.tissue or args.slide or args.all:
        output_dir = Path(args.output).resolve()
        print(f"Output directory: {output_dir}")

        # Fetch file sizes from the BioStudies API
        print("Fetching file sizes from BioStudies API ...")
        try:
            sizes = fetch_file_sizes()
            print(f"  Found {len(sizes)} file size entries")
        except requests.RequestException as e:
            print(f"[WARN] Could not fetch file sizes: {e}", file=sys.stderr)
            print("[WARN] Proceeding without size validation", file=sys.stderr)
            sizes = {}

        if args.tissue:
            download_tissue(args.tissue, output_dir, sizes)
        elif args.slide:
            slide_ids = [s.strip() for s in args.slide.split(",") if s.strip()]
            download_slides(slide_ids, output_dir, sizes)
        elif args.all:
            metadata = load_metadata()
            all_ids = sorted(metadata.keys())
            download_slides(all_ids, output_dir, sizes)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
