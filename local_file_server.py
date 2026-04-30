import argparse
import sys
import os
import re
from functools import partial
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

class ZarrCompatibleServer(SimpleHTTPRequestHandler):
    def __init__(self, *args, directory=None, default_file=None, **kwargs):
        self.default_file = default_file
        super().__init__(*args, directory=directory, **kwargs)

    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Range, Content-Type")
        self.send_header("Access-Control-Expose-Headers", "Content-Range, Content-Length, Accept-Ranges")
        self.send_header("Accept-Ranges", "bytes")
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(200)
        self.end_headers()

    def do_GET(self):
        if self.default_file and self.path in ("", "/"):
            self.path = f"/{self.default_file}"

        # Get the physical path of the file
        path = self.translate_path(self.path)
        
        # If it's not a file or no Range header, use default behavior
        if not os.path.isfile(path) or "Range" not in self.headers:
            return super().do_GET()

        # Parse Range header: "bytes=start-end"
        range_match = re.match(r'bytes=(\d+)-(\d+)?', self.headers['Range'])
        if not range_match:
            return super().do_GET()

        first, last = range_match.groups()
        file_size = os.path.getsize(path)
        
        start = int(first)
        end = int(last) if last else file_size - 1

        # Validation: ensure range is within file bounds
        if start >= file_size:
            self.send_error(416, "Requested Range Not Satisfiable")
            return

        length = end - start + 1
        
        # Read only the requested chunk
        with open(path, 'rb') as f:
            f.seek(start)
            content = f.read(length)

        # Send 206 Partial Content
        self.send_response(206)
        self.send_header("Content-type", "application/octet-stream")
        self.send_header("Content-Range", f"bytes {start}-{end}/{file_size}")
        self.send_header("Content-Length", str(length))
        self.end_headers()
        self.wfile.write(content)

def run_server(path: str, port: int = 8000) -> None:
    target = Path(path).expanduser().resolve()
    if not target.exists():
        print(f"Error: path does not exist: {target}", file=sys.stderr)
        raise SystemExit(1)

    directory = target.parent if target.is_file() else target
    default_file = target.name if target.is_file() else None

    handler = partial(ZarrCompatibleServer, directory=str(directory), default_file=default_file)
    print(f"Serving Zarr (Range + CORS) on port {port} from {directory}...")
    httpd = HTTPServer(("", port), handler)
    httpd.serve_forever()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Path to serve")
    parser.add_argument("--port", type=int, default=8000, help="Port")
    args = parser.parse_args()
    run_server(args.path, args.port)