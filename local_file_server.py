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

    def send_range_not_satisfiable(self, file_size):
        # RFC 7233 asks for the current length so the client can retry sensibly.
        self.send_response(416)
        self.send_header("Content-Range", f"bytes */{file_size}")
        self.send_header("Content-Length", "0")
        self.end_headers()

    def do_GET(self):
        if self.default_file and self.path in ("", "/"):
            self.path = f"/{self.default_file}"

        # Get the physical path of the file
        path = self.translate_path(self.path)

        # If it's not a file or no Range header, use default behavior
        if not os.path.isfile(path) or "Range" not in self.headers:
            return super().do_GET()

        # Parse Range header. Three forms are valid (RFC 7233 §2.1), and zarr clients
        # use all of them: sharded arrays store the chunk index at the end of the shard
        # and fetch it with a suffix range, then fetch the chunk itself with an explicit
        # range.
        #   bytes=start-end   explicit range
        #   bytes=start-      from start to EOF
        #   bytes=-suffix     the last `suffix` bytes
        range_match = re.fullmatch(r'bytes=(\d*)-(\d*)', self.headers['Range'].strip())
        # A comma means a multipart range, which needs a multipart/byteranges body we do
        # not build. Serving the whole file is the honest fallback.
        if not range_match or "," in self.headers['Range']:
            return super().do_GET()

        first, last = range_match.groups()
        if not first and not last:
            # "bytes=-" specifies nothing at all.
            return super().do_GET()

        file_size = os.path.getsize(path)

        if not first:
            # Suffix range: the last N bytes. Asking for more than the file holds is
            # allowed and yields the whole file.
            suffix_length = int(last)
            if suffix_length == 0:
                return self.send_range_not_satisfiable(file_size)
            start = max(0, file_size - suffix_length)
            end = file_size - 1
        else:
            start = int(first)
            # An end past EOF is clamped rather than rejected, otherwise Content-Length
            # would promise more bytes than the body carries and the client would hang.
            end = min(int(last), file_size - 1) if last else file_size - 1

        # Validation: ensure range is within file bounds
        if start >= file_size or start > end:
            return self.send_range_not_satisfiable(file_size)

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