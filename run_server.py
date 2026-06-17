import http.server
import socketserver
import webbrowser
import threading
import time
import sys

PORT = 8000

class NoCacheHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')
        self.send_header('Pragma', 'no-cache')
        self.send_header('Expires', '0')
        super().end_headers()

def open_browser():
    time.sleep(1.2)
    url = f"http://localhost:{PORT}"
    print(f"-> Browser opening: {url}")
    webbrowser.open(url)

if __name__ == '__main__':
    threading.Thread(target=open_browser, daemon=True).start()
    
    Handler = NoCacheHTTPRequestHandler
    
    socketserver.TCPServer.allow_reuse_address = True
    
    try:
        with socketserver.TCPServer(("", PORT), Handler) as httpd:
            print(f"-> HTTP Server started. Port: {PORT}")
            print("-> Press Ctrl+C to stop.")
            httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n-> Server stopped.")
        sys.exit(0)
    except Exception as e:
        print(f"\n-> Server could not be started: {e}")
        sys.exit(1)
