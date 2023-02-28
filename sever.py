from http.server import HTTPServer
from http.server import BaseHTTPRequestHandler
import MolDisplay
import sys


class Myserver(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/":
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            html = """
                <html>
                <head>
                <title> File Upload </title>
                </head>
                <body>
                <h1> File Upload </h1>
                <form action="molecule" enctype="multipart/form-data" method="post">
                <p>
                <input type="file" id="sdf_file" name="filename"/>
                </p>
                <p>
                <input type="submit" value="Upload"/>
                </p>
                </form>
                </body>
                </html>
            """
            self.send_header("Content-length", len(html))
            self.end_headers()

            self.wfile.write(html.encode())

        else:
            self.send_response(404)
            self.end_headers()
            self.wfile.write(bytes("404: not found", "utf-8"))

    def do_POST(self):
        if self.path == "/molecule":
            length = int(self.headers.get("content-length"))
            post_data = self.rfile.read(length).decode().split('\n')[4:]
            # print(post_data)

            mol = MolDisplay.Molecule()
            mol.parse(post_data)
            mol.sort()
            svg_data = mol.svg().encode()

            self.send_response(200)
            self.send_header("Content-type", "image/svg+xml")
            self.end_headers()
            self.wfile.write(svg_data)
        else:
            self.send_error(404)


# Student Id for port no
httpd = HTTPServer(('localhost', int(sys.argv[1])), Myserver)
httpd.serve_forever()
