from http.server import HTTPServer
from http.server import BaseHTTPRequestHandler
import sys
import MolDisplay


class Myserver(BaseHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/":
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            self.send_header("Content-length", len(home_page))
            self.end_headers()
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
            self.wfile.write(html.encode())

        else:
            self.send_response(404)
            self.end_headers()
            self.wfile.write(bytes("404: not found", "utf-8"))

    def do_POST(self):
        if self.path == "/molecule":
            # skip 4 lines of header information
            
            length = int(self.headers.get("content-length"))
            # for i in range(4):
            #     print(self.rfile.readline())
            post_data = self.rfile.read(length).decode().split('\n')[4:]
            # print(post_data)
            # Skip 4 lines of header information
            # lines = post_data.split(b"\n")[4:]
            # sdf_data = b"\n".join(lines)
            # print(sdf_data)
            mol = MolDisplay.Molecule()
            mol.parse(post_data)
            # Call sort() method on the molecule
            mol.sort()
            svg_data = mol.svg().encode()
            self.send_response(200)
            self.send_header("Content-type", "image/svg+xml")
            self.end_headers()
            self.wfile.write(svg_data)
        else:
            self.send_error(404)


home_page = """
<html>
  <head>
    <title> Hello, world! </title>
  </head>
  <body>
    <h1> Hello, world!\n </h1>
    <p> Welcome to the world of "hello"! </p>
  </body>
</html>
"""

httpd = HTTPServer(('localhost', int(sys.argv[1])), Myserver)
httpd.serve_forever()
