import os
import base64

def image_to_datauri(path):
    ext = os.path.splitext(path)[1].lower().lstrip('.')
    mime = 'image/png' if ext=='png' else ('image/svg+xml' if ext=='svg' else 'application/octet-stream')
    with open(path, 'rb') as fh:
        data = fh.read()
    return f"data:{mime};base64,{base64.b64encode(data).decode('ascii')}"
