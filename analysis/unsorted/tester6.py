import os
import math
import xml.etree.ElementTree as ET
from svgutils.compose import Figure, SVG, Text

# Configuration
input_dir = "figure4_svg"
output_svg = os.path.join(input_dir, "grid.svg")

def get_svg_dimensions(filepath):
    try:
        tree = ET.parse(filepath)
        root = tree.getroot()
        width = float(root.attrib.get("width", "1").replace("px", ""))
        height = float(root.attrib.get("height", "1").replace("px", ""))
        return width, height
    except Exception as e:
        print(f"[!] Failed to parse dimensions from {filepath}: {e}")
        return 1.0, 1.0

# Grid layout
cols = 7
cell_width = 70
cell_height = 95
label_height = 10
font_size = 9

# Get SVG files (exclude previous grid if present)
svg_files = sorted([
    f for f in os.listdir(input_dir)
    if f.endswith(".svg") and not f.startswith("grid")
])

# Compute grid size
rows = math.ceil(len(svg_files) / cols)
width = cols * cell_width
height = rows * cell_height

# Compose the grid
elements = []

for i, filename in enumerate(svg_files):
    col = i % cols
    row = i // cols
    x = col * cell_width
    y = row * cell_height

    svg_path = os.path.join(input_dir, filename)
    # Extract EVODEX ID from filename (remove .svg extension)
    evodex_id = os.path.splitext(filename)[0]
    evodex_id = evodex_id.replace("EVODEX.1-", "")

    orig_width, orig_height = get_svg_dimensions(svg_path)
    scale_w = (cell_width - 5) / orig_width
    scale_h = (cell_height - label_height - 5) / orig_height
    scale_factor = min(scale_w, scale_h)

    try:
        scaled_height = orig_height * scale_factor
        center_offset = ((cell_height - label_height) - scaled_height) / 2
        elements.append(
            SVG(svg_path).scale(scale_factor).move(x, y + label_height + center_offset + -5)  
        )
        elements.append(
            Text(evodex_id, x + cell_width / 2, y + label_height - 5, size=font_size, font="Arial", anchor="middle")
        )
    except Exception as e:
        print(f"[!] Failed to inline SVG {filename}: {e}")

# Save the composed SVG
Figure(width, height, *elements).save(output_svg)
print(f"[✓] Composite grid saved as {output_svg}")