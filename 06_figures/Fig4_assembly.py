"""
Fig4_assembly.py — Combine Fig 4: woodland (a) + openhabitat (b) + wetland (c) + boxplot (d)
Vertical stack, 180mm width, NEE submission spec.
"""

import fitz  # PyMuPDF
from pathlib import Path
import os, sys

FIG_DIR = Path(r"D:\OneDrive - University of Leeds\biodiversity_manscript\figures")

panel_files = {
    "a": FIG_DIR / "Fig4a_woodland_180mm_revised_cairo.pdf",
    "b": FIG_DIR / "Fig4b_openhabitat_180mm_revised_cairo.pdf",
    "c": FIG_DIR / "Fig4c_wetland_180mm_revised_cairo.pdf",
    "d": FIG_DIR / "Fig4d_guild_boxplot_180mm.pdf",
}

OUTPUT = FIG_DIR / "Fig4_combined_180mm.pdf"

MM_TO_PT = 2.834
WIDTH_MM = 180
WIDTH_PT = WIDTH_MM * MM_TO_PT
GAP_MM = 3
GAP_PT = GAP_MM * MM_TO_PT
LABEL_FONT = "hebo"    # Helvetica Bold
LABEL_SIZE = 10
LABEL_X_OFFSET = 10
LABEL_Y_OFFSET = 10

# Check files + read dimensions
print("=" * 70)
print("  Fig 4 Assembly — 180mm x variable height")
print("=" * 70)

panel_dims = {}
for label, path in panel_files.items():
    if not path.exists():
        raise FileNotFoundError(f"  Missing: {path}")
    doc = fitz.open(path)
    w_pt, h_pt = doc[0].rect.width, doc[0].rect.height
    panel_dims[label] = (w_pt, h_pt)
    doc.close()
    print(f"  Panel ({label}): {w_pt/MM_TO_PT:6.2f} x {h_pt/MM_TO_PT:6.2f} mm  [{path.name}]")

total_h_pt = sum(h for _, h in panel_dims.values()) + GAP_PT * (len(panel_files) - 1)
print(f"\n  Total: {WIDTH_MM} x {total_h_pt/MM_TO_PT:.2f} mm | Gap: {GAP_MM} mm")

# Create merged PDF
out_doc = fitz.open()
out_page = out_doc.new_page(width=WIDTH_PT, height=total_h_pt)

y_cursor = 0
for label in ["a", "b", "c", "d"]:
    src = fitz.open(panel_files[label])
    src_w, src_h = panel_dims[label]
    scale = WIDTH_PT / src_w
    target_h = src_h * scale
    target_rect = fitz.Rect(0, y_cursor, WIDTH_PT, y_cursor + target_h)
    out_page.show_pdf_page(target_rect, src, 0)
    src.close()

    out_page.insert_text(
        fitz.Point(LABEL_X_OFFSET, y_cursor + LABEL_Y_OFFSET),
        label, fontname=LABEL_FONT, fontsize=LABEL_SIZE, color=(0, 0, 0))

    print(f"  Panel ({label}) at y = {y_cursor/MM_TO_PT:6.2f} mm (h = {target_h/MM_TO_PT:.2f} mm)")
    y_cursor += target_h + GAP_PT

PREVIEW = FIG_DIR / "Fig4_combined_180mm_PREVIEW.pdf"
out_doc.save(PREVIEW, garbage=4, deflate=True)
out_doc.close()
print(f"\n  Saved: {PREVIEW.name}")
print(f"  Final: {WIDTH_MM} x {(y_cursor - GAP_PT)/MM_TO_PT:.2f} mm")

if sys.platform == "win32":
    os.startfile(str(PREVIEW))
