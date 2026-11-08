"""Shared style of the schematic figures (pipeline_overview, protocol_flow, cgc_rules).

Journal-figure style: white background, near-black hairlines, square corners, one sans-serif family and
three text sizes. Colour only where it encodes something: the CGC gene types in cgc_rules (the palette of
the run_dbcan plots) and, if USE_ACCENT is True, one muted accent for the step numbers and flow arrows.

The coordinates of every figure are inches, origin top left, y pointing down.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, Rectangle

OUT_DIR = Path(__file__).resolve().parents[1]
DPI = 150

# Variant B of pipeline_overview (monochrome + one muted accent) is the docs default; set to False for the
# strictly monochrome variant A in all three figures.
USE_ACCENT = True


def _available(candidates):
    have = {f.name for f in font_manager.fontManager.ttflist}
    return [c for c in candidates if c in have] or [candidates[-1]]


SANS = _available(["Inter", "Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"])
MONO = _available(["JetBrains Mono", "Liberation Mono", "DejaVu Sans Mono"])
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": SANS,
    "font.monospace": MONO,
    "savefig.facecolor": "white",
})

# text sizes (pt): titles / step names, body, notes and file names
L, M, S = 10, 9, 8

INK = "#1A1A1A"        # text and structure lines
GREY = "#5C5C5C"       # secondary text, file names (contrast 6.7:1 on white)
RULE = "#BFBFBF"       # light separators, gene-row backbone
ACCENT_COLOUR = "#2F5D86"   # muted steel blue
LW = 0.8               # hairline width (pt)


def accent(use=None):
    """The accent colour, or INK when the figure is monochrome."""
    return ACCENT_COLOUR if (USE_ACCENT if use is None else use) else INK


def canvas(width, height):
    fig = plt.figure(figsize=(width, height), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")
    return fig, ax


def text(ax, x, y, s, size=M, color=INK, weight="normal", mono=False, **kw):
    kw.setdefault("va", "center")
    kw.setdefault("ha", "left")
    return ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
                   family="monospace" if mono else "sans-serif", **kw)


def box(ax, x, y, w, h, ec=INK, fc="white", lw=LW, **kw):
    ax.add_patch(Rectangle((x, y), w, h, ec=ec, fc=fc, lw=lw, joinstyle="miter", **kw))


def hline(ax, x0, x1, y, color=RULE, lw=LW, **kw):
    ax.plot([x0, x1], [y, y], color=color, lw=lw, solid_capstyle="butt", **kw)


def vline(ax, x, y0, y1, color=INK, lw=LW, **kw):
    ax.plot([x, x], [y0, y1], color=color, lw=lw, solid_capstyle="butt", **kw)


def arrow(ax, xy0, xy1, color=INK, lw=LW, head=7, **kw):
    ax.add_patch(FancyArrowPatch(xy0, xy1, arrowstyle="-|>", mutation_scale=head, color=color, lw=lw,
                                 shrinkA=0, shrinkB=0, joinstyle="miter", capstyle="butt", **kw))


def save(fig, name):
    out = OUT_DIR / name
    fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    try:  # smaller file: drop the alpha channel (the background is opaque white)
        from PIL import Image
        Image.open(out).convert("RGB").save(out, optimize=True)
    except ImportError:
        pass
    print(out, out.stat().st_size, "bytes")
    return out


def runs(ax, x, y, parts, **common):
    """Draw text pieces side by side on one line: parts = [(string, {text() keywords}), ...].
    Returns the x after the last piece (inches)."""
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    for s, kw in parts:
        t = text(ax, x, y, s, **{**common, **kw})
        x += t.get_window_extent(renderer).width / fig.dpi
    return x
