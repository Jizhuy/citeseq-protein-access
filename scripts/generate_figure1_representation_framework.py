"""Generate manuscript Figure 1: CITE-seq representation and comparison framework.

Original conceptual schematic of how the evaluated representations are constructed
(RNA PCA, protein PCA, concat RNA-protein PCA, scVI_matched, totalVI), how they feed a
common downstream classifier, and how the two descriptive contrasts are defined.
Preprocessing and capacity settings follow the manuscript Methods. No data are plotted
and no published architecture figure is reproduced.

Outputs SVG (editable text), PDF, and 600-dpi PNG to figures/source/.

Usage:
    python scripts/generate_figure1_representation_framework.py
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures" / "source"
STEM = "Figure1_representation_framework"

NAVY = "#1B3A5F"
BLUE = "#2F6FAE"
MBLUE = "#9DC3E6"
LBLUE = "#EAF2FB"
PALE = "#F5F8FC"
EDGE = "#B8CCE4"
GRAY = "#5A5A5A"
LGRAY = "#D9DEE7"
WHITE = "#FFFFFF"

MIN_PT = 7.0
W, H = 7.0, 7.08

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

G_ACCESS = r"$G_{\mathrm{access}}$"
G_ASSOC = r"$G_{\mathrm{assoc}}$"

TAG_X, TAG_W = 1.08, 0.64
COL_W = [0.46, 0.46, 0.46, 0.46, 0.66, 0.66]
COL_GAP = 0.10
COL_X0 = 1.82
CHIP_X, CHIP_W = 5.54, 1.18
BUS_X = 6.86


def col(i):
    x = COL_X0 + sum(COL_W[:i]) + COL_GAP * i
    return x, COL_W[i]


def span(i, j):
    x0, _ = col(i)
    xj, wj = col(j)
    return x0, xj + wj - x0


def box(ax, x, y, w, h, fc=LBLUE, ec=EDGE, lw=0.8, r=0.04, z=1, ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, zorder=z, linestyle=ls))


def text(ax, x, y, s, size=7.5, color=NAVY, weight="normal", ha="center", va="center",
         style="normal", z=6, ls=1.1):
    assert size >= MIN_PT, f"font below {MIN_PT} pt: {s!r}"
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va,
            fontstyle=style, zorder=z, linespacing=ls)


def arrow(ax, x1, y1, x2, y2, color=BLUE, lw=0.9, ms=6, z=4, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=ms,
                                 color=color, lw=lw, zorder=z, shrinkA=0, shrinkB=0,
                                 linestyle=ls))


def step(ax, x, y, w, h, label, fc=WHITE, ec=NAVY, tc=NAVY, size=7.0, weight="normal", lw=0.8,
         ls="-"):
    box(ax, x, y, w, h, fc=fc, ec=ec, lw=lw, r=0.03, z=2, ls=ls)
    text(ax, x + w / 2, y + h / 2, label, size=size, color=tc, weight=weight, ls=1.0)


def tag(ax, x, y, w, h, kind):
    spec = {"rna": ("RNA counts", MBLUE, NAVY), "adt": ("ADT counts", NAVY, WHITE),
            "batch": ("batch", LGRAY, NAVY)}[kind]
    box(ax, x, y, w, h, fc=spec[1], ec=spec[1] if kind != "batch" else "#AEB7C6", lw=0.6,
        r=0.03, z=2)
    text(ax, x + w / 2, y + h / 2, spec[0], size=7, color=spec[2], weight="bold")


def letter(ax, x, y, s):
    text(ax, x, y, s, size=10, weight="bold", ha="left")


def grid_icon(ax, x, y, nx, ny, cell, shades):
    for i in range(nx):
        for j in range(ny):
            ax.add_patch(Rectangle((x + i * cell, y + j * cell), cell, cell,
                                   fc=shades[(i * 3 + j * 5) % len(shades)], ec=WHITE, lw=0.4,
                                   zorder=3))


def encoder_shape(ax, x, y, w, h, label, fc=LBLUE, ec=NAVY, tc=NAVY):
    inset = h * 0.22
    ax.add_patch(Polygon([(x, y), (x + w, y + inset), (x + w, y + h - inset), (x, y + h)],
                         closed=True, fc=fc, ec=ec, lw=0.9, zorder=2))
    text(ax, x + w / 2, y + h / 2, label, size=7.2, color=tc, ls=1.05)


def panel_a(ax, x, y, w, h):
    box(ax, x, y, w, h, fc=WHITE, ec=EDGE, lw=0.9, r=0.06, z=0)
    letter(ax, x + 0.07, y + h - 0.16, "A")
    text(ax, x + w / 2 + 0.06, y + h - 0.40, "Paired\nCITE-seq\nmeasurements", size=7.5,
         weight="bold", ls=1.05)
    cx, cy, r = x + w / 2, y + h - 0.98, 0.19
    ax.add_patch(Circle((cx, cy), r, fc=LBLUE, ec=BLUE, lw=1.0, zorder=3))
    ax.add_patch(Circle((cx + 0.02, cy - 0.01), 0.07, fc=MBLUE, ec=BLUE, lw=0.7, zorder=3))
    for k in range(8):
        ang = 2 * math.pi * k / 8 + 0.3
        bx, by = cx + r * math.cos(ang), cy + r * math.sin(ang)
        dx, dy = math.cos(ang), math.sin(ang)
        ex, ey = bx + dx * 0.05, by + dy * 0.05
        ax.plot([bx, ex], [by, ey], color=NAVY, lw=0.8, zorder=4)
        for s in (1, -1):
            ax.plot([ex, ex + (dx - s * dy) * 0.035], [ey, ey + (dy + s * dx) * 0.035],
                    color=NAVY, lw=0.8, zorder=4)
    text(ax, cx, cy - r - 0.14, "one cell", size=7, color=GRAY)

    ry = cy - 0.78
    arrow(ax, cx, cy - r - 0.24, cx, ry + 0.33, lw=0.8, ms=5)
    grid_icon(ax, cx - 0.21, ry, 6, 4, 0.07, [MBLUE, BLUE, LBLUE, "#5C8FC7"])
    text(ax, cx, ry - 0.20, "RNA raw\ncounts", size=7.5, weight="bold", ls=1.0)
    text(ax, cx, ry - 0.47, "(genes)", size=7, color=GRAY)

    py = ry - 1.10
    grid_icon(ax, cx - 0.21, py, 6, 3, 0.07, [NAVY, "#3E5F86", MBLUE, BLUE])
    text(ax, cx, py - 0.20, "ADT raw\ncounts", size=7.5, weight="bold", ls=1.0)
    text(ax, cx, py - 0.49, "(surface\nproteins)", size=7, color=GRAY, ls=1.0)

    by = py - 1.15
    box(ax, x + 0.10, by, w - 0.20, 0.22, fc=LGRAY, ec="#AEB7C6", lw=0.6, r=0.03, z=2)
    text(ax, cx, by + 0.11, "batch", size=7, weight="bold")
    text(ax, cx, by - 0.25, "covariate for\ndeep models", size=7, color=GRAY, ls=1.0)


def band(ax, y, h, title, fc):
    box(ax, 1.00, y, BUS_X + 0.08 - 1.00, h, fc=fc, ec=EDGE, lw=0.8, r=0.06, z=0)
    text(ax, 1.10, y + h - 0.13, title, size=8.5, weight="bold", color=BLUE, ha="left")


def chip(ax, yc, name, sub, kind):
    h = 0.40
    y = yc - h / 2
    if kind == "totalvi":
        fc, ec, tc, sc, lw, ls = NAVY, NAVY, WHITE, LBLUE, 1.0, "-"
    elif kind == "concat":
        fc, ec, tc, sc, lw, ls = MBLUE, NAVY, NAVY, NAVY, 1.4, "-"
    elif kind == "secondary":
        fc, ec, tc, sc, lw, ls = WHITE, "#8FA3BC", GRAY, GRAY, 0.8, "--"
    else:
        fc, ec, tc, sc, lw, ls = WHITE, NAVY, NAVY, GRAY, 1.0, "-"
    box(ax, CHIP_X, y, CHIP_W, h, fc=fc, ec=ec, lw=lw, r=0.04, z=2, ls=ls)
    text(ax, CHIP_X + CHIP_W / 2, yc + 0.075, name, size=7.5, color=tc, weight="bold")
    text(ax, CHIP_X + CHIP_W / 2, yc - 0.095, sub, size=7, color=sc)
    dashed = kind == "secondary"
    ax.plot([CHIP_X + CHIP_W, BUS_X], [yc, yc], color=("#8FA3BC" if dashed else NAVY), lw=0.9,
            ls=("--" if dashed else "-"), zorder=3)


def build():
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    m = 0.05
    lin_h, deep_h = 2.13, 2.65
    strip_h, g_h = 0.60, 1.32
    y_g = m
    y_strip = y_g + g_h + 0.10
    y_deep = y_strip + strip_h + 0.10
    y_lin = y_deep + deep_h + 0.10
    top = y_lin + lin_h

    panel_a(ax, m, y_deep, 0.88, top - y_deep)

    # ---------------- linear band ----------------
    band(ax, y_lin, lin_h, "Linear (PCA) representations", PALE)
    text(ax, CHIP_X + CHIP_W / 2, top - 0.13, "Representation", size=8.5, weight="bold",
         color=BLUE)
    bh = 0.40

    r1_title = top - 0.40
    letter(ax, 0.96 + 0.0, r1_title, "B")
    text(ax, TAG_X + 0.04, r1_title, "RNA PCA branch", size=8, weight="bold", ha="left")
    r1 = r1_title - 0.16 - bh
    tag(ax, TAG_X, r1 + 0.09, TAG_W, 0.22, "rna")
    labels = ["2,000\nHVGs", "Normalize\n10$^4$/cell", "log1p", "z-score\nclip ±10", "PCA\n10 PCs"]
    for i, lab in enumerate(labels):
        x, w = col(i)
        step(ax, x, r1, w, bh, lab, fc=(LBLUE if i == 4 else WHITE),
             weight=("bold" if i == 4 else "normal"))
        arrow(ax, (TAG_X + TAG_W if i == 0 else col(i - 1)[0] + col(i - 1)[1]) + 0.01,
              r1 + bh / 2, x - 0.01, r1 + bh / 2, ms=5)
    x5, w5 = col(4)
    arrow(ax, x5 + w5 + 0.01, r1 + bh / 2, CHIP_X - 0.01, r1 + bh / 2, ms=6)
    chip(ax, r1 + bh / 2, "RNA PCA", "10-D · RNA only", "plain")

    r2_title = r1 - 0.20
    letter(ax, 0.96, r2_title, "C")
    text(ax, TAG_X + 0.04, r2_title, "Protein PCA and concat RNA–protein PCA", size=8,
         weight="bold", ha="left")
    r2 = r2_title - 0.14 - bh
    xc, wc = col(4)
    step(ax, xc, r2, wc, bh, "Concatenate\n10 + 10 PCs", fc=MBLUE, weight="bold")
    xs, ws = col(5)
    step(ax, xs, r2, ws, bh, "Standard-\nScaler", fc=MBLUE, weight="bold")
    arrow(ax, xc + wc + 0.01, r2 + bh / 2, xs - 0.01, r2 + bh / 2, ms=5)
    arrow(ax, xs + ws + 0.01, r2 + bh / 2, CHIP_X - 0.01, r2 + bh / 2, ms=6)
    arrow(ax, xc + wc / 2, r1 - 0.01, xc + wc / 2, r2 + bh + 0.01, ms=5)
    lx, lw_ = span(0, 3)
    box(ax, TAG_X, r2 + 0.02, lx + lw_ - TAG_X, bh - 0.04, fc=WHITE, ec=MBLUE, lw=0.8, r=0.03,
        z=1)
    text(ax, (TAG_X + lx + lw_) / 2, r2 + bh / 2,
         "Simple multimodal feature-access reference:\nsame RNA + protein access as totalVI, linear",
         size=7.2, color=NAVY, ls=1.1)
    chip(ax, r2 + bh / 2, "Concat PCA", "20-D · RNA + protein", "concat")

    r3 = r2 - 0.12 - bh
    tag(ax, TAG_X, r3 + 0.09, TAG_W, 0.22, "adt")
    x0, w01 = span(0, 1)
    step(ax, x0, r3, w01, bh, "CLR\npseudocount 1")
    x2, w23 = span(2, 3)
    step(ax, x2, r3, w23, bh, "z-score\n(train mean, SD)")
    xp, wp = col(4)
    step(ax, xp, r3, wp, bh, "PCA\n10 PCs", fc=LBLUE, weight="bold")
    arrow(ax, TAG_X + TAG_W + 0.01, r3 + bh / 2, x0 - 0.01, r3 + bh / 2, ms=5)
    arrow(ax, x0 + w01 + 0.01, r3 + bh / 2, x2 - 0.01, r3 + bh / 2, ms=5)
    arrow(ax, x2 + w23 + 0.01, r3 + bh / 2, xp - 0.01, r3 + bh / 2, ms=5)
    arrow(ax, xp + wp / 2, r3 + bh + 0.01, xp + wp / 2, r2 - 0.01, ms=5)
    arrow(ax, xp + wp + 0.01, r3 + bh / 2, CHIP_X - 0.01, r3 + bh / 2, ms=6, color="#8FA3BC",
          ls="--")
    chip(ax, r3 + bh / 2, "Protein PCA", "10-D · secondary", "secondary")

    # ---------------- deep band ----------------
    band(ax, y_deep, deep_h, "Deep generative latent representations (established scvi-tools models)",
         LBLUE)
    dtop = y_deep + deep_h

    r4_title = dtop - 0.40
    letter(ax, 0.96, r4_title, "D")
    text(ax, TAG_X + 0.04, r4_title, "scVI_matched — capacity-matched RNA-only scVI [2]",
         size=8, weight="bold", ha="left")
    mh = 0.44
    r4 = r4_title - 0.16 - mh
    tag(ax, TAG_X, r4 + mh / 2 + 0.02, TAG_W, 0.2, "rna")
    tag(ax, TAG_X, r4 + mh / 2 - 0.22, TAG_W, 0.2, "batch")
    ex, ew = span(0, 1)
    encoder_shape(ax, ex, r4, ew, mh, "Encoder\n2 layers × 256", fc=WHITE)
    px, pw = span(2, 3)
    step(ax, px, r4 + 0.05, pw, mh - 0.10, r"$q_\phi(z \mid x)$" + "\nlatent posterior",
         fc=WHITE)
    zx, zw = col(4)
    step(ax, zx, r4 + 0.05, zw, mh - 0.10, r"$z \in \mathbb{R}^{20}$", fc=LBLUE, weight="bold",
         size=7.5)
    arrow(ax, TAG_X + TAG_W + 0.01, r4 + mh / 2, ex - 0.01, r4 + mh / 2, ms=5)
    arrow(ax, ex + ew + 0.01, r4 + mh / 2, px - 0.01, r4 + mh / 2, ms=5)
    arrow(ax, px + pw + 0.01, r4 + mh / 2, zx - 0.01, r4 + mh / 2, ms=5)
    arrow(ax, zx + zw + 0.01, r4 + mh / 2, CHIP_X - 0.01, r4 + mh / 2, ms=6)
    chip(ax, r4 + mh / 2, "scVI_matched", "20-D · RNA only", "plain")
    dy4 = r4 - 0.30
    dx, dw = span(4, 5)
    step(ax, dx, dy4, dw, 0.22, "RNA decoder → count likelihood", fc=WHITE, ec="#8FA3BC",
         tc=GRAY, size=7, ls="--", lw=0.7)
    arrow(ax, zx + zw / 2, r4 + 0.04, zx + zw / 2, dy4 + 0.23, color="#8FA3BC", ls="--", ms=5,
          lw=0.8)
    text(ax, dx - 0.06, dy4 + 0.11, "generative (training)", size=7, color=GRAY, style="italic",
         ha="right")

    r5_title = dy4 - 0.20
    letter(ax, 0.96, r5_title, "E")
    text(ax, TAG_X + 0.04, r5_title,
         "totalVI — established multimodal deep generative model [3]", size=8, weight="bold",
         ha="left")
    mh5 = 0.58
    r5 = r5_title - 0.16 - mh5
    for k, kind in enumerate(["rna", "adt", "batch"]):
        tag(ax, TAG_X, r5 + mh5 - 0.19 - k * 0.20, TAG_W, 0.18, kind)
    encoder_shape(ax, ex, r5, ew, mh5, "Joint encoder\n2 layers × 256", fc=WHITE)
    step(ax, px, r5 + 0.08, pw, mh5 - 0.16,
         r"$q_\phi(z \mid x_{\mathrm{RNA}}, x_{\mathrm{ADT}})$" + "\nshared posterior", fc=WHITE)
    step(ax, zx, r5 + 0.08, zw, mh5 - 0.16, "shared\n" + r"$z \in \mathbb{R}^{20}$", fc=NAVY,
         ec=NAVY, tc=WHITE, weight="bold", size=7.5)
    arrow(ax, TAG_X + TAG_W + 0.01, r5 + mh5 / 2, ex - 0.01, r5 + mh5 / 2, ms=5)
    arrow(ax, ex + ew + 0.01, r5 + mh5 / 2, px - 0.01, r5 + mh5 / 2, ms=5)
    arrow(ax, px + pw + 0.01, r5 + mh5 / 2, zx - 0.01, r5 + mh5 / 2, ms=5)
    arrow(ax, zx + zw + 0.01, r5 + mh5 / 2, CHIP_X - 0.01, r5 + mh5 / 2, ms=6)
    chip(ax, r5 + mh5 / 2, "totalVI", "20-D · RNA + protein", "totalvi")
    dy5 = r5 - 0.30
    gx1, gw1 = span(2, 3)
    step(ax, gx1, dy5, gw1, 0.22, "RNA decoder", fc=WHITE, ec="#8FA3BC", tc=GRAY,
         size=7, ls="--", lw=0.7)
    step(ax, dx, dy5, dw, 0.22, "Protein decoder", fc=WHITE, ec="#8FA3BC", tc=GRAY,
         size=7, ls="--", lw=0.7)
    arrow(ax, zx + zw * 0.35, r5 + 0.07, gx1 + gw1 * 0.75, dy5 + 0.23, color="#8FA3BC", ls="--",
          ms=5, lw=0.8)
    text(ax, gx1 - 0.06, dy5 + 0.11, "generative (training)", size=7, color=GRAY, style="italic",
         ha="right")
    arrow(ax, zx + zw * 0.65, r5 + 0.07, zx + zw * 0.65, dy5 + 0.23, color="#8FA3BC", ls="--",
          ms=5, lw=0.8)

    # ---------------- common classifier (F) ----------------
    fx, fw = CHIP_X - 0.04, W - m - (CHIP_X - 0.04)
    fy, fh = m, y_strip + strip_h - m
    box(ax, fx, fy, fw, fh, fc=LBLUE, ec=BLUE, lw=1.1, r=0.06, z=1)
    ax.plot([BUS_X, BUS_X], [r1 + bh / 2, fy + fh + 0.10], color=NAVY, lw=1.1, zorder=3)
    arrow(ax, BUS_X, fy + fh + 0.11, BUS_X, fy + fh - 0.02, color=NAVY, lw=1.1, ms=8)
    letter(ax, fx + 0.07, fy + fh - 0.17, "F")
    text(ax, fx + 0.30, fy + fh - 0.36, "Same downstream\nclassifier for all\nrepresentations",
         size=7.5, weight="bold", ha="left", ls=1.05)
    lx0, lw0 = fx + 0.10, fw - 0.20
    ly = fy + fh - 0.98
    step(ax, lx0, ly, lw0, 0.30, "Logistic regression", fc=NAVY, ec=NAVY, tc=WHITE,
         weight="bold", size=7.5)
    py = ly - 0.40
    step(ax, lx0, py, lw0, 0.26, "Cell-type prediction", fc=WHITE)
    my = py - 0.38
    step(ax, lx0, my, lw0, 0.26, "Macro-F1", fc=WHITE, weight="bold")
    arrow(ax, lx0 + lw0 / 2, ly - 0.01, lx0 + lw0 / 2, py + 0.27, color=NAVY, ms=6)
    arrow(ax, lx0 + lw0 / 2, py - 0.01, lx0 + lw0 / 2, my + 0.27, color=NAVY, ms=6)

    # ---------------- strict train-only strip ----------------
    sx, sw = m, CHIP_X - 0.14 - m
    box(ax, sx, y_strip, sw, strip_h, fc=WHITE, ec=NAVY, lw=0.9, r=0.05, z=1)
    text(ax, sx + 0.10, y_strip + strip_h - 0.12,
         "Strict train-only information access (primary internal analysis)", size=7.5,
         weight="bold", ha="left")
    text(ax, sx + 0.10, y_strip + strip_h - 0.28,
         "All unsupervised transforms and representation models are fitted on training cells only.",
         size=7, ha="left")
    flow = ["held-out cells", "frozen transform / encoder", "representation",
            "classifier prediction"]
    widths = [0.90, 1.45, 1.00, 1.22]
    fx0 = sx + 0.10
    fyc = y_strip + 0.13
    for i, (lab, w) in enumerate(zip(flow, widths)):
        step(ax, fx0, fyc - 0.09, w, 0.18, lab, fc=LBLUE if i else LGRAY, ec=EDGE, size=7)
        if i < len(flow) - 1:
            arrow(ax, fx0 + w + 0.02, fyc, fx0 + w + 0.16, fyc, ms=5)
        fx0 += w + 0.18

    # ---------------- contrasts (G) ----------------
    gx, gw = m, sw
    box(ax, gx, y_g, gw, g_h, fc=WHITE, ec=EDGE, lw=0.9, r=0.06, z=0)
    letter(ax, gx + 0.07, y_g + g_h - 0.16, "G")
    text(ax, gx + 0.30, y_g + g_h - 0.16, "Primary descriptive contrasts", size=8.5,
         weight="bold", ha="left")
    text(ax, gx + gw - 0.10, y_g + g_h - 0.16,
         "Descriptive contrasts; not causal architecture effects.", size=7, color=GRAY,
         style="italic", ha="right")
    rows = [
        (("RNA PCA", WHITE, NAVY), ("Concat", MBLUE, NAVY), "Protein-access gain",
         G_ACCESS + r" = F1$_{\mathrm{concat}}$ − F1$_{\mathrm{RNA\ PCA}}$",
         "adds paired protein to a simple linear representation"),
        (("Concat", MBLUE, NAVY), ("totalVI", NAVY, WHITE), "Residual totalVI-associated contrast",
         G_ASSOC + r" = F1$_{\mathrm{totalVI}}$ − F1$_{\mathrm{concat}}$",
         "same RNA + protein access, different construction"),
    ]
    ry = y_g + g_h - 0.36
    for (a, b, lab, eq, sub) in rows:
        cy0 = ry - 0.40
        for k, (nm, fc, tc) in enumerate((a, b)):
            bx = gx + 0.15 + k * 0.92
            box(ax, bx, cy0 + 0.06, 0.72, 0.28, fc=fc, ec=NAVY, lw=0.8, r=0.03, z=2)
            text(ax, bx + 0.36, cy0 + 0.20, nm, size=7.5, color=tc, weight="bold")
        arrow(ax, gx + 0.15 + 0.73, cy0 + 0.20, gx + 0.15 + 0.91, cy0 + 0.20, ms=5)
        tx = gx + 2.05
        text(ax, tx, cy0 + 0.35, lab, size=7.5, weight="bold", color=BLUE, ha="left")
        text(ax, tx, cy0 + 0.19, eq, size=7.5, ha="left")
        text(ax, tx, cy0 + 0.03, sub, size=7, color=GRAY, ha="left", style="italic")
        ry -= 0.46
    return fig


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = build()
    for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 600})):
        out = OUT_DIR / f"{STEM}.{ext}"
        fig.savefig(out, facecolor="white", **kw)
        print(f"wrote {out.relative_to(ROOT)}")
    plt.close(fig)


if __name__ == "__main__":
    main()
