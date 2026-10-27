"""Generate Supplementary Figure S1: study-design overview schematic (panels A-G).

Original vector schematic (no data plotted, no third-party artwork). Outputs SVG, PDF, and 600-dpi PNG
to figures/Supplementary_Figure_S1.{svg,pdf,png}.

Usage:
    python scripts/generate_supplementary_figure_s1.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
STEM = "Supplementary_Figure_S1"

NAVY = "#1B3A5F"
BLUE = "#2F6FAE"
MBLUE = "#9DC3E6"
LBLUE = "#EAF2FB"
PALE = "#F5F8FC"
EDGE = "#B8CCE4"
GRAY = "#5A5A5A"
LGRAY = "#D9DEE7"
WHITE = "#FFFFFF"

W, H = 7.0, 9.45

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


def box(ax, x, y, w, h, fc=LBLUE, ec=EDGE, lw=0.8, r=0.05, z=1, ls="-"):
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h,
            boxstyle=f"round,pad=0,rounding_size={r}",
            fc=fc, ec=ec, lw=lw, zorder=z, linestyle=ls,
        )
    )


def text(ax, x, y, s, size=6.8, color=NAVY, weight="normal", ha="center", va="center",
         style="normal", z=5, ls=1.15):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va,
            fontstyle=style, zorder=z, linespacing=ls)


def arrow(ax, x1, y1, x2, y2, color=BLUE, lw=1.0, ms=7, z=4, style="-|>"):
    ax.add_patch(
        FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=ms,
                        color=color, lw=lw, zorder=z, shrinkA=0, shrinkB=0)
    )


def panel(ax, x, y, w, h, letter, title):
    box(ax, x, y, w, h, fc=WHITE, ec=EDGE, lw=0.9, r=0.07, z=0)
    text(ax, x + 0.10, y + h - 0.15, letter, size=10, weight="bold", ha="left")
    text(ax, x + 0.32, y + h - 0.15, title, size=8, weight="bold", ha="left")


def bracket(ax, x1, x2, y, depth=0.07, color=NAVY, lw=0.9):
    ax.plot([x1, x1, x2, x2], [y + depth, y, y, y + depth], color=color, lw=lw,
            solid_capstyle="butt", zorder=3)


def grid_icon(ax, x, y, nx, ny, cell, shades):
    for i in range(nx):
        for j in range(ny):
            c = shades[(i * 3 + j * 5) % len(shades)]
            ax.add_patch(Rectangle((x + i * cell, y + j * cell), cell, cell, fc=c,
                                   ec=WHITE, lw=0.4, zorder=3))


def antibody(ax, x, y, ang, size=0.07, color=NAVY):
    import math

    dx, dy = math.cos(ang), math.sin(ang)
    px, py = -dy, dx
    bx, by = x + dx * size, y + dy * size
    ax.plot([x, bx], [y, by], color=color, lw=0.9, zorder=4)
    for s in (1, -1):
        ax.plot([bx, bx + (dx + s * px) * size * 0.6], [by, by + (dy + s * py) * size * 0.6],
                color=color, lw=0.9, zorder=4)


def panel_a(ax, x, y, w, h):
    panel(ax, x, y, w, h, "A", "Paired CITE-seq input")
    import math

    cx, cy, r = x + 0.48, y + h - 0.95, 0.27
    ax.add_patch(Circle((cx, cy), r, fc=LBLUE, ec=BLUE, lw=1.0, zorder=3))
    ax.add_patch(Circle((cx + 0.02, cy - 0.02), 0.10, fc=MBLUE, ec=BLUE, lw=0.7, zorder=3))
    for k in range(8):
        ang = 2 * math.pi * k / 8 + 0.2
        antibody(ax, cx + r * math.cos(ang), cy + r * math.sin(ang), ang)
    text(ax, cx, cy - r - 0.17, "single cell", size=6.3, color=GRAY)

    mx = x + 1.02
    arrow(ax, cx + r + 0.08, cy + 0.12, mx - 0.04, y + h - 0.62, lw=0.8, ms=6)
    arrow(ax, cx + r + 0.08, cy - 0.12, mx - 0.04, y + h - 1.22, lw=0.8, ms=6)
    grid_icon(ax, mx, y + h - 0.78, 5, 3, 0.075, [MBLUE, BLUE, LBLUE, "#5C8FC7"])
    text(ax, mx + 0.44, y + h - 0.66, "RNA\ncounts", size=6.5, ha="left", weight="bold", ls=1.0)
    text(ax, mx, y + h - 0.88, "transcriptome", size=6.0, color=GRAY, ha="left")
    grid_icon(ax, mx, y + h - 1.38, 5, 3, 0.075, [NAVY, "#3E5F86", MBLUE, BLUE])
    text(ax, mx + 0.44, y + h - 1.26, "ADT\ncounts", size=6.5, ha="left", weight="bold", ls=1.0)
    text(ax, mx, y + h - 1.56, "surface proteins\n(antibody-derived tags)", size=6.0, color=GRAY,
         ha="left", ls=1.05)

    cy2 = y + 0.42
    for i, lab in enumerate(["PBMC10k", "PBMC5k"]):
        bx = x + 0.14 + i * 0.97
        box(ax, bx, cy2 - 0.12, 0.90, 0.24, fc=NAVY, ec=NAVY, r=0.04, z=2)
        text(ax, bx + 0.45, cy2, lab, size=7, color=WHITE, weight="bold")
    text(ax, x + w / 2, y + 0.15, "Paired RNA and protein measured\nin the same cells",
         size=6.0, color=GRAY, ls=1.05)


def panel_b(ax, x, y, w, h):
    panel(ax, x, y, w, h, "B", "Representation framework")
    top = y + h - 0.34
    gh = 1.20
    gw = (w - 0.36) / 2
    groups = [
        ("RNA only", [("RNA PCA", "RNA-only\nlinear\nreference"),
                      ("scVI_matched", "capacity-\nmatched\nRNA-only VAE")]),
        ("RNA + protein (multimodal)", [("Concat\nRNA–protein PCA", "simple multimodal\nfeature-access\nreference"),
                                        ("totalVI", "existing\nmultimodal model;\nGayoso et al.")]),
    ]
    centers = []
    for gi, (gname, items) in enumerate(groups):
        gx = x + 0.12 + gi * (gw + 0.12)
        gy = top - gh
        box(ax, gx, gy, gw, gh, fc=PALE, ec=EDGE, r=0.05, z=1)
        text(ax, gx + gw / 2, top - 0.13, gname, size=7.2, weight="bold", color=BLUE)
        iw = (gw - 0.30) / 2
        for ii, (name, sub) in enumerate(items):
            ix = gx + 0.10 + ii * (iw + 0.10)
            iy = gy + 0.10
            ih = gh - 0.36
            deep = name in ("scVI_matched", "totalVI")
            fc = NAVY if name == "totalVI" else (MBLUE if gi == 1 else LGRAY)
            if name == "scVI_matched":
                fc = "#C5D3E3"
            box(ax, ix, iy, iw, ih, fc=fc, ec=NAVY, lw=0.8, r=0.04, z=2)
            tc = WHITE if name == "totalVI" else NAVY
            text(ax, ix + iw / 2, iy + ih - 0.18, name, size=6.9, weight="bold", color=tc,
                 ls=1.0)
            text(ax, ix + iw / 2, iy + 0.33, sub, size=6.0,
                 color=(LBLUE if name == "totalVI" else GRAY), ls=1.05)
            text(ax, ix + iw / 2, iy + 0.08, "deep generative" if deep else "linear / PCA",
                 size=5.8, color=(MBLUE if name == "totalVI" else BLUE), style="italic")
            centers.append(ix + iw / 2)

    cly = y + 0.16
    clh = 0.30
    box(ax, x + 0.55, cly, w - 1.10, clh, fc=LBLUE, ec=BLUE, lw=0.9, r=0.05, z=2)
    text(ax, x + w / 2, cly + clh / 2,
         "Common downstream classifier: logistic regression on each representation → cell-type labels",
         size=6.5, weight="bold")
    for cxx in centers:
        arrow(ax, cxx, top - gh - 0.01, cxx, cly + clh + 0.01, color=BLUE, lw=0.8, ms=5)


def panel_c(ax, x, y, w, h):
    panel(ax, x, y, w, h, "C", "Comparison logic")
    bw, bh = 1.75, 0.46
    by = y + h - 0.36 - bh
    xs = [x + 0.25, x + w / 2 - bw / 2, x + w - 0.25 - bw]
    spec = [
        ("RNA PCA", "RNA-only reference", LGRAY, NAVY),
        ("Concat RNA–protein PCA", "simple multimodal reference", MBLUE, NAVY),
        ("totalVI", "complex multimodal representation", NAVY, WHITE),
    ]
    for bx, (name, sub, fc, tc) in zip(xs, spec):
        box(ax, bx, by, bw, bh, fc=fc, ec=NAVY, lw=0.9, r=0.05, z=2)
        text(ax, bx + bw / 2, by + bh - 0.15, name, size=7.2, weight="bold", color=tc)
        text(ax, bx + bw / 2, by + 0.13, sub, size=6.2, color=(LBLUE if tc == WHITE else GRAY))
    for i in range(2):
        arrow(ax, xs[i] + bw + 0.05, by + bh / 2, xs[i + 1] - 0.05, by + bh / 2,
              color=BLUE, lw=1.2, ms=8)

    c = [bx + bw / 2 for bx in xs]
    yb = by - 0.08
    bracket(ax, c[0], c[1] - 0.04, yb, depth=0.06)
    bracket(ax, c[1] + 0.04, c[2], yb, depth=0.06)
    text(ax, (c[0] + c[1]) / 2, yb - 0.12, "Protein-access gain", size=7, weight="bold", color=BLUE)
    text(ax, (c[0] + c[1]) / 2, yb - 0.29,
         G_ACCESS + r" = F1$_{\mathrm{concat}}$ − F1$_{\mathrm{RNA\ PCA}}$", size=7)
    text(ax, (c[1] + c[2]) / 2, yb - 0.12, "Residual totalVI-associated contrast", size=7,
         weight="bold", color=BLUE)
    text(ax, (c[1] + c[2]) / 2, yb - 0.29,
         G_ASSOC + r" = F1$_{\mathrm{totalVI}}$ − F1$_{\mathrm{concat}}$", size=7)

    qy = y + 0.09
    box(ax, x + 0.25, qy, w - 0.50, 0.36, fc=LBLUE, ec=EDGE, r=0.04, z=1)
    text(ax, x + w / 2, qy + 0.245,
         "Central question: how much of the observed multimodal gain is associated with protein access, and does a",
         size=6.5, weight="bold")
    text(ax, x + w / 2, qy + 0.11,
         "residual totalVI-associated contrast remain beyond a simple RNA–protein representation?",
         size=6.5, weight="bold")
    text(ax, x + w - 0.12, y + h - 0.15, "Descriptive contrasts; not causal architecture effects.",
         size=6.3, color=GRAY, style="italic", ha="right")


def panel_d(ax, x, y, w, h):
    panel(ax, x, y, w, h, "D", "Internal evaluation protocol")
    bx, bw = x + 0.15, w - 0.30
    ty = y + h - 0.62
    text(ax, x + w / 2, ty + 0.30, "9,494 high-confidence PBMC cells (development cohort)",
         size=6.4, color=GRAY)
    frac = 7573 / 9494
    box(ax, bx, ty, bw * frac - 0.02, 0.22, fc=BLUE, ec=NAVY, lw=0.7, r=0.03, z=2)
    box(ax, bx + bw * frac + 0.02, ty, bw * (1 - frac) - 0.02, 0.22, fc=LGRAY, ec=NAVY, lw=0.7,
        r=0.03, z=2)
    text(ax, bx + bw * frac / 2, ty + 0.11, "Training cells (n = 7,573)", size=6.6, color=WHITE,
         weight="bold")
    text(ax, bx + bw * frac + bw * (1 - frac) / 2, ty + 0.11, "Held-out\n(n = 1,921)", size=5.9,
         weight="bold", ls=1.0)

    steps = [
        "All unsupervised fits (HVG, scaling, PCA,\nscVI_matched, totalVI) use training cells only",
        "Held-out cells are encoded after fitting",
        "Logistic regression is trained on training labels",
        "Macro-F1 is evaluated on held-out cells",
    ]
    sy = ty - 0.14
    hs = [0.34, 0.22, 0.22, 0.22]
    for i, (s, hh) in enumerate(zip(steps, hs)):
        sy -= hh + 0.05
        box(ax, bx, sy, bw, hh, fc=PALE, ec=EDGE, r=0.04, z=1)
        ax.add_patch(Circle((bx + 0.13, sy + hh / 2), 0.075, fc=NAVY, ec=NAVY, zorder=3))
        text(ax, bx + 0.13, sy + hh / 2, str(i + 1), size=6.2, color=WHITE, weight="bold")
        text(ax, bx + 0.28, sy + hh / 2, s, size=6.3, ha="left")
    box(ax, bx + 0.25, y + 0.11, bw - 0.50, 0.22, fc=NAVY, ec=NAVY, r=0.05, z=2)
    text(ax, x + w / 2, y + 0.22, "Strict train-only information access (primary)", size=6.8,
         color=WHITE, weight="bold")


def panel_e(ax, x, y, w, h):
    panel(ax, x, y, w, h, "E", "Bidirectional transfer")
    dy = y + h - 0.46
    cw_chip = 1.05
    for i, (lab, s, t) in enumerate([("Direction A", "PBMC10k", "PBMC5k"),
                                      ("Direction B", "PBMC5k", "PBMC10k")]):
        yy = dy - i * 0.28
        text(ax, x + 0.15, yy, lab, size=6.8, weight="bold", ha="left")
        sx = x + 0.92
        box(ax, sx, yy - 0.10, cw_chip, 0.20, fc=NAVY, ec=NAVY, r=0.04, z=2)
        text(ax, sx + cw_chip / 2, yy, f"{s} source", size=6.2, color=WHITE, weight="bold")
        tx = sx + cw_chip + 0.42
        arrow(ax, sx + cw_chip + 0.05, yy, tx - 0.05, yy, color=BLUE, lw=1.0, ms=6)
        box(ax, tx, yy - 0.10, cw_chip, 0.20, fc=LBLUE, ec=BLUE, r=0.04, z=2)
        text(ax, tx + cw_chip / 2, yy, f"{t} target", size=6.2, weight="bold")

    strategies = [
        ("Inductive concat PCA", "source-only fit; target\ntransformed unchanged"),
        ("Feature-access-matched\ntransductive concat", "unsupervised fit on source +\nunlabeled target features"),
        ("scVI_matched + scArches", "unsupervised query\nadaptation"),
        ("totalVI + scArches", "unsupervised query\nadaptation"),
    ]
    cw, chh = (w - 0.40) / 2, 0.44
    gy0 = dy - 0.28 - 0.18
    for k, (name, sub) in enumerate(strategies):
        r_, c_ = divmod(k, 2)
        cx0 = x + 0.15 + c_ * (cw + 0.10)
        cy0 = gy0 - (r_ + 1) * chh - r_ * 0.06
        fc = NAVY if name.startswith("totalVI") else PALE
        box(ax, cx0, cy0, cw, chh, fc=fc, ec=NAVY if fc == NAVY else EDGE, r=0.04, z=1)
        tc = WHITE if fc == NAVY else NAVY
        nlines = name.count("\n") + 1
        text(ax, cx0 + cw / 2, cy0 + chh - 0.06 - 0.05 * nlines, name, size=6.5, weight="bold",
             color=tc, ls=1.0)
        text(ax, cx0 + cw / 2, cy0 + (0.10 if nlines == 2 else 0.13), sub, size=5.9,
             color=(LBLUE if fc == NAVY else GRAY), ls=1.0)
    text(ax, x + 0.15, y + 0.33, "Target labels are never used during adaptation.", size=6.2,
         color=NAVY, style="italic", ha="left", weight="bold")
    text(ax, x + 0.15, y + 0.14,
         "Feature-access-matched concat matches unlabeled target-feature access,\n"
         "not the scArches learning objective or neural adaptation algorithm.",
         size=5.9, color=GRAY, style="italic", ha="left", va="center", ls=1.05)


def mini_axes_curve(ax, x, y, w, h, kind):
    ax.plot([x, x, x + w], [y + h, y, y], color=GRAY, lw=0.6, zorder=3)
    import numpy as np

    t = np.linspace(0, 1, 30)
    if kind == "decline":
        ax.plot(x + t * w, y + h * (0.85 - 0.55 * t ** 1.3), color=NAVY, lw=1.0, zorder=4)
        ax.plot(x + t * w, y + h * (0.55 - 0.05 * t), color=MBLUE, lw=0.9, ls="--", zorder=4)
    else:
        ax.plot(x + t * w, y + h * (0.10 + 0.75 * t ** 2.2), color=NAVY, lw=1.0, zorder=4)
        ax.plot(x + t * w, y + h * (0.30 + 0.62 * t ** 1.6), color=MBLUE, lw=0.9, zorder=4)


def panel_f(ax, x, y, w, h):
    panel(ax, x, y, w, h, "F", "Reliability and robustness analyses")
    cw = (w - 0.40) / 2
    ch = 0.53
    top = y + h - 0.32
    cards = [
        ("Protein corruption", "training-time (retrain) and\npost-adaptation (target only)", "decline"),
        ("Predictive uncertainty",
         r"$U_{\mathrm{pred},i} = 1 - \max_k\, p(y_i = k \mid z_i)$" + "\nfrom the downstream classifier",
         None),
        ("Latent posterior dispersion", "representation-level dispersion;\nnot formal epistemic uncertainty",
         None),
        ("Risk–coverage curves", "selective prediction: defer\nuncertain cells; track retention", "rc"),
    ]
    for k, (name, sub, kind) in enumerate(cards):
        r_, c_ = divmod(k, 2)
        cx0 = x + 0.15 + c_ * (cw + 0.10)
        cy0 = top - (r_ + 1) * ch - r_ * 0.08
        box(ax, cx0, cy0, cw, ch, fc=PALE, ec=EDGE, r=0.04, z=1)
        tx = cx0 + 0.10
        text(ax, tx, cy0 + ch - 0.13, name, size=6.8, weight="bold", ha="left")
        text(ax, tx, cy0 + 0.19, sub, size=5.9, color=GRAY, ha="left", ls=1.1)
        if kind:
            mini_axes_curve(ax, cx0 + cw - 0.48, cy0 + ch - 0.30, 0.38, 0.22, kind)
    sy = y + 0.12
    sh = top - 2 * ch - 0.08 - 0.10 - sy
    box(ax, x + 0.15, sy, w - 0.30, sh, fc=LBLUE, ec=BLUE, lw=0.8, r=0.04, z=1)
    text(ax, x + 0.25, sy + sh - 0.12, "Training stochasticity", size=6.8, weight="bold", ha="left")
    labels = ["source subsets", "model seeds", "run replicates", "variance\ndecomposition"]
    bx = x + 0.25
    bw = (w - 0.50 - 3 * 0.20) / 4
    for i, lab in enumerate(labels):
        xx = bx + i * (bw + 0.20)
        fc = NAVY if i == 3 else WHITE
        box(ax, xx, sy + 0.06, bw, 0.25, fc=fc, ec=NAVY, lw=0.7, r=0.03, z=2)
        text(ax, xx + bw / 2, sy + 0.185, lab, size=6.0, color=(WHITE if i == 3 else NAVY),
             weight="bold", ls=1.0)
        if i < 3:
            sym = "×" if i < 2 else "→"
            text(ax, xx + bw + 0.10, sy + 0.185, sym, size=7.5, weight="bold", color=BLUE)


def donor_icon(ax, x, y, s=0.085, color=BLUE):
    ax.add_patch(Circle((x, y + s * 1.25), s * 0.55, fc=color, ec="none", zorder=3))
    ax.add_patch(FancyBboxPatch((x - s * 0.85, y - s * 0.55), s * 1.7, s * 1.15,
                                boxstyle=f"round,pad=0,rounding_size={s * 0.5}",
                                fc=color, ec="none", zorder=3))


def panel_g(ax, x, y, w, h):
    panel(ax, x, y, w, h, "G", "External donor validation")
    text(ax, x + w / 2, y + h - 0.40, "Lawlor blood CITE-seq: 10 biological donors", size=6.6,
         weight="bold")
    gx0 = x + 0.25
    step = (w - 0.50) / 9
    for i in range(10):
        donor_icon(ax, gx0 + i * step, y + h - 0.70, s=0.07,
                   color=(NAVY if i % 2 == 0 else BLUE))
    items = [
        "Leave-two-donors-out folds; each\ndonor held out exactly once",
        "Aggregation unit: biological donor\n(n = 10)",
        "5 representations: RNA PCA, protein\nPCA, concat PCA, scVI_matched, totalVI",
    ]
    yy = y + h - 1.00
    for it in items:
        nl = it.count("\n") + 1
        ax.add_patch(Circle((x + 0.20, yy), 0.025, fc=NAVY, zorder=3))
        text(ax, x + 0.30, yy - 0.055 * (nl - 1), it, size=6.1, ha="left", va="center", ls=1.1)
        yy -= 0.12 + 0.12 * nl
    box(ax, x + 0.15, y + 0.11, w - 0.30, 0.22, fc=LBLUE, ec=EDGE, r=0.04, z=1)
    text(ax, x + w / 2, y + 0.22, "Lawlor author labels are protein-gated.", size=6.2,
         color=GRAY, style="italic")


def outputs_strip(ax, x, y, w, h):
    box(ax, x, y, w, h, fc=NAVY, ec=NAVY, r=0.05, z=1)
    text(ax, x + 0.14, y + h / 2, "Outputs", size=7.2, color=WHITE, weight="bold", ha="left")
    s = ("Macro-F1  ·  " + G_ACCESS + " / " + G_ASSOC
         + " with paired-bootstrap intervals  ·  transfer contrasts  ·  degradation trajectories  ·"
         "  risk–coverage  ·  donor-level differences")
    text(ax, x + 0.80, y + h / 2, s, size=6.4, color=WHITE, ha="left")


def build():
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    m, gap = 0.05, 0.10
    y_out, h_out = m, 0.30
    y4, h4 = y_out + h_out + gap, 2.25
    y3, h3 = y4 + h4 + gap, 2.35
    y2, h2 = y3 + h3 + gap, 1.80
    y1 = y2 + h2 + gap
    h1 = H - m - y1
    aw = 2.15
    panel_a(ax, m, y1, aw, h1)
    panel_b(ax, m + aw + gap, y1, W - 2 * m - aw - gap, h1)
    panel_c(ax, m, y2, W - 2 * m, h2)
    dw = 3.05
    panel_d(ax, m, y3, dw, h3)
    panel_e(ax, m + dw + gap, y3, W - 2 * m - dw - gap, h3)
    fw = 4.05
    panel_f(ax, m, y4, fw, h4)
    panel_g(ax, m + fw + gap, y4, W - 2 * m - fw - gap, h4)
    outputs_strip(ax, m, y_out, W - 2 * m, h_out)
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--stem", default=STEM)
    args = ap.parse_args()
    out_dir = args.out_dir if args.out_dir.is_absolute() else ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = build()
    for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 600})):
        out = out_dir / f"{args.stem}.{ext}"
        fig.savefig(out, facecolor="white", **kw)
        print(f"wrote {out.relative_to(ROOT)}")
    plt.close(fig)


if __name__ == "__main__":
    main()
