#!/usr/bin/env python3
"""Generate Figure 1 for the JCM manuscript: v1→v2→v3 error-mode arc.

Pure stdlib — no external dependencies.
Usage: python3 scripts/figure1_error_arc.py [output_path]
Default output: work/figure_1.svg
"""
import os
import sys

# ── data ──────────────────────────────────────────────────────────────────────
VERSIONS = ["v1", "v2", "v3"]

# (point, Wilson-lower, Wilson-upper) in percent
DATA = {
    "v1": {"vme": (10.0,  4.0, 23.1), "me": ( 2.1,  0.4, 11.1), "abs": 11.2},
    "v2": {"vme": ( 0.0,  0.0,  7.9), "me": (13.5,  6.7, 25.3), "abs":  1.0},
    "v3": {"vme": ( 0.0,  0.0,  8.4), "me": ( 2.2,  0.4, 11.3), "abs": 10.2},
}

VME_BAR_POINT = 3.0   # pass bar: VME point estimate ≤ 3%
VME_BAR_UPPER = 15.0  # pass bar: VME Wilson-upper ≤ 15%
ME_BAR_POINT  = 5.0   # pass bar: ME point estimate ≤ 5%

# ── layout ────────────────────────────────────────────────────────────────────
FIG_W = 680
FIG_H = 370
ML, MR, MT, MB = 65, 120, 50, 70  # margins (MR wide for threshold labels)
PW = FIG_W - ML - MR               # plot width  = 590
PH = FIG_H - MT - MB               # plot height = 265
Y_MAX = 30.0

BAR_W = 40
GAP = 10   # gap between bars in a group

# ── coordinate helpers ────────────────────────────────────────────────────────
def px(x_data):
    return ML + PW * x_data / len(VERSIONS)

def py(y_pct):
    return MT + PH - (y_pct / Y_MAX) * PH

GROUP_CENTERS = [i + 0.5 for i in range(len(VERSIONS))]

def bar_left(group_i, bar_i):
    cx = px(GROUP_CENTERS[group_i])
    span = 3 * BAR_W + 2 * GAP
    return cx - span / 2 + bar_i * (BAR_W + GAP)

# ── colors ────────────────────────────────────────────────────────────────────
C_VME   = "#c0392b"
C_ME    = "#2980b9"
C_ABS   = "#aab7c0"
C_GRID  = "#e0e4e8"
C_AXIS  = "#444444"
C_FAIL  = "#c0392b"
C_PASS  = "#27ae60"

# ── SVG element builders ──────────────────────────────────────────────────────
def rect(x, y, w, h, fill, opacity=1.0, rx=2):
    return (f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'fill="{fill}" opacity="{opacity}" rx="{rx}"/>')

def seg(x1, y1, x2, y2, stroke, sw=1.0, dash=""):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return (f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{sw}"{d}/>')

def txt(x, y, s, anchor="middle", size=11, fill="#222", bold=False, rotate=None):
    w = ' font-weight="bold"' if bold else ""
    rot = (f' transform="rotate({rotate[0]},{rotate[1]:.1f},{rotate[2]:.1f})"'
           if rotate else "")
    return (f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" '
            f'font-size="{size}" fill="{fill}"{w}{rot}>{s}</text>')

# ── assemble ──────────────────────────────────────────────────────────────────
E = []

E.append(f'<svg xmlns="http://www.w3.org/2000/svg" '
         f'width="{FIG_W}" height="{FIG_H}" '
         f'font-family="Helvetica Neue, Arial, sans-serif">')
E.append(f'<rect width="{FIG_W}" height="{FIG_H}" fill="white"/>')

# figure title
E.append(txt(ML + PW / 2, MT - 28,
             "FKS1 caller: verdict accuracy across three pre-registered measurements",
             size=12, bold=True))

# y gridlines + tick labels
for yt in range(0, 31, 5):
    yp = py(yt)
    E.append(seg(ML, yp, ML + PW, yp, C_GRID, sw=0.8))
    E.append(txt(ML - 8, yp + 4, str(yt), anchor="end", size=10, fill=C_AXIS))

# y-axis label (rotated)
lx = ML - 50
ly = MT + PH / 2
E.append(txt(lx, ly, "Rate (%)", size=11, fill=C_AXIS, rotate=(-90, lx, ly)))

# axis frame
E.append(seg(ML, MT, ML, MT + PH, C_AXIS, sw=1.2))
E.append(seg(ML, MT + PH, ML + PW, MT + PH, C_AXIS, sw=1.2))

# threshold lines (drawn behind bars; labels anchored inside left)
def threshold_line(y_pct, stroke, dash, label, label_dy=0):
    yp = py(y_pct)
    E.append(seg(ML, yp, ML + PW, yp, stroke, sw=1.2, dash=dash))
    E.append(txt(ML + PW + 5, yp - 4 + label_dy, label,
                 anchor="start", size=9, fill=stroke))

threshold_line(VME_BAR_POINT,  C_VME, "5,3",  f"VME ≤3%",  label_dy=12)
threshold_line(VME_BAR_UPPER,  C_VME, "3,4",  f"VME upper ≤15%")
threshold_line(ME_BAR_POINT,   C_ME,  "5,3",  f"ME ≤5%",   label_dy=-1)

# bars + error bars
for gi, ver in enumerate(VERSIONS):
    d = DATA[ver]

    for bi, (key, color) in enumerate([("vme", C_VME), ("me", C_ME), ("abs", C_ABS)]):
        bx   = bar_left(gi, bi)
        bcx  = bx + BAR_W / 2
        ybot = py(0)

        if key == "abs":
            val = d["abs"]
            ytop = py(val)
            h = ybot - ytop
            E.append(rect(bx, ytop, BAR_W, h, color, opacity=0.65))
        else:
            pt, lo, hi = d[key]
            ytop = py(pt)
            h = ybot - ytop
            if h < 0.5:
                h = 0.5   # zero-height placeholder
                ytop = ybot - 0.5
            E.append(rect(bx, ytop, BAR_W, h, color))
            # error bar (Wilson CI)
            ylo_px = py(lo)
            yhi_px = py(hi)
            cap = 5
            E.append(seg(bcx, ylo_px, bcx, yhi_px, "#333", sw=1.3))
            E.append(seg(bcx - cap/2, ylo_px, bcx + cap/2, ylo_px, "#333", sw=1.3))
            E.append(seg(bcx - cap/2, yhi_px, bcx + cap/2, yhi_px, "#333", sw=1.3))

    # version label below axis
    gcx = px(GROUP_CENTERS[gi])
    yaxis = py(0)
    E.append(txt(gcx, yaxis + 18, ver, size=13, bold=True, fill=C_AXIS))

    # PASS/FAIL verdict
    vme_pt = DATA[ver]["vme"][0]
    vme_hi = DATA[ver]["vme"][2]
    me_pt  = DATA[ver]["me"][0]
    passed = (vme_pt <= VME_BAR_POINT
              and vme_hi <= VME_BAR_UPPER
              and me_pt  <= ME_BAR_POINT)
    label  = "PASS" if passed else "FAIL"
    color  = C_PASS  if passed else C_FAIL
    E.append(txt(gcx, yaxis + 35, label, size=11, bold=True, fill=color))

# legend (top-right inside plot)
legend_items = [
    (C_VME, "VME (very-major error)"),
    (C_ME,  "ME (major error)"),
    (C_ABS, "Abstention rate"),
]
lx_start = ML + PW - 370
ly_base  = MT + 12
for i, (c, label) in enumerate(legend_items):
    lx = lx_start + i * 127
    op = 0.65 if c == C_ABS else 1.0
    E.append(rect(lx, ly_base - 9, 14, 10, c, opacity=op, rx=1))
    E.append(txt(lx + 17, ly_base, label, anchor="start", size=10, fill=C_AXIS))

E.append("</svg>")

# ── write ─────────────────────────────────────────────────────────────────────
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "work", "figure_1.svg")
if len(sys.argv) > 1:
    out = sys.argv[1]

os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
with open(out, "w", encoding="utf-8") as fh:
    fh.write("\n".join(E) + "\n")

print(f"Written: {out}")
