"""Map: CZ-2D Y108 powered ascent, NOTAM zone, T0+6 min sighting geometry and the observer locus.

SatObserver-MX look via ~/sda/starship-ITF-14/work/ift14_style.py (Blue Marble NG base, Menlo, halo labels).
Launch-plot export rules from ~/.claude/CLAUDE.md: 300 dpi, 3960 px long side, title + UTC upper left, two-row
legend, credit below the map, solid Sun -12 deg line (cyan) at the reference epoch (T0+6 min); the Sun 0 deg line
lies outside the map.

Run: ~/.venvs/astro313/bin/python make_map.py   -> CZ-2D_Y108_observer_map.png
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.home() / "sda/starship-ITF-14/work"))
import ift14_style as st  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Polygon  # noqa: E402

HERE = Path(__file__).resolve().parent / "missions" / "Y108_JSLC_20261111"
EXT = (90.0, 105.0, 30.5, 42.5)
INSET = (86.0, 106.0, 2.0, 44.0)
REF = np.datetime64("2026-11-11T11:47:00")
TOWNS = [("Jiuquan", 39.744, 98.494, "left"), ("Jiayuguan", 39.773, 98.289, "right"), ("Jinta", 39.983, 98.902, "left"),
         ("Zhangye", 38.926, 100.451, "left"), ("Xining", 36.617, 101.778, "left"), ("Lanzhou", 36.061, 103.834, "left"),
         ("Golmud", 36.402, 94.903, "left"), ("Delingha", 37.370, 97.370, "right"), ("Dunhuang", 40.142, 94.662, "left"),
         ("Yushu", 33.004, 97.007, "right"), ("Ejina", 41.958, 101.069, "left"), ("Lenghu", 38.743, 93.338, "left")]
C_S1, C_S2, C_VN, C_LOC, C_OBS = "#ff7043", "#ffb84f", "#ffe7a3", "#7ad97a", "#ff5252"


def gc_points(a, b, n=60):
    la1, lo1, la2, lo2 = map(np.radians, (a[0], a[1], b[0], b[1]))
    p = np.array([np.cos(la1) * np.cos(lo1), np.cos(la1) * np.sin(lo1), np.sin(la1)])
    q = np.array([np.cos(la2) * np.cos(lo2), np.cos(la2) * np.sin(lo2), np.sin(la2)])
    w = np.arccos(np.clip(p @ q, -1, 1))
    s = np.linspace(0, 1, n)
    v = (np.sin((1 - s)[:, None] * w) * p + np.sin(s[:, None] * w) * q) / np.sin(w)
    return np.degrees(np.arcsin(v[:, 2])), np.degrees(np.arctan2(v[:, 1], v[:, 0]))


def night(ax, ext):
    x0, x1, y0, y1 = ext
    gl, gt = np.meshgrid(np.linspace(x0, x1, 901), np.linspace(y0, y1, 721))
    alt = st.sun_alt(REF, gl, gt)
    sh = np.zeros(gl.shape + (4,))
    sh[..., 2] = 10 / 255
    sh[..., 3] = np.where(alt < -12, 0.48, np.where(alt < 0, 0.32, 0.0))
    ax.imshow(sh, extent=ext, aspect="auto", origin="lower", zorder=1, interpolation="bilinear")
    if alt.max() > 0 > alt.min():
        ax.contour(gl, gt, alt, levels=[0], colors=[st.SUN], linewidths=1.5, linestyles="solid", zorder=3)
    if alt.max() > -12 > alt.min():
        ax.contour(gl, gt, alt, levels=[-12], colors=[st.ACCENT], linewidths=1.3, linestyles="solid", zorder=3)


def main():
    sol = json.load(open(HERE / "CZ-2D_Y108_solution.json"))
    obs = json.load(open(HERE / "observer_solution.json"))
    rows = list(csv.DictReader(open(HERE / "CZ-2D_Y108_traj.csv")))
    loc = list(csv.DictReader(open(HERE / "observer_locus.csv")))
    ev = {e["event"]: e for e in sol["events"]}
    t_pick = obs["sunlit_from_met"]["geometric"]

    x0, x1, y0, y1 = EXT
    LONG, top, bot, side = 13.2, 0.42, 0.66, 0.16
    kx = np.cos(np.radians((y0 + y1) / 2))
    ax_h = LONG - top - bot
    ax_w = ax_h * (x1 - x0) * kx / (y1 - y0)
    W, H = ax_w + 2 * side, LONG
    if W > LONG:
        W = LONG
        ax_w = W - 2 * side
        ax_h = ax_w * (y1 - y0) / ((x1 - x0) * kx)
        H = ax_h + top + bot
    fig = plt.figure(figsize=(W, H), facecolor=st.BG)
    ax = fig.add_axes([side / W, bot / H, ax_w / W, ax_h / H])
    ax.set_facecolor("#05070a")
    for s in ax.spines.values():
        s.set_edgecolor("#2a323c")
    title = "CZ-2D Y108 Ascent and observer"
    fig.text(side / W, 1 - 0.5 * top / H, title, color="#cfd6dd", fontsize=13, fontweight="bold", va="center")
    fig.text((side + 0.112 * len(title) + 0.25) / W, 1 - 0.5 * top / H, "2026-11-11 11:41 UTC", color=st.FG,
             fontsize=12, va="center")
    st.draw_base(ax, EXT, REF, grat=5, night=False, sun=False, label_step=5)
    night(ax, EXT)

    # NOTAM zone, pad, impacts
    ring = sol["notam"]["ring"]
    ax.add_patch(Polygon([(p[1], p[0]) for p in ring], closed=True, fc=(0.7, 0.7, 0.7, 0.35), ec="#d0d4d8", lw=1.0,
                         zorder=4))
    st.star(ax, sol["pad"][1], sol["pad"][0], st.WARN, ms=15)
    st.label(ax, sol["pad"][1] + 0.18, sol["pad"][0] + 0.12, "JSLC", size=9)
    la1, lo1 = sol["stage1"]["impact"]
    st.label(ax, lo1 - 0.3, la1 - 0.05, f"{sol['notam']['id']}\nstage 1", size=8.5, color="#d0d4d8", ha="right")
    laf, lof = sol["fairing_impact"]
    ax.plot(lof, laf, "s", ms=6, mfc="none", mec="#d0d4d8", mew=1.2, zorder=6)
    st.label(ax, lof + 0.2, laf, "fairing", size=8.5, color="#d0d4d8")

    # powered track by phase; vernier dashed in shadow, solid once sunlit
    def seg(cond):
        r = [x for x in rows if cond(x)]
        return [float(x["lon"]) for x in r], [float(x["lat"]) for x in r]
    lo, la = seg(lambda r: r["phase"] in ("vertical rise", "stage 1"))
    ax.plot(lo, la, color=C_S1, lw=2.4, zorder=5)
    lo, la = seg(lambda r: r["phase"] == "stage 2 main+vernier")
    ax.plot(lo, la, color=C_S2, lw=2.4, zorder=5)
    lo, la = seg(lambda r: r["phase"] == "stage 2 vernier" and float(r["met"]) <= t_pick)
    ax.plot(lo, la, color=C_VN, lw=2.0, ls=(0, (4, 3)), zorder=5)
    lo, la = seg(lambda r: r["phase"] == "stage 2 vernier" and float(r["met"]) >= t_pick)
    ax.plot(lo, la, color=C_VN, lw=2.2, zorder=5)
    m2 = ev["main-engine cutoff (MECO-2)"]
    ax.plot(m2["lon"], m2["lat"], "o", ms=5, mfc=C_S2, mec="#101418", zorder=7)
    st.label(ax, m2["lon"] + 0.3, m2["lat"] + 0.1, f"MECO-2 T+{m2['met']:.0f} s", size=8.5)

    # rocket at T+6 min (sunlight pickup) and line of sight to the observer
    b = obs["best"]
    rk = (b["rocket_lat"], b["rocket_lon"])
    gla, glo = gc_points((b["lat"], b["lon"]), rk)
    ax.plot(glo, gla, color="#ffffff", lw=1.0, ls=(0, (2, 3)), alpha=0.85, zorder=5)
    st.star(ax, rk[1], rk[0], C_VN, ms=17)
    st.label(ax, rk[1] + 0.28, rk[0] - 0.02, f"T+6:00  {b['rocket_alt_km']:.0f} km\nsunlit from here", size=9)

    # observer locus with 10 s ticks
    L = [(float(p["met"]), float(p["lat"]), float(p["lon"])) for p in loc]
    ax.plot([p[2] for p in L if p[0] >= t_pick], [p[1] for p in L if p[0] >= t_pick], color=C_LOC, lw=2.2, zorder=6)
    ax.plot([p[2] for p in L if p[0] <= t_pick], [p[1] for p in L if p[0] <= t_pick], color=C_LOC, lw=1.6,
            ls=(0, (1.5, 2)), zorder=6)
    for m, la_, lo_ in L:
        if abs(m % 10) < 1e-6 and m != 360:
            ax.plot(lo_, la_, "o", ms=3.2, mfc=C_LOC, mec="#101418", mew=0.6, zorder=7)
            if m in (300.0, 330.0, 390.0, 420.0):
                st.label(ax, lo_ + 0.16, la_, f"T+{m:.0f}", size=8, color=C_LOC)
    ax.plot(b["lon"], b["lat"], "o", ms=9, mfc=C_OBS, mec="#ffffff", mew=1.2, zorder=8)
    st.label(ax, b["lon"] + 0.22, b["lat"] - 0.16,
             f"Observer {b['lat']:.3f}°N {b['lon']:.3f}°E\naz {b['az']:.1f}° el {b['el']:.1f}°  {b['range_km']:.0f} km",
             size=9, va="top")

    for name, la_, lo_, side_ in TOWNS:
        if x0 < lo_ < x1 and y0 < la_ < y1:
            ax.plot(lo_, la_, "o", ms=3.0, mfc="#e8eaed", mec="#101418", mew=0.6, zorder=6)
            dx = 0.12 if side_ == "left" else -0.12
            st.label(ax, lo_ + dx, la_ + 0.08, name, size=8.5, ha="left" if side_ == "left" else "right",
                     color="#e8eaed")

    # inset: whole ascent to insertion
    ix0, ix1, iy0, iy1 = INSET
    ih = 0.30 * ax_h
    iw = ih * (ix1 - ix0) * np.cos(np.radians(23)) / (iy1 - iy0)
    ia = fig.add_axes([(side + 0.55) / W, (bot + ax_h - ih - 0.25) / H, iw / W, ih / H])
    st.draw_base(ia, INSET, REF, grat=10, night=False, sun=False, label_step=20)
    night(ia, INSET)
    for t in ia.texts:
        t.set_fontsize(6.5)
    for s in ia.spines.values():
        s.set_edgecolor("#5a6470")
    for ph, col in (("stage 1", C_S1), ("stage 2 main+vernier", C_S2), ("stage 2 vernier", C_VN)):
        lo, la = seg(lambda r, ph=ph: r["phase"] == ph or (ph == "stage 1" and r["phase"] == "vertical rise"))
        ia.plot(lo, la, color=col, lw=1.6, zorder=5)
    lo, la = seg(lambda r: r["phase"] == "orbit")
    ia.plot(lo, la, color="#c8cdd2", lw=1.0, ls=(0, (2, 2)), zorder=5)
    ins = ev["vernier cutoff / insertion"]
    ia.plot(ins["lon"], ins["lat"], "o", ms=4, mfc=C_VN, mec="#101418", zorder=7)
    st.label(ia, ins["lon"] + 1.0, ins["lat"], f"insertion\nT+{ins['met'] / 60:.0f} min\n600 km", size=7)
    ia.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="#ffffff", lw=0.8, zorder=8))

    hud = [("Sighting  RA 21h33m32.7s  Dec −31°11′16″", st.FG, True),
           (f"T+6:00 {b['utc'][11:19]} UTC  ·  ground {b['h_km'] * 1e3:.0f} m  ·  Sun {obs['sun_alt_deg']:.0f}°", st.FGDIM, False),
           ("Locus: ±10 s ≈ ±40 km along track", st.FGDIM, False)]
    st.hud(ax, hud, x0=0.685, y0=0.045, lh=0.02, width=0.305, fontsize=8.5)

    handles = [Line2D([], [], color=C_S1, lw=2.4, label="Stage 1"),
               Line2D([], [], color=C_S2, lw=2.4, label="Stage 2 main"),
               Line2D([], [], color=C_VN, lw=2.2, label="Vernier (sunlit)"),
               Line2D([], [], color=C_VN, lw=2.0, ls=(0, (4, 3)), label="Vernier (shadow)"),
               Patch(fc=(0.7, 0.7, 0.7, 0.35), ec="#d0d4d8", label="NOTAM zone"),
               Line2D([], [], color=C_LOC, lw=2.2, label="Observer locus"),
               Line2D([], [], color="#ffffff", lw=1.0, ls=(0, (2, 3)), label="Line of sight"),
               Line2D([], [], color=st.ACCENT, lw=1.3, label="Sun −12°")]
    lg = fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(side / W, 0.06 / H), ncol=4, frameon=False,
                    fontsize=8.5, handlelength=2.4, columnspacing=1.6, labelcolor=st.FG)
    lg.set_zorder(12)
    fig.text(1 - side / W, 0.10 / H, "Illustration: Mickey · base map NASA Blue Marble", color=st.FGDIM,
             fontsize=7.5, ha="right", va="bottom")
    out = HERE / "CZ-2D_Y108_observer_map.png"
    dpi = 3960 / max(W, H)
    fig.savefig(out, dpi=dpi, facecolor=st.BG)
    print(f"wrote {out.name}  {W:.2f} x {H:.2f} in at {dpi:.0f} dpi")


if __name__ == "__main__":
    main()
