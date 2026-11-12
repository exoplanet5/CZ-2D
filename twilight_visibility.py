"""Morning-twilight ("jellyfish") geometry and launch-visibility map for a build_traj.py mission.

Reads <stem>_traj.csv and <stem>_solution.json from missions/<dir>/ (build_traj.MISSIONS) and writes there.  For every second of flight the Sun-ray tangent height is the
lowest altitude (spherical 6371 km Earth) of the ray from the vehicle toward the Sun: >= 0 km is geometrically
sunlit, >= 20 km is clear of the opaque lower atmosphere (same convention as observer_locate.py).  The Sun is
evaluated at each point's own epoch (astropy).  The visibility zone is every ground point that sees the sunlit
vehicle during the main-engine burn (pickup -> MECO-2) at >= --min-elev, clipped to the night side (Sun < 0 deg) at
the reference epoch = middle of that sunlit main-engine arc.  Solar contours on the map use the same epoch.

Map: SatObserver-MX look (ift14_style) with the launch-plot export rules of ~/.claude/CLAUDE.md.

Run: ~/.venvs/astro313/bin/python twilight_visibility.py [--mission tslc1118] [--min-elev 10]
Writes <stem>_twilight.json, <stem>_visibility.png
"""
import argparse
import csv
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import ITRS, get_sun
from astropy.time import Time
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path.home() / "sda/starship-ITF-14/work"))
import ift14_style as st  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patheffects as pe  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
A_WGS, E2 = 6378.137, 6.69437999014e-3
CITIES = [("Beijing", 39.90, 116.40, "l"), ("Taiyuan", 37.87, 112.55, "l"), ("Xi'an", 34.27, 108.94, "r"),
          ("Zhengzhou", 34.75, 113.63, "l"), ("Wuhan", 30.59, 114.31, "l"), ("Chongqing", 29.56, 106.55, "r"),
          ("Chengdu", 30.66, 104.07, "r"), ("Lanzhou", 36.06, 103.83, "r"), ("Changsha", 28.23, 112.94, "l"),
          ("Guiyang", 26.65, 106.63, "r"), ("Yinchuan", 38.49, 106.23, "r"), ("Hohhot", 40.84, 111.75, "l"),
          ("Kunming", 25.04, 102.71, "r"), ("Nanning", 22.82, 108.32, "r"), ("Guangzhou", 23.13, 113.26, "l"),
          ("Shanghai", 31.23, 121.47, "l"), ("Jinan", 36.65, 117.12, "l"), ("Nanjing", 32.06, 118.80, "r"),
          ("Xining", 36.62, 101.78, "r"), ("Shijiazhuang", 38.04, 114.51, "l"), ("Yan'an", 36.59, 109.49, "r"),
          ("Haikou", 20.04, 110.34, "l"), ("Hanoi", 21.03, 105.85, "r"), ("Fuzhou", 26.07, 119.30, "r"),
          ("Nanchang", 28.68, 115.86, "l"), ("Baotou", 40.66, 109.84, "r"), ("Ulaanbaatar", 47.92, 106.92, "l")]
C_SH, C_LIT, C_VN, C_VIS = st.WARN, "#ffffff", "#ffe7a3", st.OK


def geodetic_ecef(lat, lon, h=0.0):
    la, lo = np.radians(lat), np.radians(lon)
    n = A_WGS / np.sqrt(1 - E2 * np.sin(la) ** 2)
    return np.stack([(n + h) * np.cos(la) * np.cos(lo), (n + h) * np.cos(la) * np.sin(lo),
                     (n * (1 - E2) + h) * np.sin(la)], axis=-1)


def up_vec(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def sun_itrs(t0, mets):
    t = Time(t0, scale="utc") + np.asarray(mets, float) * u.s
    s = get_sun(t).transform_to(ITRS(obstime=t)).cartesian.xyz.to_value(u.km).T
    return s / np.linalg.norm(s, axis=1)[:, None]


def tangent_km(X, S):
    """Lowest altitude of the ray from X toward the unit Sun vector S (rows)."""
    k = np.maximum(-(X * S).sum(1), 0.0)
    return np.linalg.norm(X + k[:, None] * S, axis=1) - 6371.0


def sun_alt_at(lat, lon, S):
    return np.degrees(np.arcsin(up_vec(lat, lon) @ S))


def look(lat, lon, X):
    """Elevation, azimuth (deg) and range (km) of ECEF points X from a sea-level site."""
    O, up = geodetic_ecef(lat, lon), up_vec(lat, lon)
    d = X - O
    rho = np.linalg.norm(d, axis=-1)
    e = np.array([-np.sin(np.radians(lon)), np.cos(np.radians(lon)), 0.0])
    n = np.cross(up, e)
    return np.degrees(np.arcsin(d @ up / rho)), np.degrees(np.arctan2(d @ e, d @ n)) % 360, rho


def first_cross(T, y, thr):
    i = np.argmax(y >= thr)
    if y[i] < thr:
        return None
    if i == 0:
        return float(T[0])
    return float(T[i - 1] + (thr - y[i - 1]) / (y[i] - y[i - 1]) * (T[i] - T[i - 1]))


def analyse(out_dir, stem, min_elev):
    sol = json.load(open(out_dir / f"{stem}_solution.json"))
    rows = list(csv.DictReader(open(out_dir / f"{stem}_traj.csv")))
    t0 = datetime.fromisoformat(sol["t0_utc"])
    ev = {e["event"]: e for e in sol["events"]}
    T = np.array([float(r["met"]) for r in rows])
    X = np.array([[float(r[k]) for k in ("x_ecef", "y_ecef", "z_ecef")] for r in rows])
    lat = np.array([float(r["lat"]) for r in rows])
    lon = np.array([float(r["lon"]) for r in rows])
    alt = np.array([float(r["alt_km"]) for r in rows])
    S = sun_itrs(t0, T)
    tan = tangent_km(X, S)
    t_meco1, t_meco2 = ev["MECO-1 / hot staging"]["met"], ev["main-engine cutoff (MECO-2)"]["met"]
    t_ins = ev["vernier cutoff / insertion"]["met"]
    t_geo, t_clear = first_cross(T, tan, 0.0), first_cross(T, tan, 20.0)

    def at(t):
        return dict(met=round(t, 1), utc=(t0 + timedelta(seconds=t)).isoformat(timespec="seconds"),
                    lat=float(np.interp(t, T, lat)), lon=float(np.interp(t, T, lon)),
                    alt_km=float(np.interp(t, T, alt)), sun_tangent_km=float(np.interp(t, T, tan)))

    pick, clear = at(t_geo), at(t_clear)
    t_ref = 0.5 * (max(t_clear, t_meco1) + t_meco2)
    S_ref = sun_itrs(t0, [t_ref])[0]
    ref_utc = t0 + timedelta(seconds=t_ref)
    ground_sun = np.array([sun_alt_at(lat[i], lon[i], S[i]) for i in range(len(T))])

    sel = (T >= t_clear) & (T <= t_meco2)                  # sunlit main-engine arc (stage-1 tail included if lit)
    selv = (T > t_meco2) & (T <= t_ins)
    cities = []
    for name, la, lo, _ in CITIES:
        el, az, rho = look(la, lo, X[sel])
        i = int(np.argmax(el))
        elv = look(la, lo, X[selv])[0]
        cities.append(dict(city=name, lat=la, lon=lo, sun_alt=round(float(sun_alt_at(la, lo, S_ref)), 1),
                           max_el=round(float(el[i]), 1), az=round(float(az[i])), range_km=round(float(rho[i])),
                           met_max_el=float(T[sel][i]), el_meco1=round(float(look(la, lo, X[T == np.floor(t_meco1)])[0][0]), 1),
                           el_meco2=round(float(el[-1]), 1), az_meco2=round(float(az[-1])),
                           vernier_max_el=round(float(elv.max()), 1)))
    out = dict(stem=stem, t0_utc=sol["t0_utc"], sun_at_pad_T0=round(float(sun_alt_at(*sol["pad"], S[0])), 2),
               pickup_geometric=pick, pickup_clear_20km=clear, ref_epoch_met=round(t_ref, 1),
               ref_epoch_utc=ref_utc.isoformat(timespec="seconds"), min_elev=min_elev,
               ground_sun_alt=dict(MECO1=round(float(np.interp(t_meco1, T, ground_sun)), 2),
                                   MECO2=round(float(np.interp(t_meco2, T, ground_sun)), 2),
                                   insertion=round(float(np.interp(t_ins, T, ground_sun)), 2)),
               sun_tangent_km=dict(MECO1=round(float(np.interp(t_meco1, T, tan)), 1),
                                   MECO2=round(float(np.interp(t_meco2, T, tan)), 1)), cities=cities)
    return sol, rows, dict(T=T, X=X, lat=lat, lon=lon, alt=alt, tan=tan, sel=sel, S_ref=S_ref, ref_utc=ref_utc,
                           t_meco1=t_meco1, t_meco2=t_meco2, t_ins=t_ins, t_geo=t_geo, t_clear=t_clear), out


def vis_grid(ext, D, min_elev, step=0.02, dt=0.25):
    """Max elevation of the sunlit main-engine arc and Sun altitude on a step-deg grid.  The arc is resampled
    every dt s (cubic) so the envelope is not scalloped between the 1 s trajectory rows."""
    x0, x1, y0, y1 = ext
    gl, gt = np.meshgrid(np.arange(x0, x1 + 1e-9, step), np.arange(y0, y1 + 1e-9, step))
    O, up = geodetic_ecef(gt, gl).astype(np.float32), up_vec(gt, gl).astype(np.float32)
    tt = np.arange(D["t_clear"], D["t_meco2"] + 1e-9, dt)
    Xs = CubicSpline(D["T"], D["X"])(tt).astype(np.float32)
    best = np.full(gl.shape, -1.0, dtype=np.float32)          # max sin(elevation)
    for x in Xs:
        d = x - O
        np.maximum(best, (d * up).sum(-1) / np.sqrt((d * d).sum(-1)), out=best)
    sun = np.degrees(np.arcsin(up.astype(float) @ D["S_ref"]))
    return gl, gt, np.degrees(np.arcsin(best.astype(float))), sun


def render(out_dir, stem, sol, rows, D, out, min_elev, ext, inset_ext):
    x0, x1, y0, y1 = ext
    LONG, top, bot, side = 13.2, 0.42, 0.62, 0.16
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
    title = "CZ-2D Launch visibility"
    fig.text(side / W, 1 - 0.5 * top / H, title, color="#cfd6dd", fontsize=13, fontweight="bold", va="center")
    fig.text((side + 0.112 * len(title) + 0.25) / W, 1 - 0.5 * top / H,
             datetime.fromisoformat(sol["t0_utc"]).strftime("%Y-%m-%d %H:%M UTC"), color=st.FG, fontsize=12,
             va="center")
    ref = np.datetime64(D["ref_utc"].isoformat())

    def base(a, e, grat, lstep):
        st.draw_base(a, e, ref, grat=grat, night=False, sun=False, label_step=lstep)
        gl, gt = np.meshgrid(np.linspace(e[0], e[1], 901), np.linspace(e[2], e[3], 721))
        sa = np.degrees(np.arcsin(up_vec(gt, gl) @ D["S_ref"]))
        sh = np.zeros(gl.shape + (4,))
        sh[..., 2] = 10 / 255
        sh[..., 3] = np.where(sa < -12, 0.48, np.where(sa < 0, 0.32, 0.0))
        a.imshow(sh, extent=e, aspect="auto", origin="lower", zorder=1, interpolation="bilinear")
        if sa.max() > 0 > sa.min():
            a.contour(gl, gt, sa, levels=[0], colors=[st.SUN], linewidths=1.5, linestyles="solid", zorder=3)
        if sa.max() > -12 > sa.min():
            a.contour(gl, gt, sa, levels=[-12], colors=[st.ACCENT], linewidths=1.3, linestyles="solid", zorder=3)

    base(ax, ext, 5, 5)
    gl, gt, el, sun = vis_grid(ext, D, min_elev)
    # zone = {elevation >= min_elev} and {Sun < 0}: the zero level of this continuous field is its exact boundary
    zone = np.minimum(el - min_elev, -sun)
    ax.contourf(gl, gt, zone, levels=[0.0, 1e3], colors=[(1, 1, 1, 0.10)], zorder=2)
    ax.contour(gl, gt, zone, levels=[0.0], colors=[C_VIS], linewidths=1.4, linestyles=[(0, (5, 3))], zorder=4)

    T, lat, lon = D["T"], D["lat"], D["lon"]

    def track(a, lo_t, hi_t, **kw):
        m = (T >= lo_t - 0.5) & (T <= hi_t + 0.5)
        a.plot(lon[m], lat[m], **kw)

    def draw_track(a, k=1.0):
        track(a, D["t_meco2"], D["t_ins"], color=C_VN, lw=1.8 * k, ls=(0, (4, 2.5)), zorder=6)
        track(a, 0, D["t_clear"], color=C_SH, lw=2.6 * k, solid_capstyle="round", zorder=7)
        track(a, D["t_clear"], D["t_meco2"], color=C_LIT, lw=3.6 * k, solid_capstyle="round", zorder=8,
              path_effects=[pe.withStroke(linewidth=7 * k, foreground=(1, 1, 1, 0.25))])

    for q in sol["zones"]:
        ring = q["ring"] + q["ring"][:1]
        ax.fill([p[1] for p in ring], [p[0] for p in ring], color=(1, 82 / 255, 82 / 255, 0.25), zorder=5)
        ax.plot([p[1] for p in ring], [p[0] for p in ring], color=st.DANGER, lw=1.3, zorder=5)
    draw_track(ax)
    for la_, lo_ in (sol["stage1"]["impact"], sol["fairing_impact"]):
        ax.plot(lo_, la_, "x", ms=6, mew=1.6, color=st.DANGER, zorder=9)
    pad = sol["pad"]
    st.star(ax, pad[1], pad[0], st.ACCENT, ms=15)
    st.label(ax, pad[1] + 0.35, pad[0] + 0.3, "TSLC", color=st.ACCENT, size=9.5)

    for cname, la_, lo_, sd in CITIES:
        if not (x0 + 0.4 < lo_ < x1 - 0.4 and y0 + 0.4 < la_ < y1 - 0.4):
            continue
        ax.plot(lo_, la_, "o", ms=3.2, mfc="#d8dde3", mec=(0, 0, 0, 0.7), mew=0.6, zorder=6)
        st.label(ax, lo_ + (0.22 if sd == "l" else -0.22), la_ + 0.05, cname, color="#d8dde3", size=8,
                 ha="left" if sd == "l" else "right")

    def mark(t, text, dx=0.35, dy=0.0, ha="left", **kw):
        la_, lo_ = np.interp(t, T, lat), np.interp(t, T, lon)
        ax.plot(lo_, la_, "o", ms=3.8, mfc=st.BG, mec="#ffffff", mew=1.0, zorder=10)
        st.label(ax, lo_ + dx, la_ + dy, text, color="#e8eaed", size=8, ha=ha, **kw)

    mark(D["t_meco2"], f"MECO-2 T+{D['t_meco2']:.0f} s", dy=-0.05)
    for tm in (480.0, 600.0, 720.0):
        if y0 + 0.5 < np.interp(tm, T, lat):
            mark(tm, f"T+{tm / 60:.0f} min")
    pk = out["pickup_clear_20km"]
    ax.plot(pk["lon"], pk["lat"], "s", ms=8, mfc=st.OK, mec=(0, 0, 0, 0.8), mew=1, zorder=11)
    st.label(ax, pk["lon"] - 0.4, pk["lat"] - 0.05, f"Sunlight pickup\nT+{pk['met']:.0f} s · {pk['alt_km']:.0f} km",
             color=st.OK, size=9, ha="right", va="center")

    # inset: whole powered corridor to insertion
    ix0, ix1, iy0, iy1 = inset_ext
    ih = 0.30 * ax_h
    iw = ih * (ix1 - ix0) * np.cos(np.radians((iy0 + iy1) / 2)) / (iy1 - iy0)
    ia = fig.add_axes([(side + 0.14) / W, (bot + ax_h - ih - 0.14) / H, iw / W, ih / H])
    base(ia, inset_ext, 10, 10)
    for t in ia.texts:
        t.set_fontsize(6.5)
    for s in ia.spines.values():
        s.set_edgecolor("#5a6470")
    draw_track(ia, 0.6)
    ia.plot(pad[1], pad[0], "*", ms=8, mfc=st.ACCENT, mec=(0, 0, 0, 0.8), mew=0.6, zorder=9)
    la_i, lo_i = np.interp(D["t_ins"], T, lat), np.interp(D["t_ins"], T, lon)
    ia.plot(lo_i, la_i, "o", ms=4, mfc=C_VN, mec="#101418", zorder=9)
    st.label(ia, lo_i + 1.2, la_i + 0.3, f"Insertion\nT+{D['t_ins']:.0f} s\n{sol['target']['alt_km']:.0f} km", size=7)
    ia.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="#ffffff", lw=0.8, zorder=8))

    handles = [Line2D([], [], color=C_SH, lw=2.6, label="Ascent, in shadow"),
               Line2D([], [], color=C_LIT, lw=3.4, label="Sunlit plume"),
               Line2D([], [], color=C_VN, lw=1.8, ls=(0, (4, 2.5)), label="Vernier climb"),
               Line2D([], [], marker="s", ls="", ms=7, mfc=st.OK, mec="k", label="Sunlight pickup"),
               Patch(facecolor=(1, 1, 1, 0.12), edgecolor=C_VIS, ls=(0, (5, 3)), label=f"{min_elev:.0f}° visibility"),
               Line2D([], [], color=st.SUN, lw=1.5, label="Sun 0°"),
               Line2D([], [], color=st.ACCENT, lw=1.3, label="Sun −12°"),
               Patch(facecolor=(1, 82 / 255, 82 / 255, 0.3), edgecolor=st.DANGER, label="NOTAM zone"),
               Line2D([], [], marker="x", ls="", ms=6, mew=1.6, color=st.DANGER, label="Stage 1 / fairing")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(side / W, 0.06 / H), ncol=(len(handles) + 1) // 2,
               frameon=False, fontsize=8.5, labelcolor=st.FG, handlelength=2.2, columnspacing=1.6,
               handletextpad=0.6, borderaxespad=0.0)
    fig.text(1 - side / W, 0.10 / H, "Illustration: Mickey · base map NASA Blue Marble", color=st.FGDIM,
             fontsize=7.5, ha="right", va="bottom")
    png = out_dir / f"{stem}_visibility.png"
    dpi = 3960 / max(W, H)
    fig.savefig(png, dpi=dpi, facecolor=st.BG)
    plt.close(fig)
    print(f"wrote {png.name}  {W:.2f} x {H:.2f} in at {dpi:.0f} dpi")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mission", default="tslc1118")
    ap.add_argument("--min-elev", type=float, default=10.0)
    ap.add_argument("--extent", type=float, nargs=4, default=[98.0, 120.0, 22.0, 45.5], metavar=("W", "E", "S", "N"))
    ap.add_argument("--inset", type=float, nargs=4, default=[96.0, 124.0, 1.0, 43.0])
    ap.add_argument("--no-map", action="store_true")
    a = ap.parse_args()
    from build_traj import MISSIONS, mission_dir
    out_dir, stem = mission_dir(a.mission), MISSIONS[a.mission]["stem"]
    sol, rows, D, out = analyse(out_dir, stem, a.min_elev)
    p, c = out["pickup_geometric"], out["pickup_clear_20km"]
    print(f"T0 {out['t0_utc']} UTC; Sun at the pad {out['sun_at_pad_T0']:+.2f} deg")
    print(f"sunlight pickup: geometric T+{p['met']:.0f} s at {p['alt_km']:.1f} km; clear (ray >= 20 km) T+{c['met']:.0f} s"
          f" at {c['alt_km']:.1f} km, {c['lat']:.3f} N {c['lon']:.3f} E")
    print(f"Sun-ray tangent at MECO-1 {out['sun_tangent_km']['MECO1']:+.1f} km, MECO-2 {out['sun_tangent_km']['MECO2']:+.1f} km;"
          f" ground Sun below the vehicle: MECO-1 {out['ground_sun_alt']['MECO1']:+.1f}, MECO-2 "
          f"{out['ground_sun_alt']['MECO2']:+.1f}, insertion {out['ground_sun_alt']['insertion']:+.1f} deg")
    print(f"reference epoch T+{out['ref_epoch_met']:.0f} s = {out['ref_epoch_utc']} UTC")
    print(f"{'city':13s} {'Sun':>6s} {'max el':>7s} {'az':>4s} {'km':>5s} {'T+':>5s} {'el@MECO-2':>10s} {'az':>4s} {'vernier el':>10s}")
    for q in sorted(out["cities"], key=lambda q: -q["max_el"]):
        print(f"{q['city']:13s} {q['sun_alt']:+6.1f} {q['max_el']:7.1f} {q['az']:4d} {q['range_km']:5d} {q['met_max_el']:5.0f}"
              f" {q['el_meco2']:10.1f} {q['az_meco2']:4d} {q['vernier_max_el']:10.1f}")
    json.dump(out, open(out_dir / f"{stem}_twilight.json", "w"), indent=1, ensure_ascii=False)
    print(f"wrote {stem}_twilight.json")
    if not a.no_map:
        render(out_dir, stem, sol, rows, D, out, a.min_elev, tuple(a.extent), tuple(a.inset))


if __name__ == "__main__":
    main()
