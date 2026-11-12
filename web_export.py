#!/usr/bin/env python3
"""Export a build_traj.py mission to the 3D page's mission data (launch-page template, PROFILE_FORMAT.md).

Resamples the built trajectory onto uniform grids; it fits and integrates nothing.  The Earth-fixed track holds for
every T0 on the page; only the Sun moves.  Sun and Moon vectors come from the shared DE440s kernel (~/JPLDE), read
in place.  Also writes <stem>_page_ref.json: the sunlight pickup and the 10 deg zone areas computed here with the
page's own definitions (jelly-grid.js, profile.js), for the page-vs-Python check.

Run after build_traj.py:  ~/.venvs/astro313/bin/python web_export.py [--mission tslc1118|y108|all] [--out FILE]
Default output: docs/data/<page id>.json (the page must exist; use --out for the first new_page.py run).
"""
import argparse
import csv
import json
import math
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from skyfield.api import Loader, load_file
from skyfield.framelib import itrs

from build_traj import MISSIONS, mission_dir   # before cz2d: cz2d puts sibling projects (with their own build_traj) first on sys.path
import cz2d as C
from nrol97 import ecef_to_geodetic

HERE = Path(__file__).resolve().parent
DOCS = HERE / "docs"
REPO = "https://github.com/exoplanet5/CZ-2D"
A_WGS, B_WGS = 6378.137, 6356.752314245
RAAN_RATE = 360.98564736629 / 1440.0          # deg of RAAN per minute of T0

# Rocket cartoon: 41.06 m, 3.35 m core.  Stage 1 27.9 m, stage 2 with the fairing overlap 6.2 m, fairing 7.0 m.
SEGMENTS = dict(s1=dict(len=27.9, w=3.35, kind="liquid", band=True), s2=dict(len=6.2, w=3.35, kind="liquid"),
                fairing=dict(len=7.0, w=3.35, kind="fairing"), payload=dict(len=3.2, w=2.2, kind="payload"))

PAGE = dict(
    tslc1118=dict(
        id="cz2d-tslc-1118", label="CZ-2D · Taiyuan · 18 Nov", title="CZ-2D Launch visibility", vehicle="CZ-2D",
        mission="Taiyuan · SSO 526 km", pad_name="Taiyuan LC-9", t0_min="22:43:00", t0_max="23:05:00",
        presets=[("Open 22:43", "22:43:00"), ("Y104 pattern 22:50", "22:50:00"), ("Model 22:51", "22:51:00")],
        met_default=240, orbit_label="526 km SSO", badge="Pre-launch model", home="Xi'an",
        about=["Pre-launch model fitted to the two NOTAM zones; T0 and the orbit are inferred, not announced.",
               "WGS84 geometry · JPL DE440s sunlight."],
        views=[["Global", 108, 20, 3.6], ["Pad", 111.59, 37.9, 1.12], ["Jelly zone", 109.5, 34, 1.5],
               ["Corridor", 107, 20, 2.3]],
        extent=[95, 123, -12, 47], step=0.125,
        cities=[["Xi'an", 108.94, 34.27], ["Yan'an", 109.49, 36.59], ["Taiyuan", 112.55, 37.87],
                ["Zhengzhou", 113.63, 34.75], ["Luoyang", 112.45, 34.62], ["Wuhan", 114.31, 30.59],
                ["Shijiazhuang", 114.51, 38.04], ["Chongqing", 106.55, 29.56], ["Changsha", 112.94, 28.23],
                ["Chengdu", 104.07, 30.66], ["Lanzhou", 103.83, 36.06], ["Yinchuan", 106.23, 38.49],
                ["Jinan", 117.12, 36.65], ["Baotou", 109.84, 40.66], ["Hohhot", 111.75, 40.84],
                ["Guiyang", 106.63, 26.65], ["Xining", 101.78, 36.62], ["Beijing", 116.40, 39.90],
                ["Yulin", 109.73, 38.29], ["Baoji", 107.24, 34.36], ["Yichang", 111.29, 30.69],
                ["Xiangyang", 112.12, 32.01], ["Tianshui", 105.72, 34.58], ["Datong", 113.30, 40.08]],
        extra_cities=[["Vientiane", 102.63, 17.97, 1], ["Da Nang", 108.22, 16.05, 1], ["Bangkok", 100.50, 13.76, 1],
                      ["Phnom Penh", 104.92, 11.56, 1], ["Ho Chi Minh", 106.63, 10.82, 1],
                      ["Kuala Lumpur", 101.69, 3.14, 1], ["Singapore", 103.82, 1.35, 1]],
        zone_desc={0: "Stage 1 zone, Lanzhou / Wuhan FIR; the same polygon as Y99 (2024-12-16) and Y104 (2026-09-19)",
                   1: "Fairing zone, Wuhan FIR; the modelled fairing lands 15 km west of its axis"},
        pickup_note="The vehicle reaches sunlight before staging, so the whole main-engine burn is lit over ground "
                    "still in twilight. Clear sunlight: the Sun's ray passes at least 20 km above the surface. "
                    "The vernier plume after main-engine cutoff has 6 % of the thrust and is far fainter.",
        viewers_title="Best views of the main-engine burn",
        viewers_note="Sunlit main-engine burn only, sea-level sites. Bold: Sun at or below −6°, the darker sky. "
                     "Dim: Sun already up.",
        notes=[
            "Pre-launch model, not telemetry and not an official ephemeris. T0 22:51:00 UTC (06:51 Beijing, 19 November) is assumed: stage-1 NOTAM start + 8 min, the pattern of CZ-2D Y99 and Y105 from Taiyuan (Y104 was + 7 min). One minute of T0 moves the RAAN by 0.25 deg.",
            "Mission identification is inferred. The stage-1 zone A5897/26 is corner for corner the zone of Y99 (PIESAT-2 09-12) and Y104 (PIESAT-2 13-16, 2026-218), and both went to 518 x 533 km at 97.49 deg. The model flies a 526 km circle at that inclination with a 1200 kg payload (four ~270 kg satellites plus dispenser, an estimate).",
            "Vehicle: published CZ-2D figures (CGWIC / SAST): stage 1 182 t of propellant, 2961.6 kN at sea level; stage 2 52.7 t, YF-22C main engine 742 kN plus YF-23C vernier 47.1 kN. The stage-2 dry mass of 4.5 t is calibrated to 1.9 t to 500 km SSO; the 700 km capacity is then within 50 kg of the published 1.3 t.",
            "Profile: stage 1 is fitted to the centre of its zone. Stage 2 burns main engine + vernier, then the vernier alone in a long climb to direct circular insertion; the main-engine cutoff is the propellant-optimal one. No CZ-2D event timeline is published, so MECO-2 and insertion times are model results.",
            "Plane: launch azimuth 187.66 deg gives 88.4 deg at MECO-1, so stage 2 flies a 9 deg dogleg to 97.49 deg. The RAAN carries +0.09 deg, the mean offset of Y99 and Y104 first TLEs from the fuel-optimal model plane.",
            "Fairing: jettison at the 1135 W/m2 heat-flux point. The modelled fairing lands 4 km long and 15 km west of the A5898/26 axis, inside the 72 x 35 km zone. The model starts its yaw at stage-2 ignition; the real ground track may run up to about 15 km east of it early in stage 2.",
            "Sunlight: geometric WGS84 occultation of the Sun's centre with DE440s, 30 s cadence. Near pickup the ray grazes the ground, so the plume starts deep red; 'clear sunlight' means the ray passes at least 20 km above the surface.",
            "On this page the zones' night edge follows the live Sun 0 line, so it moves during playback; the area figures use one reference epoch, the sunlight pickup. Glow is an illustration. Plume size and brightness, weather, extinction and terrain are not modeled, and the sky is in civil twilight under most of the track. Geometric visibility is not a claim of naked-eye detection.",
            "Shifting T0 keeps the Earth-fixed path and moves the Sun. Launch times, pad coordinates and satellite masses: data from GCAT (J. McDowell, planet4589.org/space/gcat). First TLEs of the earlier flights: Space-Track.org.",
        ]),
    y108=dict(
        id="cz2d-y108", label="CZ-2D Y108 · Jiuquan · 11 Nov", title="CZ-2D Y108 Launch visibility",
        vehicle="CZ-2D Y108", mission="Origins Pathfinder · SSO 600 km", pad_name="Jiuquan LC-43/94",
        t0_min="11:41:00", t0_max="11:41:59", presets=[("Official 11:41", "11:41:00")], met_default=360,
        orbit_label="600 km SSO", badge="Post-launch model, orbit assumed", home="Y108 observer",
        about=["Model of a flown launch: official lift-off minute, one NOTAM zone, and an assumed 600 km SSO. "
               "Not telemetry; no catalogued orbit was available when it was built.",
               "WGS84 geometry · JPL DE440s sunlight."],
        views=[["Global", 97, 22, 3.6], ["Pad", 100.30, 39.96, 1.12], ["Observer", 98.5, 36.8, 1.4],
               ["Corridor", 95, 22, 2.3]],
        extent=[70, 122, -14, 50], step=0.25,
        cities=[["Jiuquan", 98.49, 39.73], ["Jiayuguan", 98.29, 39.77], ["Zhangye", 100.45, 38.93],
                ["Dunhuang", 94.66, 40.14], ["Golmud", 94.90, 36.40], ["Delingha", 97.36, 37.37],
                ["Xining", 101.78, 36.62], ["Lanzhou", 103.83, 36.06], ["Yushu", 97.01, 33.00],
                ["Chamdo", 97.17, 31.14], ["Lhasa", 91.13, 29.65], ["Nyingchi", 94.36, 29.65],
                ["Chengdu", 104.07, 30.66], ["Kunming", 102.71, 25.04], ["Mandalay", 96.09, 21.96],
                ["Yangon", 96.17, 16.84], ["Dhaka", 90.41, 23.81]],
        extra_cities=[["Y108 observer", 98.616, 39.675, 1], ["Jinta", 98.90, 39.98, 3], ["Delingha", 97.36, 37.37, 2],
                      ["Thimphu", 89.64, 27.47, 2], ["Guwahati", 91.74, 26.14, 2], ["Dhaka", 90.41, 23.81, 1],
                      ["Mandalay", 96.09, 21.96, 2], ["Yangon", 96.17, 16.84, 1], ["Bangkok", 100.50, 13.76, 1],
                      ["Port Blair", 92.74, 11.62, 2]],
        zone_desc={0: "Stage 1 zone, Lanzhou FIR; the only NOTAM of this launch"},
        pickup_note="Dusk launch: the whole main-engine burn is in the Earth's shadow. The vehicle climbs into "
                    "sunlight at about T+6:00, in the vernier climb, which is when it was seen from near Jiuquan. "
                    "Clear sunlight: the Sun's ray passes at least 20 km above the surface.",
        viewers_title="Main-engine burn: in shadow",
        viewers_note="No city sees a sunlit main-engine burn. The green zone is where the sunlit vernier climb "
                     "rises 10° or more; open the sky chart at 'Y108 observer' for the T+6:00 sighting.",
        notes=[
            "Model of a flown launch, not telemetry and not an official ephemeris. Lift-off was announced as 19:41 Beijing time (11:41 UTC) on 11 November 2026; the seconds are unknown and 11:41:00 is used.",
            "The 600 km Sun-synchronous orbit is an assumption. The launch azimuth of 191.3 deg is 3.5-4 deg east of the standard CZ-2D SSO corridor from Jiuquan, and by itself points to an inclination near 94-95 deg; reaching SSO from it needs about 3.4 deg of extra dogleg. No catalogued orbit was available when the model was built.",
            "Vehicle: the release describes a lightened configuration (3.35 m composite fairing, lightened stage 2, no fins). Its mass savings are estimates (fairing -250 kg, stage 2 -150 kg, fins -300 kg). The payload is set to that vehicle's capacity on this corridor, 1811 kg.",
            "Published CZ-2D figures (CGWIC / SAST): stage 1 182 t of propellant, 2961.6 kN at sea level; stage 2 52.7 t, YF-22C main engine 742 kN plus YF-23C vernier 47.1 kN. The stage-2 dry mass of 4.5 t is calibrated to 1.9 t to 500 km SSO.",
            "Profile: stage 1 is fitted to the centre of NOTAM A5728/26. Stage 2 burns main engine + vernier, then the vernier alone in a long climb to direct circular insertion; the main-engine cutoff is the propellant-optimal one. No CZ-2D event timeline is published.",
            "Observer: a sighting at RA 21h33m32.7s, Dec -31d11m16s at T0+6 min back-projects to 39.675 N 98.616 E, 12 km ESE of Jiuquan, 755 km from the vehicle. A second of timing error moves that point about 4 km along a NNE-SSW line. It is marked 'Y108 observer'; the page's sky chart there is at sea level, the solution on terrain at 1449 m.",
            "Sunlight: geometric WGS84 occultation of the Sun's centre with DE440s, 30 s cadence. The vehicle leaves the Earth's shadow at about T+6:00; 'clear sunlight' means the ray passes at least 20 km above the surface.",
            "The zone's night edge follows the live Sun 0 line; the area figures use one reference epoch, the sunlight pickup. Glow is an illustration. Plume size and brightness, weather, extinction and terrain are not modeled.",
            "Launch times, pad coordinates and satellite masses of the validation flights: data from GCAT (J. McDowell, planet4589.org/space/gcat). Their first TLEs: Space-Track.org.",
        ]),
)


def load_xyz(path):
    rows = list(csv.DictReader(open(path)))
    return (np.array([float(r["met"]) for r in rows]),
            np.array([[float(r[k]) for k in ("x_ecef", "y_ecef", "z_ecef")] for r in rows]))


def uniform(spl, t_a, t_b, dt_max):
    """Uniform samples t_a..t_b ending exactly at t_b. -> track dict without style"""
    n = max(1, int(math.ceil((t_b - t_a) / dt_max - 1e-9)))
    ts = np.linspace(t_a, t_b, n + 1)
    g = np.array([ecef_to_geodetic(x) for x in spl(ts)])
    return dict(met0=float(t_a), dt=float((t_b - t_a) / n), lon=[round(float(v), 6) for v in g[:, 1]],
                lat=[round(float(v), 6) for v in g[:, 0]], alt=[round(max(float(v), 0.0), 4) for v in g[:, 2]])


def fall_track(rows, color, label):
    a = np.array(rows, float)                              # t, lat, lon, alt
    t = a[:, 0]
    n = max(1, int(math.ceil((t[-1] - t[0]) / 5.0)))
    ts = np.linspace(t[0], t[-1], n + 1)
    return dict(met0=float(t[0]), dt=float((t[-1] - t[0]) / n),
                lon=[round(float(x), 4) for x in np.interp(ts, t, a[:, 2])],
                lat=[round(float(x), 4) for x in np.interp(ts, t, a[:, 1])],
                alt=[round(max(float(x), 0.0), 3) for x in np.interp(ts, t, a[:, 3])], color=color, label=label,
                dashed=True)


def body_table(body, start, end, model):
    ts = Loader(str(HERE), verbose=False).timescale(builtin=True)
    eph = load_file(os.path.expanduser("~/JPLDE/de440s.bsp"))
    n = int((end - start).total_seconds() // 30) + 2
    utc = [start + timedelta(seconds=30 * k) for k in range(n)]
    vec = eph["earth"].at(ts.from_datetimes(utc)).observe(eph[body]).frame_xyz(itrs).km.T
    return dict(start_utc=start.strftime("%Y-%m-%dT%H:%M:%SZ"), dt=30,
                vectors=[[round(float(v), 1) for v in p] for p in vec], model=model)


# ---------------------------------------------------------------- the page's own jelly definitions, in Python
def clearance(p, sun):
    """profile.js sunlightClearance: > 0 when the segment vehicle -> Sun centre clears the WGS84 ellipsoid."""
    ax = np.array([A_WGS, A_WGS, B_WGS])
    r, d = p / ax, (sun - p) / ax
    t = np.clip(-(r @ d) / (d @ d), 0.0, 1.0)
    return float(((r + t * d) ** 2).sum() - 1.0)


def site(lon, lat):
    la, lo = np.radians(lat), np.radians(lon)
    n = A_WGS / np.sqrt(1 - (1 - (B_WGS / A_WGS) ** 2) * np.sin(la) ** 2)
    P = np.stack([n * np.cos(la) * np.cos(lo), n * np.cos(la) * np.sin(lo), n * (B_WGS / A_WGS) ** 2 * np.sin(la)], -1)
    U = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], -1)
    return P, U


def page_reference(track, sun, t0, meta):
    """Pickup and zone areas as computed by the page (jelly-grid.js) from the exported track and Sun table."""
    n = len(track["lon"])
    T = track["met0"] + track["dt"] * np.arange(n)
    P = np.array([C.Z.geodetic_to_ecef(la, lo, h) for la, lo, h in zip(track["lat"], track["lon"], track["alt"])])
    s0 = datetime.strptime(sun["start_utc"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
    SV = np.array(sun["vectors"])
    off = (t0 - s0).total_seconds()

    def sun_at(met):
        x = (off + met) / sun["dt"]
        i = int(np.clip(np.floor(x), 0, len(SV) - 2))
        return SV[i] + (SV[i + 1] - SV[i]) * (x - i)

    def pos(met):
        x = (met - T[0]) / track["dt"]
        i = int(np.clip(np.floor(x), 0, n - 2))
        return P[i] + (P[i + 1] - P[i]) * (x - i)

    lit = np.array([clearance(P[i], sun_at(T[i])) > 0 for i in range(n)])
    pickup = None
    if lit[0]:
        pickup = float(T[0])
    else:
        i = int(np.argmax(lit))
        if lit[i]:
            lo, hi = T[i - 1], T[i]
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if clearance(pos(mid), sun_at(mid)) > 0:
                    hi = mid
                else:
                    lo = mid
            pickup = float(0.5 * (lo + hi))
    out = dict(pickup_met=pickup)
    if pickup is None:
        return out
    g = ecef_to_geodetic(pos(pickup))
    out.update(pickup_alt_km=float(g[2]), pickup_lat=float(g[0]), pickup_lon=float(g[1]))
    v = meta["visibility"]
    lon0, lon1, lat0, lat1 = v["extent"]
    step = v["grid_step_deg"]
    nx, ny = round((lon1 - lon0) / step) + 1, round((lat1 - lat0) / step) + 1
    gl, gt = np.meshgrid(lon0 + step * np.arange(nx), lat0 + step * np.arange(ny))
    O, U = site(gl, gt)
    end = T[0] + track["dt"] * (n - 1)
    count = math.ceil((end - T[0]) / v["time_step_s"])
    solid, s4, c4 = [], [], 0
    for k in range(count + 1):
        met = T[0] + (end - T[0]) * k / count
        p = pos(met)
        if clearance(p, sun_at(met)) <= 0:
            continue
        if met <= meta["solid_end_met"]:
            solid.append(p)
        else:
            if c4 % 3 == 0:
                s4.append(p)
            c4 += 1

    def max_el(pts):
        best = np.full(gl.shape, -1.0)
        for p in pts:
            d = p - O
            np.maximum(best, (d * U).sum(-1) / np.linalg.norm(d, axis=-1), out=best)
        return np.degrees(np.arcsin(np.clip(best, -1, 1))) - 10.0

    a, b = max_el(solid), max_el(s4)
    d = sun_at(pickup) - O
    sun_el = np.degrees(np.arcsin((d * U).sum(-1) / np.linalg.norm(d, axis=-1)))
    cell = 6371.0 ** 2 * np.radians(step) ** 2 * np.cos(np.radians(gt))
    night = sun_el <= 0
    out.update(area_solid_km2=float(cell[night & (a > 0)].sum()),
               area_prime_km2=float(cell[(a > 0) & (sun_el <= meta["jelly"]["prime_sun_deg"])].sum()),
               area_s4_km2=float(cell[night & ~(a > 0) & (b > 0)].sum()), samples=len(solid) + len(s4))
    return out


def export(key, out_path):
    M, G = MISSIONS[key], PAGE[key]
    mdir, stem = mission_dir(key), M["stem"]
    sol = json.loads((mdir / f"{stem}_solution.json").read_text())
    ev = {e["event"]: e for e in sol["events"]}
    t_meco1, t_fair = ev["MECO-1 / hot staging"]["met"], sol["stage2"]["t_fair"]
    t_meco2, t_ins = sol["stage2"]["t_meco2"], sol["stage2"]["t_ins"]
    t_meco1 = sol["stage1"]["t_meco"]
    period, revs, t_end = (sol["orbit_track"][k] for k in ("period_s", "revolutions", "t_end"))
    t0 = datetime.fromisoformat(sol["t0_utc"]).replace(tzinfo=timezone.utc)
    alt = sol["target"]["alt_km"]

    ta, Xa = load_xyz(mdir / f"{stem}_traj.csv")
    to, Xo = load_xyz(mdir / f"{stem}_orbit.csv")
    ascent = uniform(CubicSpline(ta, Xa), 0.0, t_ins, 1.0)
    orbit = uniform(CubicSpline(to, Xo), t_ins, t_end, 5.0)
    for k in ("lon", "lat", "alt"):                        # adjacent segments share their endpoint exactly
        ascent[k][-1] = orbit[k][0]
    ascent.update(color="#d95926", label="Ascent")
    orbit.update(color="#c77dff", label="Orbit")
    falls = sol["falls"]
    tracks = dict(ascent=ascent, orbit=orbit, s1fall=fall_track(falls["S1"]["rows"], "#ff8a50", "Stage 1 fall"),
                  fairfall=fall_track(falls["Fairing"]["rows"], "#cfd6dd", "Fairing fall"))

    pad = tuple(sol["pad"])
    zspec = [("S1", "#ff8a50", falls["S1"]), ("Fairing", "#cfd6dd", falls["Fairing"])]
    zones = []
    for i, q in enumerate(sol["zones"]):
        role, color, imp = zspec[i]
        ring = [[round(lo, 5), round(la, 5)] for la, lo in q["ring"]]
        ring.append(ring[0])
        cen = C.zone_centroid(q["ring"])
        zones.append(dict(key=q["id"].split("/")[0], id=q["id"], role=role, color=color, t0=q["t0"], t1=q["t1"],
                          ring=ring, centroid=[round(cen[1], 4), round(cen[0], 4)],
                          dist_km=round(C.gc_dist(pad, cen), 1), bearing=round(C.bearing(pad, cen), 2),
                          desc=G["zone_desc"][i], impact=dict(met=round(imp["t"], 1), lon=round(imp["lon"], 4),
                                                              lat=round(imp["lat"], 4))))

    def e(met, label, major=False, track=None, short=None):
        o = dict(met=float(met), label=label, major=major, kind="ship")
        if track:
            o["track"] = track
        if short:
            o["short"] = short
        return o

    events = [e(0, "Liftoff", True), e(C.STD.t_rise, "Pitch-over"),
              e(t_meco1, "MECO-1, hot staging, stage 2 ignition", True, short="MECO-1 · staging"),
              e(t_fair, "Fairing jettison"),
              e(t_meco2, "Main-engine cutoff (MECO-2), vernier continues", True, short="MECO-2"),
              e(t_ins, "Vernier cutoff, orbit insertion", True, short="Insertion"),
              e(t_ins + period, "1 revolution"), e(t_end, "1.5 revolutions"),
              e(falls["S1"]["t"], "Stage 1 impact", True, "s1fall"),
              e(falls["Fairing"]["t"], "Fairing impact", False, "fairfall")]
    if key == "y108":
        events.append(e(360.0, "Sighting from near Jiuquan (T0 + 6 min)", False))
    events.sort(key=lambda q: q["met"])

    meta = dict(
        title=G["title"], vehicle=G["vehicle"], mission=G["mission"], date=sol["t0_utc"][:10],
        pad=[pad[1], pad[0]], pad_alt=M["site"].pad_h, pad_name=G["pad_name"],
        t0_default=sol["t0_utc"][11:19], t0_min=G["t0_min"], t0_max=G["t0_max"],
        t0_presets=[dict(label=a, time=b) for a, b in G["presets"]], met_default=G["met_default"],
        insertion_met=t_ins, period_s=period, orbit_revolutions=revs, inc=sol["orbit"]["i"],
        perigee_km=sol["orbit"]["hp"], apogee_km=sol["orbit"]["ha"], az=sol["stage1"]["az"],
        s4_cutoff_met=t_ins, solid_end_met=t_meco2, plume_names=["main engine", "vernier"],
        cutoff_label="Insertion", model_badge=G["badge"], about=G["about"],
        readme_url=f"{REPO}/blob/main/missions/{M['dir']}/README.md#assumptions-and-uncertainty",
        stack=dict(total_m=41.1, burn_end=t_ins, segments=SEGMENTS,
                   states=[dict(until=t_meco1 + C.STD.dt_sep, segs=["s1", "s2", "fairing"], zoom=1, burn="liquid"),
                           dict(until=t_fair, segs=["s2", "fairing"], zoom=1.9, burn="liquid"),
                           dict(until=t_ins, segs=["s2", "payload"], zoom=2.4, burn="liquid"),
                           dict(segs=["s2", "payload"], zoom=2.4, burn=None)]),
        observer_home=G["home"], extra_cities=G["extra_cities"], flight_segments=["ascent", "orbit"],
        powered_intervals=[[0.0, t_ins]],
        phases=[dict(until=0, label="Pre-launch"), dict(until=t_meco1, label="Stage 1"),
                dict(until=t_meco2, label="Stage 2 main engine + vernier"), dict(until=t_ins, label="Vernier climb"),
                dict(until=t_end + 1, label=G["orbit_label"])],
        views=G["views"], visibility=dict(extent=G["extent"], grid_step_deg=G["step"], time_step_s=2, levels=[10]),
        jelly=dict(prime_sun_deg=-6, clear_ray_km=20.0, high_km=60,
                   stages=[dict(label="Stage 1", a=0.0, b=t_meco1), dict(label="Main engine", a=t_meco1, b=t_meco2),
                           dict(label="Vernier", a=t_meco2, b=t_ins)],
                   cities=G["cities"], solid_label="main-engine burn", pickup_note=G["pickup_note"],
                   viewers_title=G["viewers_title"], viewers_note=G["viewers_note"]),
        links=[dict(label="Model, validation and source code", href=REPO),
               dict(label="Mission notes and map", href=f"{REPO}/tree/main/missions/{M['dir']}")],
        notes=G["notes"])

    day = t0.replace(hour=0, minute=0, second=0)
    hms = lambda s: timedelta(hours=int(s[:2]), minutes=int(s[3:5]), seconds=int(s[6:8]))
    t_a = day + hms(G["t0_min"]) - timedelta(minutes=11)
    t_b = day + hms(G["t0_max"]) + timedelta(seconds=t_end + 90)
    sun = body_table("sun", t_a, t_b, "DE440s geometric Sun centre; WGS84 occultation")
    moon = body_table("moon", t_a, t_b, "DE440s geometric Moon centre")
    tle = sol["tle"]
    data = dict(schema_version=1, meta=meta, tracks=tracks, zones=zones, events=events,
                tle_ref=dict(t0=sol["t0_utc"] + "Z", raan_rate_deg_per_min=RAAN_RATE,
                             sets=[dict(name=f"{G['vehicle'].upper()} MODEL - NOT OBSERVED", l1=tle[1], l2=tle[2])]),
                sun=sun, moon=moon)
    out_path.write_text(json.dumps(data, separators=(",", ":"), ensure_ascii=False))
    ref = page_reference(ascent, sun, t0, meta)
    (mdir / f"{stem}_page_ref.json").write_text(json.dumps(ref, indent=1))
    print(f"wrote {out_path} ({out_path.stat().st_size / 1e3:.0f} kB): "
          + ", ".join(f"{k} {len(v['lon'])} pts dt {v['dt']:.2f}" for k, v in tracks.items()))
    print(f"   page reference: pickup T+{ref['pickup_met']:.2f} s"
          + (f" at {ref['pickup_alt_km']:.1f} km; 10 deg zone main {ref['area_solid_km2'] / 1e6:.3f} M km2, prime "
             f"{ref['area_prime_km2'] / 1e6:.3f}, vernier-only {ref['area_s4_km2'] / 1e6:.3f}" if "area_solid_km2" in ref else ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mission", choices=list(MISSIONS) + ["all"], default="all")
    ap.add_argument("--out", type=Path, help="output file (one mission only); default docs/data/<page id>.json")
    a = ap.parse_args()
    keys = list(MISSIONS) if a.mission == "all" else [a.mission]
    for key in keys:
        C.use_site(MISSIONS[key]["site"])
        export(key, a.out if a.out and len(keys) == 1 else DOCS / "data" / f"{PAGE[key]['id']}.json")


if __name__ == "__main__":
    main()
