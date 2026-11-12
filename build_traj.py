"""CZ-2D ascent + insertion orbit from the NOTAM drop zones, for the missions in MISSIONS.

  y108      CZ-2D Y108 (探源计划探路者卫星), JSLC LC-43/94, 2026-11-11 19:41 BJT: 600 km SSO, payload = the
            lightened vehicle's capacity on the corridor (the release says the rocket was lightened "to meet the
            satellite's lift requirement").  One zone (stage 1).
  tslc1118  CZ-2D from TSLC LC-9, 2026-11-18 ~22:51 UTC (NOTAM A5897/26 stage 1, A5898/26 fairing): taken as the
            follow-up of 2026-218 (PIESAT-2 13-16), whose zones are identical; target = the 2026-218 / 2024-241
            insertion orbit, 526 km circular, i 97.49.  T0 = stage-1 NOTAM start + 8 min (Y99, Y105).

Stage 1 is fitted to the near zone; stage 2 flies main engine + vernier, then a vernier-only climb to direct
circular insertion, with the propellant-optimal main-engine cutoff and the fuel-optimal dogleg to the target
plane (cz2d.py, validated in validate_precedents.py).  A second zone is the fairing check.

Run: ~/.venvs/astro313/bin/python build_traj.py [--mission y108|tslc1118] [--t0 ISO] [--alt KM] [--inc I]
     [--payload KG] [--stem STEM]
Reads the NOTAM file from, and writes into, missions/<dir>/:
<stem>_{solution.json,traj.csv,orbit.csv,traj.kml,SSO<alt>.tle}
"""
import argparse
import csv
import json
from datetime import datetime

import numpy as np
import simplekml

import cz2d as C
from nrol97 import RE, coast, fit_tle, v_rel

MISSIONS = dict(
    y108=dict(dir="Y108_JSLC_20261111", stem="CZ-2D_Y108", name="CZ-2D Y108, 探源计划探路者卫星 (CAS Origins Pathfinder)", site=C.JSLC,
              notam="notam_jslc_sso", t0="2026-11-11T11:41:00", alt=600.0, inc=None, veh=C.Y108, payload=None,
              raan_offset=0.0, marks=[("T0+6 min observation", 360.0)], tle=(99108, "26999A", "CZ-2D Y108 PAYLOAD")),
    # raan_offset: the two flights on these exact zones (Y99, Y104; T0 known to the minute) came out 0.09 deg
    # east of the fuel-optimal model plane (validation.json), so the predicted plane is shifted by that much.
    tslc1118=dict(dir="TSLC_20261118", stem="CZ-2D_TSLC_20261118", name="CZ-2D, TSLC LC-9, 2026-11-18 (PIESAT-2 follow-up to 2026-218)",
                  site=C.TSLC, notam="notam_tslc_sso", t0="2026-11-18T22:51:00", alt=526.0, inc=97.49, veh=C.STD,
                  payload=1200.0, raan_offset=0.09, marks=[], tle=(99118, "26998A", "CZ-2D TSLC 1118 PAYLOAD")),
)


def mission_dir(key):
    """Folder holding a mission's NOTAM file and every output: missions/<dir>/."""
    return C.HERE / "missions" / MISSIONS[key]["dir"]


def phase_of(t, s1, s2):
    if t <= C.STD.t_rise:
        return "vertical rise"
    if t <= s1["t_meco"]:
        return "stage 1"
    if t <= s2["t_meco2"]:
        return "stage 2 main+vernier"
    if t <= s2["t_ins"] + 1e-6:
        return "stage 2 vernier"
    return "orbit"


def thrust_of(phase, veh):
    return {"vertical rise": veh.s1_T_sl, "stage 1": veh.s1_T_vac, "stage 2 main+vernier": veh.T_main + veh.T_vern,
            "stage 2 vernier": veh.T_vern}.get(phase, 0.0)


def state(fr, s1, s2, orb, t):
    a = s1["arc"] if t <= s1["t_meco"] else (s2["arc"] if t <= s2["t_ins"] else orb)
    y = a(t)
    r, v = y[:3], y[3:6]
    la, lo, h = fr.lla(r, t)
    vr = v_rel(r, v)
    fpa = np.degrees(np.arcsin(r @ vr / np.linalg.norm(r) / np.linalg.norm(vr)))
    re, ve = fr.to_ecef(r, t, v)
    m = y[6] if len(y) > 6 else np.nan
    return dict(t=t, la=la, lo=lo, h=h, v=np.linalg.norm(v), vrel=np.linalg.norm(vr), fpa=fpa, re=re, ve=ve, m=m,
                r=r, vi=v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mission", choices=list(MISSIONS), default="y108")
    ap.add_argument("--t0", default=None)
    ap.add_argument("--alt", type=float, default=None)
    ap.add_argument("--inc", type=float, default=None, help="default: mission value, else sun-synchronous for --alt")
    ap.add_argument("--payload", type=float, default=None, help="default: mission value, else the vehicle capacity")
    ap.add_argument("--stem", default=None)
    a = ap.parse_args()
    M = MISSIONS[a.mission]
    a.t0, a.alt, a.stem = a.t0 or M["t0"], a.alt or M["alt"], a.stem or M["stem"]
    a.inc = M["inc"] if a.inc is None else a.inc
    a.payload = M["payload"] if a.payload is None else a.payload
    site = C.use_site(M["site"])
    out = mission_dir(a.mission)

    fr = C.Frame(datetime.fromisoformat(a.t0))
    zones = sorted(C.parse_notams(out / M["notam"]), key=lambda q: C.gc_dist(C.PAD, C.zone_centroid(q["ring"])))
    z = zones[0]
    ax = C.zone_axis_from_ring(z["ring"])
    zc = C.zone_centroid(z["ring"])
    inc = C.sso_inclination(a.alt) if a.inc is None else a.inc
    veh = M["veh"]
    print(f"{M['name']}; {site.name} {site.pad}, T0 {a.t0} UTC")
    print(f"NOTAM {z['id']} {z['b']}-{z['c']}: centroid {zc[0]:.3f} N {zc[1]:.3f} E, "
          f"{C.gc_dist(C.PAD, zc):.0f} km at {C.bearing(C.PAD, zc):.2f} deg, box {ax.length:.0f} x {ax.width:.0f} km")
    print(f"target {a.alt:.0f} km circular, i {inc:.3f} deg; vehicle {veh.name}, lift-off "
          f"{veh.m0(0) / 1e3:.1f} t + payload")

    caps = {}
    for name, v in [("standard", C.STD)] + ([] if veh is C.STD else [("lightened", veh)]):
        caps[name], _ = C.capacity(fr, v, ax, inc, a.alt, tol=5.0)
        print(f"   capacity on this corridor, {name} CZ-2D: {caps[name]:.0f} kg")
    payload = caps[name] if a.payload is None else a.payload
    print(f"   payload {payload:.0f} kg")
    R = C.fly(fr, veh, payload, ax, inc, a.alt, raan_offset=M["raan_offset"], verbose=True)
    s1, s2 = R["s1"], R["s2"]
    t_ins = s2["t_ins"]
    orb = coast(s2["arc"].r_end, s2["arc"].v_end, t_ins, t_ins + 2 * 5900.0, dt=5.0)
    el = C.elements(s2["arc"].r_end, s2["arc"].v_end)
    t1, la1, lo1, f1 = C.s1_impact(fr, s1)
    _, laf, lof, ff = C.fairing_impact(fr, s2, veh)
    al1, cr1 = ax.along_cross(la1, lo1)

    def fall_rows(arc, step=5.0):
        ts = np.append(np.arange(arc.T[0], arc.t_end, step), arc.t_end)
        return [[round(float(t), 2)] + [round(float(v), 4) for v in fr.lla(arc(t)[:3], t)] for t in ts]

    falls = dict(S1=dict(t=t1, lat=la1, lon=lo1, rows=fall_rows(f1)),
                 Fairing=dict(t=ff.t_end, lat=laf, lon=lof, rows=fall_rows(ff)))
    period = 2 * np.pi * np.sqrt((RE + (el["hp"] + el["ha"]) / 2) ** 3 / C.MU)
    revs = 1.5
    t_end = t_ins + revs * period
    n0 = C.elements(s1["arc"].r_end, s1["arc"].v_end)["i"]

    events = [("lift-off", 0.0), ("pitch-over", veh.t_rise), ("MECO-1 / hot staging", s1["t_meco"]),
              ("fairing jettison", s2["t_fair"]), ("main-engine cutoff (MECO-2)", s2["t_meco2"]),
              *M["marks"], ("vernier cutoff / insertion", t_ins)]
    print(f"\nstage 1: az {s1['az']:.3f} deg, kick {s1['kick']:.3f} deg; MECO-plane i {n0:.2f} deg; impact "
          f"{la1:.3f} N {lo1:.3f} E (zone axis {al1:+.1f} / {cr1:+.1f} km) T+{t1:.0f} s")
    print(f"fairing impact {laf:.3f} N {lof:.3f} E, {C.gc_dist(C.PAD, (laf, lof)):.0f} km, T+{ff.t_end:.0f} s")
    fair_zone = None
    if len(zones) > 1:
        zf = zones[1]
        alf, crf = C.zone_axis_from_ring(zf["ring"]).along_cross(laf, lof)
        cf = C.zone_centroid(zf["ring"])
        fair_zone = dict(id=zf["id"], centroid=cf, range_km=C.gc_dist(C.PAD, cf), along_cross_km=[alf, crf])
        print(f"   fairing zone {zf['id']} {zf['b']}-{zf['c']}: centroid {cf[0]:.3f} N {cf[1]:.3f} E, "
              f"{fair_zone['range_km']:.0f} km; modelled fairing {alf:+.1f} km along / {crf:+.1f} km cross")
    ev_out = []
    for name, t in events:
        S = state(fr, s1, s2, orb, t)
        ev_out.append(dict(event=name, met=round(t, 1), utc=fr.utc(t).isoformat(timespec="seconds"),
                           lat=round(S["la"], 4), lon=round(S["lo"], 4), alt_km=round(S["h"], 1),
                           v_inertial=round(S["v"], 3), v_rel=round(S["vrel"], 3), fpa=round(S["fpa"], 2),
                           mass_kg=None if np.isnan(S["m"]) else round(S["m"]),
                           downrange_km=round(C.gc_dist(C.PAD, (S["la"], S["lo"])), 0)))
        print(f"   T+{t:6.1f} {fr.utc(t):%H:%M:%S}  {name:30s} {S['la']:7.3f} N {S['lo']:8.3f} E  h {S['h']:6.1f} km"
              f"  v {S['v']:.3f} km/s  fpa {S['fpa']:+5.1f}  {ev_out[-1]['downrange_km']:5.0f} km")
    print(f"insertion orbit (osculating): {el['hp']:.1f} x {el['ha']:.1f} km, i {el['i']:.3f}, RAAN {el['raan']:.3f}, "
          f"argp {el['argp']:.1f}, u {el['u']:.1f}")
    print(f"node: RAAN {R['raan']:.3f} (fuel-optimal {M['raan_offset']:+.2f} deg) vs plane through MECO-1 "
          f"{R['raan_through_meco']:.3f}; "
          f"main engine {s2['t_main']:.1f} s, vernier {t_ins - s2['t_ses']:.1f} s; margin {s2['margin']:.0f} kg")

    lines, rms, mx = fit_tle(fr, orb, t_ins + 60.0, t_ins + 2 * 5800.0, *M["tle"])
    (out / f"{a.stem}_SSO{a.alt:.0f}.tle").write_text("\n".join(lines) + "\n")
    print(f"TLE fit rms {rms:.2f} km (max {mx:.2f})\n" + "\n".join(lines))

    rows = []
    for t in np.arange(0.0, t_ins + 600.0, 1.0):
        S = state(fr, s1, s2, orb, t)
        ph = phase_of(t, s1, s2)
        rows.append(dict(met=round(t, 1), utc=fr.utc(t).isoformat(timespec="seconds"), phase=ph,
                         lat=round(S["la"], 5), lon=round(S["lo"], 5), alt_km=round(S["h"], 3),
                         x_ecef=round(S["re"][0], 4), y_ecef=round(S["re"][1], 4), z_ecef=round(S["re"][2], 4),
                         vx_ecef=round(S["ve"][0], 6), vy_ecef=round(S["ve"][1], 6), vz_ecef=round(S["ve"][2], 6),
                         v_inertial=round(S["v"], 4), v_rel=round(S["vrel"], 4), fpa=round(S["fpa"], 3),
                         mass_kg="" if np.isnan(S["m"]) else round(S["m"], 1),
                         thrust_kN=round(thrust_of(ph, veh) / 1e3, 1)))
    orows = []
    for t in np.append(np.arange(t_ins, t_end, 5.0), t_end):
        y = orb(t)
        la, lo, h = fr.lla(y[:3], t)
        re = fr.to_ecef(y[:3], t)
        orows.append(dict(met=round(float(t), 3), lat=round(la, 5), lon=round(lo, 5), alt_km=round(h, 3),
                          x_ecef=round(re[0], 4), y_ecef=round(re[1], 4), z_ecef=round(re[2], 4)))
    with open(out / f"{a.stem}_orbit.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(orows[0]))
        w.writeheader()
        w.writerows(orows)
    with open(out / f"{a.stem}_traj.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    sol = dict(mission=M["name"], site=site.name, t0_utc=a.t0, notam=z, zones=zones, fairing_zone=fair_zone,
               raan_offset=M["raan_offset"], pad=C.PAD, target=dict(alt_km=a.alt, inc=inc), vehicle={k: getattr(veh, k) for k in veh.__dataclass_fields__},
               capacity_kg=caps, payload_kg=round(payload), stage1=dict(az=s1["az"], kick=s1["kick"], t_meco=s1["t_meco"],
               meco_plane_i=n0, impact=[la1, lo1], impact_met=t1, zone_along_cross_km=[al1, cr1]),
               fairing_impact=[laf, lof], stage2=dict(t_main=s2["t_main"], t_meco2=s2["t_meco2"], t_ins=t_ins,
               steering=s2["x"].tolist(), margin_kg=s2["margin"], t_fair=s2["t_fair"]),
               raan=R["raan"], raan_through_meco=R["raan_through_meco"],
               orbit=dict(hp=el["hp"], ha=el["ha"], i=el["i"], raan=el["raan"], argp=el["argp"]), tle=lines,
               tle_rms_km=rms, events=ev_out, falls=falls,
               orbit_track=dict(period_s=period, revolutions=revs, t_end=t_end))
    json.dump(sol, open(out / f"{a.stem}_solution.json", "w"), indent=1, ensure_ascii=False, default=float)

    k = simplekml.Kml(name=f"{a.stem.replace('_', ' ')} SSO {a.alt:.0f} km  T0 {fr.t0:%Y-%m-%d %H:%M} UTC")
    for q in zones:
        pg = k.newpolygon(name=q["id"], outerboundaryis=[(p[1], p[0]) for p in q["ring"] + [q["ring"][0]]])
        pg.style.polystyle.color = simplekml.Color.changealphaint(80, simplekml.Color.rgb(117, 117, 117))
    k.newpoint(name=site.name, coords=[(C.PAD[1], C.PAD[0])])
    for name, pts, col in [("stage 1", [r for r in rows if r["phase"] in ("vertical rise", "stage 1")], (255, 184, 79)),
                           ("stage 2 main+vernier", [r for r in rows if r["phase"] == "stage 2 main+vernier"], (255, 213, 79)),
                           ("stage 2 vernier", [r for r in rows if r["phase"] == "stage 2 vernier"], (79, 195, 247)),
                           ("orbit", [r for r in rows if r["phase"] == "orbit"], (200, 200, 200))]:
        ls = k.newlinestring(name=name, coords=[(r["lon"], r["lat"], r["alt_km"] * 1e3) for r in pts])
        ls.altitudemode = simplekml.AltitudeMode.absolute
        ls.style.linestyle.color = simplekml.Color.rgb(*col)
        ls.style.linestyle.width = 3
    for e in ev_out:
        p = k.newpoint(name=f"{e['event']} T+{e['met']:.0f}", coords=[(e["lon"], e["lat"], e["alt_km"] * 1e3)])
        p.altitudemode = simplekml.AltitudeMode.absolute
    k.newpoint(name=f"stage 1 impact T+{t1:.0f}", coords=[(lo1, la1)])
    k.newpoint(name=f"fairing impact T+{ff.t_end:.0f}", coords=[(lof, laf)])
    k.save(str(out / f"{a.stem}_traj.kml"))
    print(f"wrote missions/{M['dir']}/{a.stem}_solution.json, _traj.csv, _orbit.csv, _traj.kml, _SSO{a.alt:.0f}.tle")


if __name__ == "__main__":
    main()
