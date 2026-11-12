"""Validate the CZ-2D model on flown JSLC and TSLC launches before trusting it for a new flight.

Tests (stage-1 zone hits are fit targets and validate nothing):
  1. capacity: published 1.9 t to 500 km SSO (Xinhua 2025-09-29) and 1.3 t to 645-700 km SSO (CGWIC / CASC),
     flown on the standard SSO corridor (Y97 zone).  s2_dry is the one calibrated knob (500 km point); the
     700 km point is the check.
  2. plane: stage 1 fitted to each flight's NOTAM zone, stage 2 flown to the catalogued altitude and inclination
     with the fuel-optimal node; predicted RAAN (J2-drifted to the first TLE epoch) vs Space-Track.
     T0 is known to the minute only (GCAT), i.e. +/-0.125 deg of RAAN.
  3. mass closure: propellant margin with the catalogued payload must be >= 0.

TSLC flights (site="TSLC") carry both NOTAM zones, so the modelled fairing impact is checked against the far zone;
Y99 and Y104 flew the exact zones of the 2026-11-18 launch.

Run: ~/.venvs/astro313/bin/python validate_precedents.py [--quick] [--site JSLC|TSLC]   -> validation.json
(--site re-runs one site and keeps the other site's flights already in validation.json)
"""
import argparse
import json
from datetime import datetime

import numpy as np

import cz2d as C
from nrol97 import RE, j2_node_rate

ST = json.load(open(C.HERE / "precedents" / "st_earliest.json"))
for _f in sorted((C.HERE / "precedents").glob("st_20??-???.json")):     # st.py launch --intldes ... --json
    for _rows in json.load(open(_f))["gp"].values():
        ST += _rows

# zone rings from ~/launch/notam/notam_data.json (stage-1 zone = the near one); payload masses from GCAT
PRECEDENTS = [
    dict(tag="Y97 2024-04-15 Siwei Gaojing 3-01", t0=datetime(2024, 4, 15, 4, 12), norad=59510, payload=600.0,
         ring=[(35.716667, 98.633333), (35.766667, 98.4), (35.433333, 98.316667), (35.383333, 98.533333)]),
    dict(tag="Y101 2025-01-17 PRSC-EO1 +2", t0=datetime(2025, 1, 17, 4, 7), norad=62728, payload=666.0,
         ring=[(35.983333, 98.366667), (35.616667, 98.25), (35.566667, 98.483333), (35.933333, 98.6)]),
    dict(tag="Y100 2025-03-15 Superview 3-02 +1", t0=datetime(2025, 3, 15, 4, 11), norad=63208, payload=1100.0,
         ring=[(35.75, 98.4), (35.383333, 98.3), (35.35, 98.516667), (35.716667, 98.633333)]),
    dict(tag="Y42 2025-06-14 Zhangheng 1-02", t0=datetime(2025, 6, 14, 7, 56), norad=64393, payload=730.0,
         ring=[(35.633333, 98.283333), (35.266667, 98.166667), (35.216667, 98.4), (35.583333, 98.516667)]),
    dict(tag="Y41 2021-12-29 Tianhui-4 x2", t0=datetime(2021, 12, 29, 11, 13), norad=50572, payload=2000.0,
         ring=[(35.203889, 99.67), (35.563056, 99.713056), (35.581111, 99.493056), (35.221944, 99.450833)],
         far=[(32.628056, 99.445833), (33.256111, 99.52), (33.286944, 99.146944), (32.658889, 99.075)]),
    # Taiyuan LC-9.  Payload = GCAT satellite masses (dispenser not counted)
    dict(tag="Y99 2024-12-16 PIESAT-2 09-12", site="TSLC", t0=datetime(2024, 12, 16, 18, 50), norad=62334,
         payload=1130.0,
         ring=[(33.566667, 110.516667), (33.533333, 110.733333), (33.166667, 110.683333), (33.183333, 110.45)],
         far=[(30.583333, 109.95), (30.55, 110.316667), (29.9, 110.216667), (29.95, 109.85)]),
    dict(tag="Y105 2026-03-25 Superview Neo 2-05/06", site="TSLC", t0=datetime(2026, 3, 25, 22, 51), norad=68377,
         payload=1200.0,
         ring=[(34.033333, 110.533333), (34.0, 110.75), (33.633333, 110.683333), (33.65, 110.466667)],
         far=[(30.166667, 109.633333), (30.116667, 110.0), (29.466667, 109.883333), (29.516667, 109.516667)]),
    dict(tag="Y104 2026-09-19 PIESAT-2 13-16", site="TSLC", t0=datetime(2026, 9, 19, 10, 50), norad=100749,
         payload=1080.0,
         ring=[(33.566667, 110.516667), (33.533333, 110.733333), (33.166667, 110.683333), (33.183333, 110.45)],
         far=[(30.583333, 109.95), (30.55, 110.316667), (29.9, 110.216667), (29.95, 109.85)]),
]


def first_tle(norad):
    rows = [r for r in ST if int(r["NORAD_CAT_ID"]) == norad]
    rows.sort(key=lambda r: r["EPOCH"])
    return rows[0]


def plane_test(p, veh):
    C.use_site(getattr(C, p.get("site", "JSLC")))
    tle = first_tle(p["norad"])
    inc = float(tle["INCLINATION"])
    h = (float(tle["PERIAPSIS"]) + float(tle["APOAPSIS"])) / 2
    fr = C.Frame(p["t0"])
    ax = C.zone_axis_from_ring(p["ring"])
    R = C.fly(fr, veh, p["payload"], ax, inc, h)
    s1, s2 = R["s1"], R["s2"]
    t_ins = s2["t_ins"]
    el = C.elements(s2["arc"].r_end, s2["arc"].v_end)
    ep = datetime.fromisoformat(tle["EPOCH"])
    dt_days = (ep - fr.utc(t_ins)).total_seconds() / 86400
    rate = j2_node_rate(RE + h, 0.0, inc)
    raan_pred = (el["raan"] + rate * dt_days) % 360
    d_raan = (raan_pred - float(tle["RA_OF_ASC_NODE"]) + 180) % 360 - 180
    # natural plane: no dogleg = plane through the stage-1 burnout point and its inertial velocity
    nat = C.elements(s1["arc"].r_end, s1["arc"].v_end)["i"]
    out = dict(tag=p["tag"], site=C.SITE.name, inc=inc, h=round(h, 1), payload=p["payload"], az=round(s1["az"], 2),
               kick=round(s1["kick"], 3), meco_plane_i=round(nat, 2), raan_pred=round(raan_pred, 3),
               raan_tle=float(tle["RA_OF_ASC_NODE"]), d_raan=round(d_raan, 3), margin_kg=round(s2["margin"]),
               t_meco1=round(s1["t_meco"], 1), t_meco2=round(s2["t_meco2"], 1), t_ins=round(t_ins, 1),
               yaw_tan=[round(float(v), 4) for v in s2["x"][2:4]])
    if "far" in p:                                     # fairing impact vs the far zone
        t, la, lo, _ = C.fairing_impact(fr, s2, veh)
        cf = C.zone_centroid(p["far"])
        al, cr = C.zone_axis_from_ring(p["far"]).along_cross(la, lo)
        out["fairing_along_cross_km"] = [round(al, 1), round(cr, 1)]
        out["fairing_impact"] = [round(la, 3), round(lo, 3)]
        out["fairing_vs_far_zone_km"] = round(C.gc_dist((la, lo), cf), 1)
        out["fairing_range_km"] = round(C.gc_dist(C.PAD, (la, lo)), 1)
        out["far_zone_range_km"] = round(C.gc_dist(C.PAD, cf), 1)
    return out


def capacity_test(veh):
    p = PRECEDENTS[0]
    C.use_site(C.JSLC)
    fr = C.Frame(p["t0"])
    ax = C.zone_axis_from_ring(p["ring"])
    out = {}
    for h, pub in [(500.0, 1900.0), (645.0, 1300.0), (700.0, 1300.0)]:
        cap, R = C.capacity(fr, veh, ax, C.sso_inclination(h), h, tol=10.0)
        out[f"{h:.0f}"] = dict(model_kg=round(cap), published_kg=pub, t_meco2=round(R["s2"]["t_meco2"], 1),
                               t_ins=round(R["s2"]["t_ins"], 1))
        print(f"   capacity {h:.0f} km SSO: model {cap:.0f} kg, published {pub:.0f} kg "
              f"(MECO-2 T+{R['s2']['t_meco2']:.0f} s, insertion T+{R['s2']['t_ins']:.0f} s)", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip the capacity runs")
    ap.add_argument("--site", choices=["JSLC", "TSLC"], help="re-run one site only")
    a = ap.parse_args()
    veh = C.STD
    out_path = C.HERE / "validation.json"
    old = json.load(open(out_path)) if out_path.exists() else {}
    res = dict(vehicle={k: getattr(veh, k) for k in veh.__dataclass_fields__})
    if a.quick or a.site == "TSLC":
        if "capacity" in old:
            res["capacity"] = old["capacity"]
    else:
        print("Capacity (Y97 corridor, standard CZ-2D):")
        res["capacity"] = capacity_test(veh)
    print("Plane / closure on flown launches:")
    todo = [p for p in PRECEDENTS if a.site in (None, p.get("site", "JSLC"))]
    tags = {p["tag"] for p in todo}
    res["flights"] = [f for f in old.get("flights", []) if f["tag"] not in tags] if a.site else []
    for p in todo:
        o = plane_test(p, veh)
        res["flights"].append(o)
        extra = (f"  fairing {o['fairing_range_km']:.0f} km vs far zone {o['far_zone_range_km']:.0f} km "
                 f"(miss {o['fairing_vs_far_zone_km']:.0f} km, along/cross {o['fairing_along_cross_km']})"
                 if "fairing_impact" in o else "")
        print(f"   {o['tag']:36s} i {o['inc']:6.2f} h {o['h']:5.0f}  az {o['az']:6.2f}  MECO-plane i {o['meco_plane_i']:6.2f}"
              f"  dRAAN {o['d_raan']:+6.3f}  margin {o['margin_kg']:+5d} kg  MECO-2 T+{o['t_meco2']:.0f}"
              f"  ins T+{o['t_ins']:.0f}{extra}", flush=True)
    json.dump(res, open(out_path, "w"), indent=1)
    print("wrote validation.json")


if __name__ == "__main__":
    main()
