"""Where was the ground observer?  Back-project an observed (RA, Dec) of the powered CZ-2D Y108 onto the terrain.

Given: at about T0+6 min (MET 360 s, +/- tens of seconds) the rocket was seen at RA 21h33m32.7s, Dec -31d11m16s.
For each candidate MET t the observer O satisfies O = r(t) - rho * u(t) with O on the ground, where r(t) is the
modelled rocket position (ECEF, CZ-2D_Y108_traj.csv) and u(t) is the observed direction rotated ICRS -> ITRS at
UTC = T0 + t.  The reported RA/Dec is taken as ICRS (J2000) astrometric, as a plate solve against a star catalogue
gives it: astropy's ICRS -> ITRS for a direction applies the same annual aberration the reference stars carry, and
refraction cancels against the stars.  The near-side intersection with WGS84 + ETOPO1 terrain (NOAA ERDDAP, cached
in terrain_etopo1.npz; orthometric heights, geoid ~ -35 m ignored: < 0.1 km in position) is the observer.

Varying t traces a 1-D locus; the model sensitivity runs (--sens) re-fly the ascent with other assumptions.

Run: ~/.venvs/astro313/bin/python observer_locate.py [--ra 21h33m32.7s --dec -31d11m16s] [--met 360 --dt 60]
                                                        [--of-date] [--sens]
Writes observer_locus.csv, observer_solution.json
"""
import argparse
import csv
import json
import urllib.request
from datetime import datetime, timedelta

import astropy.units as u
import numpy as np
from astropy.coordinates import FK5, ITRS, AltAz, EarthLocation, SkyCoord, get_sun
from astropy.time import Time
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

import cz2d as C
from nrol97 import ecef_to_geodetic, geodetic_to_ecef

OUT = C.HERE / "missions" / "Y108_JSLC_20261111"
TRAJ = OUT / "CZ-2D_Y108_traj.csv"
GRID = C.HERE / "terrain_etopo1.npz"
T0 = datetime(2026, 11, 11, 11, 41, 0)
TOWNS = [("Jiuquan", 39.744, 98.494), ("Jiayuguan", 39.773, 98.289), ("Jinta", 39.983, 98.902),
         ("Yumen (Laojunmiao)", 39.830, 97.550), ("Gaotai", 39.378, 99.819), ("Zhangye", 38.926, 100.451),
         ("Sunan (Hongwansi)", 38.837, 99.616), ("Qilian", 38.177, 100.250), ("Delingha", 37.370, 97.370),
         ("PMO Delingha obs.", 37.378, 97.729), ("Dingxin", 40.300, 99.530), ("JSLC Dongfeng", 40.960, 100.270),
         ("Ejina", 41.958, 101.069), ("Dunhuang", 40.142, 94.662), ("Xining", 36.617, 101.778),
         ("Golmud", 36.402, 94.903), ("Lenghu obs.", 38.606, 93.896)]


# ---------------------------------------------------------------- terrain
class Terrain:
    """ETOPO1 (1 arc-min) grid from NOAA ERDDAP, downloaded once into terrain_etopo1.npz, bilinear lookups.
    Heights are orthometric; the ~-35 m geoid offset here moves the observer < 0.1 km and is ignored."""
    BOX = (34.0, 43.0, 94.0, 103.0)            # lat0, lat1, lon0, lon1

    def __init__(self):
        from scipy.interpolate import RegularGridInterpolator
        if not GRID.exists():
            la0, la1, lo0, lo1 = self.BOX
            url = ("https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.csv?altitude"
                   f"%5B({la0}):1:({la1})%5D%5B({lo0}):1:({lo1})%5D")
            with urllib.request.urlopen(url, timeout=300) as r:
                lines = r.read().decode().splitlines()[2:]
            a = np.array([[float(x) for x in ln.split(",")] for ln in lines])
            lats, lons = np.unique(a[:, 0]), np.unique(a[:, 1])
            H = np.full((len(lats), len(lons)), np.nan)
            H[np.searchsorted(lats, a[:, 0]), np.searchsorted(lons, a[:, 1])] = a[:, 2]
            np.savez_compressed(GRID, lats=lats, lons=lons, h=H)
        g = np.load(GRID)
        self.f = RegularGridInterpolator((g["lats"], g["lons"]), g["h"] / 1e3)

    def h_km(self, la, lo):
        return float(self.f((la, lo)))

    def fetch(self, pts):                      # grid is local; kept for interface symmetry
        return None


# ---------------------------------------------------------------- geometry
def direction_itrs(coord, met):
    t = Time(T0 + timedelta(seconds=float(met)), scale="utc")
    v = coord.transform_to(ITRS(obstime=t)).cartesian.xyz.value
    return v / np.linalg.norm(v)


def ground_hit(r, uu, h_of, rho_max=3000.0):
    """Nearest rho > 0 where r - rho*u reaches the terrain; iterates the terrain height (3 passes)."""
    def f(rho, hg):
        return ecef_to_geodetic(r - rho * uu)[2] - hg

    grid = np.linspace(1.0, rho_max, 300)
    hg = 0.0
    rho = None
    for _ in range(4):
        vals = [f(x, hg) for x in grid]
        k = next((i for i in range(len(vals) - 1) if vals[i] > 0 >= vals[i + 1]), None)
        if k is None:
            return None
        rho = brentq(f, grid[k], grid[k + 1], args=(hg,), xtol=1e-4)
        la, lo, _ = ecef_to_geodetic(r - rho * uu)
        hn = h_of(la, lo)
        if abs(hn - hg) < 0.002:
            break
        hg = hn
    O = r - rho * uu
    la, lo, h = ecef_to_geodetic(O)
    up = geodetic_to_ecef(la, lo, 0.0)
    up /= np.linalg.norm(up)
    el = np.degrees(np.arcsin(uu @ up))
    e = np.array([-np.sin(np.radians(lo)), np.cos(np.radians(lo)), 0.0])
    n = np.cross(up, e)
    az = np.degrees(np.arctan2(uu @ e, uu @ n)) % 360
    return dict(lat=la, lon=lo, h_km=h, range_km=rho, el=el, az=az)


def spline_of(T, X):
    return CubicSpline(T, X)


def load_traj():
    rows = list(csv.DictReader(open(TRAJ)))
    T = np.array([float(r["met"]) for r in rows])
    X = np.array([[float(r[k]) for k in ("x_ecef", "y_ecef", "z_ecef")] for r in rows])
    return T, X, rows


def sun_alt(la, lo, h, met):
    t = Time(T0 + timedelta(seconds=float(met)), scale="utc")
    loc = EarthLocation(lat=la * u.deg, lon=lo * u.deg, height=h * u.km)
    return float(get_sun(t).transform_to(AltAz(obstime=t, location=loc)).alt.deg)


def sun_tangent_km(r, met):
    """Lowest altitude (km, spherical 6371 km Earth) of the ray from the rocket toward the Sun: > 0 geometrically
    sunlit, > ~20 km clear of the opaque lower atmosphere (CERES-2 convention); negative = in the Earth's shadow."""
    t = Time(T0 + timedelta(seconds=float(met)), scale="utc")
    s = get_sun(t).transform_to(ITRS(obstime=t)).cartesian.xyz.to_value(u.km)
    s /= np.linalg.norm(s)
    k = -r @ s
    return float(np.linalg.norm(r + max(k, 0.0) * s) - 6371.0)


def locus(coord, sp, mets, terrain):
    out = []
    for t in mets:
        r = sp(t)
        hit = ground_hit(r, direction_itrs(coord, t), terrain.h_km)
        if hit is None:
            continue
        la_r, lo_r, h_r = ecef_to_geodetic(r)
        hit.update(met=float(t), utc=(T0 + timedelta(seconds=float(t))).isoformat(timespec="seconds"),
                   rocket_lat=la_r, rocket_lon=lo_r, rocket_alt_km=h_r, sun_tangent_km=sun_tangent_km(r, t))
        out.append(hit)
    return out


def radec_from(site, r_rocket, met, coord):
    """Independent check in GCRS: direction site -> rocket (astropy ITRS -> GCRS of both points) against the
    observed ICRS coordinate carried to GCRS as a star (annual aberration applied). -> (RA/Dec string, sep arcsec)"""
    from astropy.coordinates import GCRS, CartesianRepresentation
    t = Time(T0 + timedelta(seconds=float(met)), scale="utc")
    loc = EarthLocation(lat=site["lat"] * u.deg, lon=site["lon"] * u.deg, height=site["h_km"] * u.km)
    rk = SkyCoord(ITRS(CartesianRepresentation(*(r_rocket * u.km)), obstime=t)).transform_to(GCRS(obstime=t))
    og, _ = loc.get_gcrs_posvel(t)
    d = (rk.cartesian - og).xyz.to_value(u.km)
    d /= np.linalg.norm(d)
    seen = SkyCoord(x=d[0], y=d[1], z=d[2], representation_type="cartesian", frame=GCRS(obstime=t))
    star = coord.transform_to(GCRS(obstime=t))
    return seen.to_string("hmsdms", precision=1), float(seen.separation(star).arcsec)


def nearest_towns(points, k=6):
    res = []
    for name, la, lo in TOWNS:
        d = [(C.gc_dist((la, lo), (p["lat"], p["lon"])), p["met"]) for p in points]
        dmin, tmin = min(d)
        res.append(dict(town=name, lat=la, lon=lo, dist_km=round(dmin, 1), met_closest=tmin))
    res.sort(key=lambda x: x["dist_km"])
    return res[:k]


# ---------------------------------------------------------------- sensitivity: re-fly the ascent
def variant_position(kw_fly, met, veh=None, payload=None, inc=None, t0=T0):
    fr = C.Frame(t0)
    z = C.parse_notams(OUT / "notam_jslc_sso")[0]
    ax = C.zone_axis_from_ring(z["ring"])
    inc = C.sso_inclination(600.0) if inc is None else inc
    R = C.fly(fr, veh or C.Y108, payload or 1811.0, ax, inc, 600.0, **kw_fly)
    s1, s2 = R["s1"], R["s2"]
    a = s1["arc"] if met <= s1["t_meco"] else s2["arc"]
    y = a(met)
    return fr.to_ecef(y[:3], met), R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ra", default="21h33m32.7s")
    ap.add_argument("--dec", default="-31d11m16s")
    ap.add_argument("--met", type=float, default=360.0)
    ap.add_argument("--dt", type=float, default=60.0, help="half-width of the MET scan, s")
    ap.add_argument("--of-date", action="store_true", help="treat RA/Dec as equinox-of-date (FK5 J2026.86)")
    ap.add_argument("--sens", action="store_true")
    a = ap.parse_args()

    if a.of_date:
        coord = SkyCoord(a.ra, a.dec, frame=FK5(equinox=Time(T0 + timedelta(seconds=a.met))))
    else:
        coord = SkyCoord(a.ra, a.dec, frame="icrs")
    T, X, _ = load_traj()
    sp = spline_of(T, X)
    terr = Terrain()
    mets = np.arange(a.met - a.dt, a.met + a.dt + 1e-9, 1.0)
    L = locus(coord, sp, mets, terr)
    with open(OUT / "observer_locus.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(L[0]))
        w.writeheader()
        w.writerows([{k: (round(v, 5) if isinstance(v, float) else v) for k, v in p.items()} for p in L])
    best = min(L, key=lambda p: abs(p["met"] - a.met))
    sa = sun_alt(best["lat"], best["lon"], best["h_km"], best["met"])
    print(f"observed RA {a.ra} Dec {a.dec} ({'of date' if a.of_date else 'ICRS/J2000'})")
    print(f"T0+{a.met:.0f} s = {best['utc']} UTC: rocket {best['rocket_lat']:.3f} N {best['rocket_lon']:.3f} E "
          f"h {best['rocket_alt_km']:.1f} km")
    print(f"OBSERVER: {best['lat']:.4f} N {best['lon']:.4f} E, ground {best['h_km'] * 1e3:.0f} m; "
          f"looking az {best['az']:.1f} el {best['el']:.1f} at {best['range_km']:.0f} km; Sun alt {sa:.1f} deg")
    print("time-uncertainty locus:")
    for p in L:
        if abs(p["met"] - a.met) % 10 < 1e-6:
            print(f"   T+{p['met']:5.0f}  {p['lat']:8.4f} N {p['lon']:9.4f} E  ground {p['h_km'] * 1e3:5.0f} m  "
                  f"el {p['el']:5.2f}  az {p['az']:6.2f}  range {p['range_km']:5.0f} km  "
                  f"({C.gc_dist((best['lat'], best['lon']), (p['lat'], p['lon'])):5.0f} km from T+{a.met:.0f})"
                  f"  sun ray {p['sun_tangent_km']:+6.1f} km")
    pick = {}
    for lab, thr in (("geometric", 0.0), ("clear (20 km)", 20.0)):
        p = next((p for p in L if p["sun_tangent_km"] >= thr), None)
        pick[lab] = None if p is None else p["met"]
        if p:
            print(f"rocket sunlit ({lab}) from T+{p['met']:.0f} s -> observer at {p['lat']:.3f} N {p['lon']:.3f} E")
    span = {}
    for dt in (10, 20, 30, 60):
        ps = [p for p in L if abs(p["met"] - a.met) <= dt + 1e-9]
        span[f"+-{dt}s"] = dict(north=[max(ps, key=lambda p: p["lat"])[k] for k in ("lat", "lon")],
                                south=[min(ps, key=lambda p: p["lat"])[k] for k in ("lat", "lon")],
                                length_km=round(C.gc_dist((ps[0]["lat"], ps[0]["lon"]),
                                                          (ps[-1]["lat"], ps[-1]["lon"])), 1))
    rt, sep = radec_from(best, sp(best["met"]), best["met"], coord)
    print(f"round trip (GCRS): rocket seen from the solved site at {rt}, {sep:.1f} arcsec from the observed position")
    near = nearest_towns(L)
    print("closest towns to the locus:")
    for t in near:
        print(f"   {t['town']:20s} {t['dist_km']:6.1f} km  (at T+{t['met_closest']:.0f})")
    res = dict(observed=dict(ra=a.ra, dec=a.dec, frame="FK5 of date" if a.of_date else "ICRS"), met_nominal=a.met,
               best=best, sun_alt_deg=sa, span=span, nearest_towns=near, sunlit_from_met=pick)

    if a.sens:
        print("model sensitivity (observer at T+%.0f s; shift from nominal):" % a.met)
        uu = direction_itrs(coord, a.met)
        variants = [("nominal (1811 kg, SSO, free vernier)", {}, {}),
                    ("payload 1500 kg", {}, dict(payload=1500.0)),
                    ("payload 2100 kg", {}, dict(payload=2100.0)),
                    ("vernier tail 300 s", dict(tail=300.0), {}),
                    ("vernier tail 450 s", dict(tail=450.0), {}),
                    ("plane RAAN -0.2 deg (validation bias)", dict(raan_offset=-0.2), {}),
                    ("i 94.5 deg (no extra dogleg)", {}, dict(inc=94.5)),
                    ("i 96.0 deg", {}, dict(inc=96.0)),
                    ("standard (heavier) vehicle", {}, dict(veh=C.STD, payload=1598.0)),
                    ("T0 11:41:30 (sky 30 s later, same MET)", {}, {})]
        nom = None
        res["sensitivity"] = []
        for name, kwf, kwv in variants:
            try:
                r, R = variant_position(kwf, a.met, **kwv)
            except Exception as e:                     # noqa: BLE001
                print(f"   {name:40s} failed: {e}")
                continue
            hit = ground_hit(r, direction_itrs(coord, a.met + 30.0) if name.startswith("T0 11:41:30") else uu,
                             terr.h_km)
            la_r, lo_r, h_r = ecef_to_geodetic(r)
            if nom is None:
                nom = hit
            d = C.gc_dist((nom["lat"], nom["lon"]), (hit["lat"], hit["lon"]))
            print(f"   {name:40s} rocket {la_r:7.3f} {lo_r:8.3f} h {h_r:6.1f} -> observer {hit['lat']:8.4f} N "
                  f"{hit['lon']:9.4f} E  shift {d:5.1f} km  (MECO-2 T+{R['s2']['t_meco2']:.0f}, margin "
                  f"{R['s2']['margin']:+.0f} kg)", flush=True)
            res["sensitivity"].append(dict(variant=name, rocket=[la_r, lo_r, h_r], observer=[hit["lat"], hit["lon"]],
                                           shift_km=d, t_meco2=R["s2"]["t_meco2"], margin=R["s2"]["margin"]))
    json.dump(res, open(OUT / "observer_solution.json", "w"), indent=1, default=float)
    print("wrote observer_locus.csv, observer_solution.json")


if __name__ == "__main__":
    main()
