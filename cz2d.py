"""CZ-2D (Long March 2D) 3-DOF ascent model for Jiuquan and Taiyuan SSO launches (use_site()).

Stage 1: 4 x YF-20C (YF-21C), zero-lift gravity turn on an Earth-fixed great circle, burned to depletion and
fitted so the empty stage falls on the NOTAM zone centre (zq2e.stage1 / s1_impact, reused unchanged).
Stage 2: YF-24C = YF-22C main engine + 4-chamber YF-23C vernier, hot-staged at MECO.  Phase A: main +
vernier; phase B: vernier only, a long low-thrust climb to direct circular insertion (the CZ-2 family's
"vernier sustainer" phase; there is no coast and no restart).  One bilinear-tangent steering law spans both
phases (optimal control gives tan(steering) bilinear in time independent of thrust level); the main-engine
cutoff time is chosen to minimise propellant, and the final plane node is the fuel-optimal dogleg.

Generic physics (frames, J2-J4 gravity, atmosphere, coast / fall, elements, fit_tle) come from
../launch_traj_NROL-97/nrol97.py; the RK4 kernel, stage-1 law and helpers from ../launch_traj_ZQ-2E_sso/zq2e.py.
Neither is copied.  Units: km, km/s, s, kg, N.  Inertial frame = TEME of date (GMST rotation of ECEF).
"""
import sys
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "launch_traj_ZQ-2E_sso"))
sys.path.insert(0, str(HERE.parent / "launch_traj_CERES-2"))
import zq2e as Z  # noqa: E402
from zq2e import (MU, RE, Arc, Frame, ZoneAxis, bearing, coast, elements, fall, gc_dist, plane_normal,  # noqa: E402,F401
                  plane_through, sso_inclination, stage1, zone_centroid)
from ceres2 import parse_notams  # noqa: E402,F401

G0_MS2 = 9.80665


@dataclass(frozen=True)
class Site:
    """Launch site: pad (GCAT lp table), pad elevation, terrain height under the stage-1 and fairing zones (km)."""
    name: str
    pad: tuple
    pad_h: float
    s1_h: float
    fair_h: float


# JSLC heights are the ones every validated JSLC number was computed with (zq2e defaults); do not change them.
JSLC = Site("JSLC LC-43/94", (40.9605, 100.2980), 1.05, 4.2, 4.2)
# TSLC LC-9: pad ~1.45 km; stage-1 zones in the Qinling south slope (Shaanxi / Hubei border, ~33.4 N 110.6 E,
# ~0.8 km), fairing zones in the west Hubei mountains (~30.2 N 110.1 E, ~1.2 km).  Heights are ESTIMATES; 0.5 km
# of terrain error moves the stage-1 impact by < 0.5 km.
TSLC = Site("TSLC LC-9", (38.8632, 111.5895), 1.45, 0.8, 1.2)
SITE = JSLC
PAD = JSLC.pad


def use_site(site):
    """Select the launch site for everything that follows (pad, pad height, impact terrain)."""
    global SITE, PAD
    SITE, PAD = site, site.pad
    Z.PAD_H = site.pad_h
    return site


@dataclass(frozen=True)
class CZ2D:
    """Published: thrusts, Isp (SAST / CGWIC), stage propellant loads 182 t / 52.7 t, lift-off ~250 t.
    ESTIMATES: stage dry masses, fairing mass, stage-1 vacuum Isp, residual reserve."""
    name: str = "CZ-2D"
    s1_T_sl: float = 2961.6e3      # N, YF-21C at sea level (SAST)
    s1_isp_sl: float = 2550.0 / G0_MS2
    s1_isp_vac: float = 2830.0 / G0_MS2   # ESTIMATE (YF-20 family ~289 s vac)
    s1_prop: float = 182.0e3
    s1_dry: float = 9.8e3         # ESTIMATE incl. interstage struts and residuals (CZ-4 stage 1: 9.5-9.8 t)
    T_main: float = 742.04e3      # N, YF-22C vac
    isp_main: float = 2942.0 / G0_MS2
    T_vern: float = 47.1e3        # N, YF-23C (4 chambers) vac
    isp_vern: float = 2834.0 / G0_MS2
    s2_prop: float = 52.7e3
    s2_dry: float = 4.5e3         # CALIBRATED (validate_precedents.py): 1.9 t to 500 km SSO.  GCAT carries 4006 kg
    #                               for the empty orbital stage; the rest is support bay, adapter, separation hardware
    s2_reserve: float = 300.0     # kg unusable / guidance reserve, ESTIMATE
    fairing: float = 1.1e3        # ESTIMATE, 3.35 m x 6.98 m aluminium fairing
    adapter: float = 0.0          # kg payload support bay counted with the payload side (0 = inside s2_dry)
    diameter: float = 3.35
    t_rise: float = 15.0          # s vertical rise before the kick; at lift-off T/W 1.20 a 10 s rise leaves only
    #                               24 m/s and the zero-lift turn is then set by Coriolis / geodetic tilt, not the kick
    a_max1: float = 1e3           # no throttling on YF-21C
    f_min1: float = 1.0
    dt_sep: float = 1.0           # hot staging: stage 1 separates ~1 s after MECO
    beta_s1: float = 350.0        # kg/m^2 tumbling empty stage 1
    beta_fair: float = 30.0
    q_fair: float = 1135.0        # W/m^2 fairing jettison criterion

    # zq2e.stage1 interface
    @property
    def s1_T_vac(self):
        return self.s1_T_sl * self.s1_isp_vac / self.s1_isp_sl

    @property
    def s1_mdot(self):
        return self.s1_T_sl / (self.s1_isp_sl * G0_MS2)

    @property
    def md_main(self):
        return self.T_main / (self.isp_main * G0_MS2)

    @property
    def md_vern(self):
        return self.T_vern / (self.isp_vern * G0_MS2)

    @property
    def area(self):
        return np.pi * (self.diameter / 2) ** 2

    def m_upper(self, payload):
        return self.s2_prop + self.s2_dry + self.fairing + self.adapter + payload

    def m0(self, payload):
        return self.s1_prop + self.s1_dry + self.m_upper(payload)


STD = CZ2D()
# Y108 (2026-11-11): "首次采用3350复合材料整流罩+减重状态二子级+无尾翼的构型，配套复材支承舱".  The savings are
# not published; ESTIMATES: composite fairing -250 kg, lightened stage 2 + composite support bay -150 kg, no
# stage-1 fins -300 kg (four fins) with a small drag cut that is ignored.
Y108 = replace(STD, name="CZ-2D Y108 (lightened)", fairing=STD.fairing - 250.0, s2_dry=STD.s2_dry - 150.0,
               s1_dry=STD.s1_dry - 300.0)


# ---------------------------------------------------------------- stage 2: main + vernier, then vernier
def _s2P(veh, nhat, x, t_ses, main_on):
    P = np.zeros(Z._NP)
    P[Z.P_NH:Z.P_NH + 3] = nhat
    P[Z.P_PI0], P[Z.P_PI1], P[Z.P_YA0], P[Z.P_YA1] = x[0], x[1], x[2], x[3]
    T = veh.T_main + veh.T_vern if main_on else veh.T_vern
    md = veh.md_main + veh.md_vern if main_on else veh.md_vern
    P[Z.P_TSES], P[Z.P_T2], P[Z.P_MD2], P[Z.P_AMAX2], P[Z.P_FMIN2] = t_ses, T, md, 1e6, 1.0
    return P


def stage2(frame, s1, nhat, x, t_main, t_fair, dt=0.5, tail=None, n_pre=None):
    """MECO -> insertion.  x = [tan-pitch0, tan-pitch1, tan-yaw0, tan-yaw1, t_total/100] with tau = (t-t_ses)/100;
    the main engine runs t_main s (or t_total - tail when tail is given: vernier-only phase of fixed length),
    the vernier the whole t_total; fairing dropped t_fair s after ignition.  n_pre: plane normal held with zero
    yaw until fairing jettison (yaw steering toward nhat starts only then); None = steer to nhat from ignition."""
    veh = s1["veh"]
    a1 = s1["arc"]
    t_ses = s1["t_meco"]
    t_tot = x[4] * 100.0
    tm = min(t_main if tail is None else t_tot - tail, t_tot)
    y = a1.Y[-1].copy()
    y[6] = veh.m_upper(s1["payload"])
    PA, PB = _s2P(veh, nhat, x, t_ses, True), _s2P(veh, nhat, x, t_ses, False)
    P0 = None if n_pre is None else _s2P(veh, n_pre, [x[0], x[1], 0.0, 0.0], t_ses, True)
    cuts = sorted({min(t_fair, t_tot), tm, t_tot})
    arcs, t, fair_state = [], t_ses, None
    for tc in cuts:
        if tc > t - t_ses + 1e-9:
            P = PA if t - t_ses < tm - 1e-9 else PB
            if P0 is not None and t - t_ses < min(t_fair, tm) - 1e-9:
                P = P0
            a = Arc(*Z._march(t, y, t_ses + tc, dt, Z.M_S2, P, -1.0))
            arcs.append(a)
            t, y = a.t_end, a.Y[-1].copy()
        if abs(tc - min(t_fair, t_tot)) < 1e-9 and fair_state is None:
            fair_state = y[:6].copy()
            y[6] -= veh.fairing
    arc = Arc.join(arcs)
    m_ins = arc.Y[-1, 6]
    used = veh.s2_prop - (m_ins - veh.s2_dry - veh.adapter - s1["payload"])
    return dict(arc=arc, t_ses=t_ses, t_meco2=t_ses + tm, t_fair=t_ses + min(t_fair, t_tot), t_ins=t_ses + t_tot,
                m_ins=m_ins, used=used, margin=veh.s2_prop - veh.s2_reserve - used, fair_state=fair_state,
                t_main=tm, x=np.array(x, float))


def _targets(arc, r_t, nhat):
    r, v = arc.r_end, arc.v_end
    rn = np.linalg.norm(r)
    vr = r @ v / rn
    vh = np.linalg.norm(np.cross(r, v)) / rn
    return [(rn - r_t) / 5.0, vr * 200.0, (vh - np.sqrt(MU / r_t)) * 200.0, (r @ nhat) / 5.0, (v @ nhat) * 1e3]


def solve_s2(frame, s1, nhat, r_t, t_main, t_fair, x0=None, tail=None, n_pre=None):
    """Steering + total burn time that end horizontal and circular at radius r_t in plane nhat."""
    starts = [x0] if x0 is not None else []
    starts += [np.array([0.35, -0.12, 0.0, 0.0, 4.2]), np.array([0.6, -0.2, 0.0, 0.0, 5.0]),
               np.array([0.2, -0.06, 0.0, 0.0, 3.6]), np.array([0.9, -0.3, 0.0, 0.0, 5.6])]
    best = None
    for s in starts:
        so = least_squares(lambda p: _targets(stage2(frame, s1, nhat, p, t_main, t_fair, tail=tail,
                                                     n_pre=n_pre)["arc"], r_t, nhat), s, diff_step=1e-6,
                           max_nfev=300)
        if best is None or so.cost < best.cost:
            best = so
        if so.cost < 1e-10:
            break
    if best.cost > 1e-7:
        raise RuntimeError(f"stage-2 targeting did not converge (cost {best.cost:.1e})")
    return stage2(frame, s1, nhat, best.x, t_main, t_fair, tail=tail, n_pre=n_pre), best.x


def heat_flux(arc, t):
    return Z.heat_flux(arc, t)


def fairing_time(s2, q_limit, step=0.5):
    a = s2["arc"]
    for t in np.arange(s2["t_ses"], s2["t_ins"], step):
        if heat_flux(a, t) < q_limit:
            return t - s2["t_ses"]
    return s2["t_ins"] - s2["t_ses"]


def fly(frame, veh, payload, zone_axis, inc, h, pad=None, s1=None, optimal_node=True, t_main=None, t_fair=None,
        cross_km=0.0, tail=None, raan_offset=0.0, hold_plane=False, verbose=False):
    """Stage 1 fitted to the zone, stage 2 to an h-km circular orbit of inclination inc.  Free: main-engine burn
    (propellant-optimal unless t_main or the vernier-only tail length is given), fairing at the heat-flux point unless t_fair given, plane node
    (fuel-optimal dogleg unless optimal_node False -> plane through the stage-2 ignition point).  hold_plane: no
    yaw before fairing jettison (stage 2 stays in the MECO-1 plane until then; see stage2 n_pre)."""
    pad = pad or SITE.pad
    s1 = s1 or fit_stage1(frame, veh, payload, zone_axis, pad=pad, cross_km=cross_km)
    r_t = RE + h
    a1 = s1["arc"]
    n0, raan0 = plane_through(a1.r_end, inc, a1.v_end)
    n_pre = None
    if hold_plane:
        n_pre = np.cross(a1.r_end, a1.v_end)
        n_pre /= np.linalg.norm(n_pre)
    st = dict(x=None, tf=150.0 if t_fair is None else t_fair, tm=185.0 if t_main is None else t_main)

    def run(O, tm):
        s2, st["x"] = solve_s2(frame, s1, plane_normal(inc, O), r_t, tm, st["tf"], st["x"], tail=tail, n_pre=n_pre)
        return s2

    def best_tm(O):
        if t_main is not None or tail is not None:
            s2 = run(O, t_main)
            st["tm"] = s2["t_main"]
            return s2
        m = minimize_scalar(lambda tm: run(O, tm)["used"], bounds=(120.0, 200.0), method="bounded",
                            options=dict(xatol=0.2))
        st["tm"] = float(m.x)
        return run(O, st["tm"])

    raan = raan0
    for _ in range(4):                       # fairing time <-> trajectory fixed point
        if optimal_node:
            m = minimize_scalar(lambda O: best_tm(O)["used"], bracket=(raan - 0.4, raan + 0.4), tol=1e-6)
            raan = float(m.x) % 360
        s2 = best_tm(raan)
        if t_fair is not None:
            break
        tf_new = fairing_time(s2, veh.q_fair)
        if verbose:
            print(f"   fairing {st['tf']:.1f} -> {tf_new:.1f} s after ignition, main {st['tm']:.1f} s, "
                  f"used {s2['used']:.0f} kg")
        done = abs(tf_new - st["tf"]) < 0.6
        st["tf"] = tf_new
        if done:
            break
    raan = (raan + raan_offset) % 360          # sensitivity: plane shifted off the fuel-optimal node
    s2 = best_tm(raan)
    return dict(s1=s1, s2=s2, raan=raan, raan_through_meco=raan0, inc=inc, h=h, frame=frame, veh=veh,
                payload=payload, t_fair_rel=st["tf"])


def capacity(frame, veh, zone_axis, inc, h, pad=None, p_lo=300.0, p_hi=2800.0, tol=5.0, **kw):
    """Payload (kg) for which the propellant margin is zero (stage 1 re-fitted to the zone for each payload)."""
    def margin(p):
        try:
            return fly(frame, veh, p, zone_axis, inc, h, pad=pad, **kw)["s2"]["margin"]
        except RuntimeError:          # stage 1 can no longer reach the zone (max drop range falls with mass)
            return -5e3

    p = brentq(margin, p_lo, p_hi, xtol=tol)
    return p, fly(frame, veh, p, zone_axis, inc, h, pad=pad, **kw)


def s1_impact(frame, s1, h_stop=None):
    """Empty stage 1 to the terrain under the active site's stage-1 zone. -> (t, lat, lon, Arc)."""
    return Z.s1_impact(frame, s1, h_stop=SITE.s1_h if h_stop is None else h_stop)


def fairing_impact(frame, s2, veh, h_stop=None):
    """Fairing half from jettison to the terrain under the active site's fairing zone. -> (t, lat, lon, Arc)."""
    st = s2["fair_state"]
    f = fall(st[:3], st[3:], s2["t_fair"], veh.beta_fair, h_stop=SITE.fair_h if h_stop is None else h_stop,
             dt_hi=1.0, dt_lo=0.5)
    la, lo, _ = frame.lla(f.r_end, f.t_end)
    return f.t_end, la, lo, f


def _fit_to_zone(frame, veh, payload, zone_axis, guess, pad, cross_km):
    """zq2e.fit_stage1_to_zone with the active site's impact terrain height."""
    def res(p):
        _, lat, lon, _ = s1_impact(frame, stage1(frame, veh, p[0], p[1], payload, pad))
        al, cr = zone_axis.along_cross(lat, lon)
        return [al / 10, (cr - cross_km) / 10]

    so = least_squares(res, guess, diff_step=1e-4, x_scale=[1.0, 0.05])
    if so.cost > 1e-6:
        raise RuntimeError(f"stage-1 fit did not converge (cost {so.cost:.2e})")
    return stage1(frame, veh, so.x[0], so.x[1], payload, pad)


def fit_stage1(frame, veh, payload, zone_axis, pad=None, cross_km=0.0, guess=None):
    """Azimuth + kick dropping stage 1 on the zone centre, on the flat (larger-kick) branch: the drop range has a
    maximum against the kick near the 600-640 km CZ-2D zones, and the steep branch lofts MECO above 80 km."""
    pad = pad or SITE.pad
    az0 = bearing(pad, zone_axis.center) - 0.9
    for g in ([guess] if guess else []) + [(az0, 0.27), (az0, 0.32), (az0, 0.23), (az0, 0.40)]:
        try:
            return _fit_to_zone(frame, veh, payload, zone_axis, g, pad, cross_km)
        except RuntimeError:
            continue
    raise RuntimeError("stage-1 zone fit failed")


def zone_axis_from_ring(ring, pad=None):
    return ZoneAxis(np.asarray(ring), pad or SITE.pad)
