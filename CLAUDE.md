# CLAUDE.md

CZ-2D Y108 (探源计划探路者卫星, CAS "Origins" Pathfinder) from JSLC LC-43/94, **2026-11-11 19:41 BJT = 11:41 UTC**
(official; seconds unknown, model uses 11:41:00). Single NOTAM A5728/26 (`notam_jslc_sso`, 11:32-12:05 UTC) = stage-1
drop zone. The user asked for a **600 km SSO** model first, then the ground-observer puzzle (RA/Dec at T0+6 min).
Results, assumptions and sources are in `missions/Y108_JSLC_20261111/README.md`.

Second mission: **CZ-2D from TSLC LC-9, 2026-11-18** (`notam_tslc_sso`: A5897/26 stage 1, A5898/26 fairing), modelled
as the PIESAT-2 follow-up of 2026-218 at 526 km / 97.49 deg, with the dawn twilight plume. See the TSLC section below
and `missions/TSLC_20261118/README.md`.

**Layout:** scripts and the shared model sit in the root; each mission's NOTAM file and every output live in
`missions/<dir>/` (`Y108_JSLC_20261111`, `TSLC_20261118`; `build_traj.mission_dir(key)`). The root `README.md` is the
repository overview. Git repo `exoplanet5/CZ-2D` (public); the 3D page in `docs/` is served by GitHub Pages at
https://exoplanet5.github.io/CZ-2D/ from `main` -> `/docs`.

Everything runs with `~/.venvs/astro313/bin/python <script>.py`.

## Pipeline (order matters)

1. `validate_precedents.py [--quick] [--site JSLC|TSLC]` (~5 min; `--site TSLC` 20 s): capacity calibration +
   plane/closure on JSLC Y97, Y101, Y100, Y42, Y41 and TSLC Y99, Y105, Y104 (these also check the fairing against
   its zone). Run it before trusting any change to `cz2d.py`. Writes `validation.json`; `--site` keeps the other
   site's rows.
2. `build_traj.py [--mission y108|tslc1118] [--t0 ISO] [--alt KM] [--inc I] [--payload KG] [--stem S]` (y108 ~2 min,
   tslc1118 40 s): writes `<stem>_{solution.json,traj.csv,orbit.csv,traj.kml,SSO<alt>.tle}` into the mission folder.
   Missions are the `MISSIONS` table (folder, site, NOTAM file, T0, target, vehicle, payload, RAAN calibration).
   y108 payload = lightened capacity. The solution JSON carries the stage-1 and fairing fall tracks for the page.
3. `observer_locate.py [--met 360 --dt 60] [--of-date] [--sens]`: back-projects the sighting onto ETOPO1 terrain.
   `--sens` re-flies 9 variants (~10 min). Writes `observer_locus.csv`, `observer_solution.json`.
4. `make_map.py`: `CZ-2D_Y108_observer_map.png` (launch-plot export rules from the global CLAUDE.md).
5. `twilight_visibility.py [--mission tslc1118] [--min-elev 10] [--extent W E S N]`: Sun-ray tangent height
   along the trajectory, sunlight pickup, per-city Sun / elevation table (`<stem>_twilight.json`) and the launch
   visibility map (`<stem>_visibility.png`). Steps 3-4 are Y108 only.
6. `web_export.py [--mission tslc1118|y108|all]`: mission data for the 3D page, `docs/data/cz2d-tslc-1118.json` and
   `docs/data/cz2d-y108.json`, plus `<stem>_page_ref.json` (pickup and 10 deg zone areas by the page's definitions).
   Page text, cities, views and notes per mission are the `PAGE` table in that file.

**3D page (`docs/`, launch-page skill):** created with `new_page.py` from `~/launch/launch-page-template`; one
profile `sso`, missions `cz2d-tslc-1118` (default) and `cz2d-y108`. `solid_end_met` = MECO-2 with
`plume_names: ['main engine', 'vernier']` (a template field added for this page): the bold lit arc and white zone are
the main-engine burn, the thin arc and green zone the vernier climb. On Y108 the main burn is in shadow, so its
white zone and viewers table are empty by design. Page vs Python on 2026-11-12: pickup equal to 0.001 s, zone areas
equal to 4 digits (TSLC 2.037 M km2 main, Y108 13.17 M km2 vernier); full `check_page.py` 22 of 22. After changing a mission: `build_traj.py`,
`web_export.py`, `check_page.py . --quick`, then `publish_page.py . --repo exoplanet5/CZ-2D --message "..."`.
In scripts, import `build_traj` before `cz2d`: cz2d puts sibling projects first on `sys.path`, and CERES-2 has its
own `build_traj.py`.

**Sites:** `cz2d.use_site(C.JSLC | C.TSLC)` sets the pad, pad height (it writes `zq2e.PAD_H`) and the terrain height
under the stage-1 and fairing zones. Call it before `Frame`-dependent work; `C.PAD` follows it. The JSLC heights
(1.05 / 4.2 / 4.2 km) reproduce the Y108 outputs bit-for-bit; do not change them.

`cz2d.py` is the core: CZ-2D vehicle (`STD`, `Y108`), two-phase stage 2 (main+vernier, then vernier only), capacity.
Reused, not copied: `../launch_traj_ZQ-2E_sso/zq2e.py` (RK4 kernel, stage-1 zero-lift turn, zone fit, plane helpers),
`../launch_traj_NROL-97/nrol97.py` (frames, gravity, atmosphere, coast/fall, elements, `fit_tle`),
`../launch_traj_CERES-2/ceres2.py` (`parse_notams`). Terrain grid `terrain_etopo1.npz` comes from NOAA ERDDAP
`etopo180` (one download, 34-43 N, 94-103 E); OpenTopoData rate-limits (HTTP 429), do not go back to per-point calls.

## Model choices (validated; do not regress)

- **Vehicle (CGWIC/SAST):** S1 182 t, YF-21C 2961.6 kN SL, 2550 m/s; S2 52.7 t, YF-22C 742.04 kN / 2942 m/s +
  YF-23C vernier 47.1 kN / 2834 m/s; lift-off ~250 t. Hot staging at MECO-1 (156.7 s burn).
- **S2 dry 4.5 t is the one calibrated knob** (1.9 t to 500 km SSO, Xinhua 2025). Check point: 700 km -> 1252 kg vs
  1.3 t published. GCAT's 4006 kg is the empty orbital stage only.
- **t_rise 15 s.** Lift-off T/W is only 1.20; after a 10 s rise (24 m/s) the zero-lift turn is driven by Coriolis /
  geodetic tilt, not the kick, and MECO-1 ends below 40 km. Stage 1 fits the flat (larger-kick) branch; the 600-640 km
  CZ-2D zones sit near the maximum drop range, which falls as payload grows (capacity() treats "zone unreachable" as
  a negative margin).
- **Main-engine cutoff is propellant-optimal -> long vernier climb** (MECO-2 ~T+317 s, insertion ~T+874 s). This is
  what reproduces the published 500 -> 700 km capacity drop (model 650 kg, published 600 kg). Heritage-length vernier
  tails (170-320 s) cannot reach 700 km and cut 500 km capacity by 0.3-1.4 t. One bilinear-tangent law spans both phases.
- **Fairing at the 1135 W/m2 heat-flux point** (~51 s after S2 ignition). Independent check: Y41 fairing falls 15 km
  from its 894 km NOTAM zone.
- **Plane:** fuel-optimal node. RAAN vs first TLE: +0.09..+0.32 deg on the five precedents (T0 known to the minute),
  i.e. the model plane sits ~0.2 deg east. Standard SSO corridors fit at az 194.6-195.6 (MECO-1 plane i 93.8-94.9 ->
  97.4); Y41 fits 184.9 (86.2 -> 89.0). **Y108 fits 191.3 (MECO-1 plane 91.4):** SSO needs ~3.4 deg more dogleg than
  usual (~200 kg). The corridor alone points to i ~ 94-95; SSO is the user's assumption.

## Y108 result (T0 11:41:00, 600 km SSO, payload 1811 kg = lightened capacity; standard vehicle 1598 kg)

MECO-1 T+156.7 63.6 km -> zone centre; fairing T+207.7; MECO-2 T+317.2 228 km; **T+360 33.518 N 98.358 E 288.6 km**;
insertion T+874.4 6.60 N 92.39 E, i 97.788, RAAN 141.03.

## Observer (RA 21h33m32.7s Dec -31d11m16s, ICRS)

- T+360: **39.675 N 98.616 E, 1449 m**, az 182.0 el 19.4, range 755 km, Sun -18 deg. 12 km ESE of Jiuquan city.
  Round-trip check (astropy GCRS) reproduces the RA/Dec to 0.0".
- Time locus runs NNE-SSW, ~3.5-4.5 km per second of MET (±10 s = 35/46 km, ±30 s = 98/124 km).
- **The rocket leaves the Earth's shadow at T+360 s** (Sun-ray tangent +0.4 km; clear sunlight +20 km at T+372). A
  plume 755 km away at night is seen only once sunlit, so times before ~T+360 (the northern half, Jinta side) are
  unlikely. South of ~39.4 N the locus climbs into the Qilian Shan at 4-4.6 km.
- Model sensitivity at T+360: payload ±300 kg 11 km; RAAN -0.2 deg 5.5 km W; vernier tail 450 s 10 km; i 94.5 32 km E;
  standard vehicle 20 km; equinox-of-date reading of RA/Dec 6.5 km NE; T0 +30 s 1.4 km.

## TSLC 2026-11-18 (T0 22:51:00 UTC assumed, 526 km, i 97.49, payload 1200 kg, standard vehicle)

- **Identification:** the stage-1 zone equals Y99 (2024-12-16) and Y104 (2026-09-19 = 2026-218) corner for corner;
  both went to 518 x 533 km, 97.49 deg. Y105 (Superview Neo 2, 2026-03-25) used other zones. GCAT lists 2026-218 at
  JQ, but Space-Track (TSC) and the zones say Taiyuan.
- **T0:** stage-1 NOTAM start + 8 min (Y99, Y105; Y104 + 7). No official time yet; update `MISSIONS` when known.
- **Result:** az 187.66, MECO-1 plane i 88.4 (9 deg dogleg); MECO-1 T+156.7 64 km; fairing T+208.2; MECO-2 T+316.7
  218 km at 33.31 N 110.46 E; insertion T+840.4 at 5.33 N 104.73 E; RAAN 328.28 (LTAN 18:17); margin 672 kg.
- **RAAN calibration +0.09 deg** (`raan_offset`): Y99 -0.09, Y104 -0.10 against first TLEs. JSLC flights sit at
  +0.09..+0.32, so the bias is corridor-specific; do not apply the TSLC value to JSLC.
- **Fairing:** model lands 4 km long, 15 km west of the A5898/26 axis (zone 72 x 35 km). `fly(hold_plane=True)`
  (no yaw before fairing jettison) lands 32 km east and shifts RAAN +0.45 deg, and it is worse on JSLC Y41 too, so
  the fuel-optimal yaw from ignition stays the default. Y105's fairing zone is missed by 123 km along track
  (different profile, unexplained).
- **Twilight:** Sun -6.1 deg at the pad at T0. Sunlit from T+115 s (geometric) / T+141 s at 50 km (ray >= 20 km),
  i.e. before staging; the whole main burn is sunlit over twilight ground. Map reference epoch = middle of the
  sunlit main-engine arc (22:54:56 UTC); the visibility zone covers pickup -> MECO-2 only (vernier plume is 6 % of
  the thrust).

