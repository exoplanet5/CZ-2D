# CZ-2D launch profiles: ascent models and 3D page

Ascent and insertion-orbit models of two Long March 2D launches, built from their NOTAM drop zones, with an
interactive 3D page. They are models, not telemetry and not official ephemerides.

**3D page: https://exoplanet5.github.io/CZ-2D/**

| Mission | Lift-off (UTC) | Orbit modelled | Notes | 3D page |
|---|---|---|---|---|
| [CZ-2D, Taiyuan LC-9](missions/TSLC_20261118/README.md) | 2026-11-18 22:51 (assumed) | 526 km, 97.49° | Pre-launch; dawn twilight plume | [open](https://exoplanet5.github.io/CZ-2D/?profile=sso&mission=cz2d-tslc-1118) |
| [CZ-2D Y108, Jiuquan LC-43/94](missions/Y108_JSLC_20261111/README.md) | 2026-11-11 11:41 (official) | 600 km SSO (assumed) | Flown; observer located from a T0+6 min sighting | [open](https://exoplanet5.github.io/CZ-2D/?profile=sso&mission=cz2d-y108) |

Each mission folder holds its NOTAM text, results, map and a README with the assumptions and uncertainty.

## Layout

```
cz2d.py                  CZ-2D vehicle and ascent model (Jiuquan and Taiyuan)
validate_precedents.py   capacity, orbit plane and fairing checks on eight flown launches -> validation.json
build_traj.py            trajectory, insertion orbit, TLE and KML for one mission
twilight_visibility.py   sunlit-plume geometry and the launch visibility map
observer_locate.py       Y108: back-projects the sighting onto terrain
make_map.py              Y108: observer map
web_export.py            mission data for the 3D page
missions/<mission>/      NOTAM file, outputs, map, README
precedents/              first TLEs of the validation flights (Space-Track.org)
docs/                    the 3D page (GitHub Pages serves main -> /docs)
```

## Model in brief

- **Vehicle:** published CZ-2D figures (CGWIC / SAST). The stage-2 dry mass is the one calibrated value, set so the
  model carries 1.9 t to 500 km SSO; it then carries 1252 kg to 700 km against 1.3 t published.
- **Stage 1:** a zero-lift turn fitted so the empty stage falls on the centre of its NOTAM zone.
- **Stage 2:** main engine plus vernier, then the vernier alone in a long climb to direct circular insertion. The
  main-engine cutoff is the propellant-optimal one, and the plane change is the fuel-optimal dogleg.
- **Validation:** on eight flown launches the predicted RAAN is within −0.10° to +0.32° of the first TLE, with
  lift-off known to the minute (0.25° per minute). Where a fairing zone exists the modelled fairing lands
  15–23 km from its centre on three flights and 123 km on one (Y105, a different profile).

## Run

```sh
python3 tools/serve.py                      # the page, at http://127.0.0.1:8770/
python validate_precedents.py --quick       # about 2 min
python build_traj.py --mission tslc1118     # or y108
python twilight_visibility.py --mission tslc1118
python web_export.py                        # both missions -> docs/data/
```

The Python needs numpy, scipy, numba, astropy, skyfield, simplekml and matplotlib. The model also imports shared
physics from sibling projects that are not in this repository (`launch_traj_ZQ-2E_sso/zq2e.py`,
`launch_traj_NROL-97/nrol97.py`, `launch_traj_CERES-2/ceres2.py`), the DE440s kernel and a plotting style module,
so the page and the results here can be read without them but the scripts cannot be rerun from this repository
alone.

## Sources

- NOTAMs: the text files in each mission folder.
- CZ-2D data: CGWIC and SAST LM-2D pages; Xinhua 2025-09-29 (1.9 t to 500 km SSO); CASC (1.3 t to 700 km SSO).
- First TLEs of the validation flights: Space-Track.org.
- Launch times, pads and satellite masses: data from GCAT (J. McDowell, planet4589.org/space/gcat), CC-BY 4.0.
- Sun and Moon: JPL DE440s. Terrain: NOAA ETOPO1. Base maps: NASA Blue Marble and Black Marble.
- Page licences: `docs/LICENSE-CZ-6A-viewer.txt`, `docs/tex/TILES-LICENSE.txt`, `docs/data/*-LICENSE.txt`.
