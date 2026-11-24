# sim3

A human history simulator. It starts around 10,000 BCE on Earth's real terrain and simulates climate, populations, cultures, languages, religions, technology, organizations and the people who lead them, with a browser viewer for watching history unfold and reviewing what happened.

Status: **Phase 1, step 4**: Earth's terrain, sea and ice from 1 January 10,000 BCE, with a monthly climate driven by sunlight and the Earth's orbit: air and sea temperatures, winds, wind-driven ocean currents, the Atlantic-style overturning circulation, sea ice, rain, snow and Köppen climate zones. A present-day check run compares the model with observed temperatures, rainfall, climate zones and sea ice. Known misses (the North Atlantic's warmth at Bergen, the monsoon, sea ice amounts) are left for a later climate polish step. See [docs/design.md](docs/design.md) for the plan.

## Running it on Windows

1. Open the [Releases page](https://github.com/djeyprott/sim3/releases) and download `sim3.exe` from the newest release.
2. Double-click `sim3.exe`. A console window opens and the viewer opens in your browser.
   - Windows may warn that the app is from an unknown publisher. Click **More info**, then **Run anyway**.
3. To stop the simulation, close the console window.

If the browser doesn't open, go to http://127.0.0.1:7340 yourself. The viewer only works on this computer; nothing is reachable from the network.

## Map links

A view can be opened directly by adding `#latitude,longitude,zoom` to the address, for example http://127.0.0.1:7340/#54,3,9 for the North Sea.

## Data sources

- Terrain today: NOAA NCEI (2022), ETOPO 2022 Global Relief Model, https://doi.org/10.25921/fd45-gt74
- Ice-age topography and ice extent: Peltier, Argus & Drummond (2015), ICE-6G_C (VM5a), J. Geophys. Res. Solid Earth 120, 450-487; Argus et al. (2014), Geophys. J. Int. 198, 537-563. Distributed by PMIP4.

- Observed temperatures and rainfall 1981-2010: CHELSA V2.1 (Karger et al. 2017), https://doi.org/10.16904/envidat.228.v2.1
- Earth's orbit: Berger & Loutre (1991), via NOAA NCEI.
- CO2: Antarctic ice-core composite (Bereiter et al. 2015) and Mauna Loa (NOAA GML).

Ice thickness in 10,000 BCE is sim3's own estimate from these data (see `tools/data/build_world.py`).

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
