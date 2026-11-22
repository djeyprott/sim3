# sim3: project rules

A human history simulator starting around 10,000 BCE on Earth's real terrain. The full design is in `docs/design.md`; read it before planning any work.

## Working with the owner

- The owner directs and approves; Claude writes all code. The owner does not read code, so explain changes in terms of what the simulation does and what they will see in the viewer.
- **Every step is planned and discussed with the owner before anything is pushed.** Propose the plan, wait for approval, build, show the result, wait for approval again, then push.
- Approved work goes straight to `main`. No pull requests.
- The owner runs Windows. Anything they need to do must work there without installing programming tools.

## Simulation rules

- **No arbitrary values.** Parameters have real-world units and a stated source (a comment citing the paper, dataset or reasoning). No made-up scores like "stability +5".
- **Derived, not declared.** Aggregate qualities (stability, prosperity, legitimacy) are computed from underlying facts, never stored as free-floating numbers.
- Every phase ends with something the owner can see in the viewer.
- Every event that matters goes into the `history` event log, with its causes.

## Code layout

- `crates/history`: event log.
- `crates/sim`: simulation core, time (`time.rs`; one tick = one day; proleptic Gregorian calendar; astronomical year numbering internally, BCE/CE only for display). Systems run at their own cadence (daily, seasonal, yearly), not all every day.
- `crates/server`: the `sim3` binary. Runs the simulation on a background thread and serves the viewer on 127.0.0.1:7340.
- `crates/world`: geodesic grid (`grid.rs`; map grid level 8 = 655,362 cells, climate grid level 6 = 40,962 cells; indices stable across levels), world data file parsing (`data.rs`), sea/land/inland-basin classification (`terrain.rs`). `gridgen` binary exports cell centres for the data pipeline.
- `viewer/`: plain HTML/CSS/JS, embedded into the binary at compile time. `map.js` draws the map with WebGL2 from a pixel-to-cell lookup image. Views can be linked with `#lat,lon,zoom`.
- New systems (world, climate, ...) get their own crate under `crates/`.

## World data

Raw datasets (ETOPO 2022, ICE-6G_C) can't be downloaded from Claude's workspace, only from GitHub's machines. `.github/workflows/data.yml` runs `tools/data/build_world.py` there and publishes `world-data-v<N>` (N in `tools/data/VERSION`); bump VERSION to publish a new file, never overwrite. The server embeds `data/world.bin` (gitignored); download it from that release before building. CI and release workflows do this; publish new data before pushing code that needs it. Job logs can be read without auth at `https://api.github.com/repos/djeyprott/sim3/actions/jobs/<id>/logs`.

## Releases

Pushing to `main` with a new `version` under `[workspace.package]` in `Cargo.toml` makes GitHub build `sim3.exe` and publish it as release `v<version>`. Bump the version whenever a change should reach the owner as a download. Tags cannot be pushed from Claude's workspace, so don't rely on them.

## Checks before proposing a push

```
cargo fmt --all
cargo clippy --all-targets -- -D warnings
cargo test
```

Run the binary and look at the viewer (e.g. a headless browser screenshot) before telling the owner it works.
