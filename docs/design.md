# Human History Simulator: Design

Copy of the working design doc. When the two differ, discuss which is right and update both.

## Core design principle

"Avoid arbitrary values" means two rules:

- **Parameters have real units and can be calibrated.** Calories per hectare of emmer wheat, birth intervals in months, the share of people a lord can feed from a harvest, travel days by terrain. Each can be checked against archaeological and historical data, and each states its source.
- **Derived, not declared.** A state's "stability", a region's "prosperity" or a ruler's "legitimacy" are never stored numbers. They are computed from the underlying facts: who holds power, whether those people are being paid, whether the harvest failed, whether the clergy backs the ruler.

## Architecture

Ten systems, each with a thin first version and a clean interface, deepened over time.

### 1. World and climate

Earth's real terrain to start; generated planets (with simplified plate tectonics) can come later and reuse the same climate model. Climate comes from physics, not painted biomes: incoming sunlight by latitude, an energy-balance model, prevailing wind belts, ocean currents, rain shadows behind mountains, and climate zones derived from all of that. Variability is layered on as real mechanisms:

- solar cycles
- ENSO-like ocean oscillations (like El Niño) and multi-decade oscillations
- volcanic winters
- slow orbital forcing

Droughts, little ice ages and monsoon failures follow from these without being scripted. Starting at 10,000 BCE means the model must also handle the end of the last ice age: the Younger Dryas cold spell ending, the remaining ice sheets melting, and sea level rising, which floods land bridges such as Doggerland.

### 2. Ecology and carrying capacity

Plants, game, soils and water set how many people a place can feed under each way of making a living: foraging, herding, farming, irrigation. Climate shifts move those limits, and that is where upheavals come from.

### 3. Demography: a hybrid model

- **Population groups**: cohorts by age, sex, culture, religion, language, occupation and location, following real life tables, fertility and migration.
- **Historical figures**: full individuals, promoted out of the population groups when they become significant.

### 4. Culture, language and religion

- **Culture**: practices and values that drift, spread on contact and get adopted when they work.
- **Language**: a real sound system plus vocabulary, with sound-change rules. Families split on their own; names come from the generated languages; mutual intelligibility is measured by word overlap.
- **Religion**: beliefs and practices, a clergy as an organization, syncretism, schisms and conversion driven by social incentives.

### 5. Individuals

Personality (HEXACO or Big Five), needs, ambitions, memories, and relationships whose opinions come from those memories. Goals are chosen by weighing utility against personality.

### 6. Technology as knowledge, not a fixed tree

Techniques have prerequisites (materials, earlier techniques, a need). Discovery chances grow with practitioners and contact. Knowledge spreads through trade and migration and can be lost when populations shrink. The order of discoveries emerges.

### 7. Organizations

One generic structure covers states, tribes, guilds, banks, churches and gangs: members, roles, resources, decision and succession rules, optional territory. Internal politics follow selectorate theory. Outsiders can back factions.

### 8. Economy

Production, trade routes, prices from supply and demand, credit and debt, taxes.

### 9. Diplomacy

Person to person (marriages, friendships, feuds, patronage) and organization to organization (treaties, tribute, alliances, wars), driven by leaders' own interests filtered through the organization's needs.

### 10. Legends viewer

Every event goes into a structured log with its causes. Browse people, places, organizations, eras, cultures and languages; view maps for any year; read generated prose. The viewer runs in the browser and loads only what is on screen, so long histories stay responsive and the simulation's speed is unaffected.

## Phases

| Phase | Deliverable |
| --- | --- |
| 0 | Design doc, tech stack, repo setup |
| 1 | Terrain and climate model, with a map viewer and climate playback over centuries |
| 2 | Ecology, carrying capacity, first population groups and migration |
| 3 | Event log and legends viewer skeleton |
| 4 | Culture, language and religion dynamics |
| 5 | Historical figures |
| 6 | Organizations and internal politics |
| 7 | Economy and technology |
| 8 | Diplomacy and conflict |
| 9+ | Deepening, calibration against real data, performance |

## Decisions

1. **Stack.** Native Rust simulation; browser viewer served locally by the same program. No terminal viewer. Claude writes the code; the owner directs and approves.
2. **World.** Earth's real terrain with our own physically based climate model, which must also work for generated worlds later.
3. **Start.** About 10,000 BCE. Runs are open-ended.
4. **Run style.** Live, with pause and inspect, and browsing of everything that has already happened.
5. **Repo.** github.com/djeyprott/sim3, commits straight to `main` after approval. License GPL-3.0-or-later.
6. **Platform.** The owner runs Windows; releases ship a ready-made `sim3.exe`.
7. **Time step.** One tick is one day, labelled with the proleptic Gregorian calendar. Each system runs at the cadence it needs: daily for weather and people's actions, seasonally or yearly for slow processes such as demography, harvests, ice sheets and sea level.
8. **Climate cadence.** The climate updates once per simulated month (decided 22 Nov for speed); daily weather is derived from it where needed (Step 5). Climate gets at most ~10% of simulation time.
9. **Climate history.** Eruptions, solar changes and other shifts are random with realistic statistics; the ice sheets respond to the simulated climate; CO2 follows the ice-core record (Bereiter et al. 2015).
10. **Present-day check** is a run the owner can start in the viewer (1931-2010, compared with CHELSA 1981-2010).
11. **No Earth-specific shortcuts.** Climate (and every later system) is driven only by the terrain data, never by knowing where Earth's continents are, so generated worlds get the same treatment. Earth is the test case.
