//! The event log.
//!
//! Every system records what happens here, so the viewer can show it later.
//! Phase 0 keeps the log in memory; Phase 3 moves it to an on-disk database
//! so long runs don't have to fit in RAM.

use serde::Serialize;

/// Days elapsed since the simulation started. This is the canonical time
/// stamp; calendar dates are derived from it.
pub type Tick = u64;

pub type EventId = u64;

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum EventKind {
    /// The simulation began at this tick.
    SimulationStarted,
    /// A volcano put this much sulfur into the stratosphere (Tg S).
    VolcanicEruption {
        lat: f32,
        lon: f32,
        sulfur_tg: f32,
    },
    /// The world's yearly average fell this far below the recent average.
    WorldCooled {
        cooling_c: f32,
    },
    /// A strong El Niño (warm) or La Niña (cold) peaked: the eastern
    /// equatorial sea anomaly (°C) and the places where the last three
    /// months' rain fell furthest below and above normal (% of normal).
    ElNino(EnsoReport),
    LaNina(EnsoReport),
}

#[derive(Debug, Clone, Serialize)]
pub struct EnsoReport {
    /// Longitude of the middle of the ocean basin.
    pub basin_lon: f32,
    pub strength_c: f32,
    pub dry_lat: f32,
    pub dry_lon: f32,
    pub dry_percent: f32,
    pub wet_lat: f32,
    pub wet_lon: f32,
    pub wet_percent: f32,
}

#[derive(Debug, Clone, Serialize)]
pub struct Event {
    pub id: EventId,
    pub tick: Tick,
    #[serde(flatten)]
    pub kind: EventKind,
    /// Events that led to this one. Empty for root events.
    pub causes: Vec<EventId>,
}

#[derive(Debug, Default)]
pub struct EventLog {
    events: Vec<Event>,
}

impl EventLog {
    pub fn new() -> Self {
        Self::default()
    }

    /// Records an event and returns its id. Ids are assigned in order, so
    /// they double as a cursor for "everything since X".
    pub fn record(&mut self, tick: Tick, kind: EventKind, causes: Vec<EventId>) -> EventId {
        let id = self.events.len() as EventId;
        self.events.push(Event {
            id,
            tick,
            kind,
            causes,
        });
        id
    }

    pub fn len(&self) -> usize {
        self.events.len()
    }

    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }

    /// Events with an id greater than or equal to `from`, at most `limit` of them.
    pub fn from_id(&self, from: EventId, limit: usize) -> &[Event] {
        let start = (from as usize).min(self.events.len());
        let end = start.saturating_add(limit).min(self.events.len());
        &self.events[start..end]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ids_are_sequential_and_paging_works() {
        let mut log = EventLog::new();
        for t in 0..5 {
            assert_eq!(log.record(t, EventKind::SimulationStarted, vec![]), t);
        }
        assert_eq!(log.from_id(3, 10).len(), 2);
        assert_eq!(log.from_id(0, 2).len(), 2);
        assert!(log.from_id(99, 10).is_empty());
    }
}
