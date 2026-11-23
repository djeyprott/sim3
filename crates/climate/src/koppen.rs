//! Köppen-Geiger climate zones from monthly temperature and rainfall,
//! following the rules of Kottek et al. (2006, Meteorol. Z. 15, 259-263).
//! The same rules are applied to the model and to the observations, so the
//! two can be compared, and the viewer can colour the map by zone.

use serde::Serialize;

/// The zones shown: the five main classes with the dry class split into
/// desert and steppe (the difference matters a great deal to people).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum Zone {
    /// A: every month at least 18 °C.
    Tropical,
    /// BW: rainfall below half the dryness threshold.
    Desert,
    /// BS: rainfall below the dryness threshold.
    Steppe,
    /// C: coldest month between -3 and 18 °C.
    Temperate,
    /// D: coldest month -3 °C or colder.
    Continental,
    /// E: warmest month below 10 °C.
    Polar,
}

impl Zone {
    pub const ALL: [Zone; 6] = [
        Zone::Tropical,
        Zone::Desert,
        Zone::Steppe,
        Zone::Temperate,
        Zone::Continental,
        Zone::Polar,
    ];

    pub fn code(self) -> &'static str {
        match self {
            Zone::Tropical => "A",
            Zone::Desert => "BW",
            Zone::Steppe => "BS",
            Zone::Temperate => "C",
            Zone::Continental => "D",
            Zone::Polar => "E",
        }
    }

    pub fn index(self) -> u8 {
        Zone::ALL.iter().position(|&z| z == self).unwrap() as u8
    }

    /// The main class letter (BW and BS are both B).
    pub fn main_class(self) -> char {
        self.code().chars().next().unwrap()
    }
}

/// Classifies one place from its 12 monthly mean temperatures (°C) and
/// rainfall totals (mm). `northern` decides which months are summer
/// (April-September in the north, October-March in the south).
pub fn classify(t: &[f32; 12], p: &[f32; 12], northern: bool) -> Zone {
    let t_max = t.iter().cloned().fold(f32::MIN, f32::max);
    let t_min = t.iter().cloned().fold(f32::MAX, f32::min);
    if t_max < 10.0 {
        return Zone::Polar;
    }
    let t_ann = t.iter().sum::<f32>() / 12.0;
    let p_ann: f32 = p.iter().sum();
    let summer: f32 = (0..12)
        .filter(|&m| (3..9).contains(&m) == northern)
        .map(|m| p[m])
        .sum();
    let winter = p_ann - summer;
    // Dryness threshold (mm): 2 T if at least 2/3 of the rain falls in
    // winter, 2 T + 28 if at least 2/3 falls in summer, else 2 T + 14.
    let p_th = if winter >= 2.0 / 3.0 * p_ann {
        2.0 * t_ann
    } else if summer >= 2.0 / 3.0 * p_ann {
        2.0 * t_ann + 28.0
    } else {
        2.0 * t_ann + 14.0
    };
    if p_ann < 10.0 * p_th {
        return if p_ann < 5.0 * p_th {
            Zone::Desert
        } else {
            Zone::Steppe
        };
    }
    if t_min >= 18.0 {
        Zone::Tropical
    } else if t_min > -3.0 {
        Zone::Temperate
    } else {
        Zone::Continental
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_places() {
        // Monthly normals 1991-2020, rounded (national weather services).
        // London: C.
        let t = [
            5.0, 5.5, 7.5, 10.0, 13.5, 16.5, 19.0, 18.5, 15.5, 12.0, 8.0, 5.5,
        ];
        let p = [
            55.0, 41.0, 41.0, 44.0, 49.0, 45.0, 45.0, 50.0, 49.0, 69.0, 59.0, 56.0,
        ];
        assert_eq!(classify(&t, &p, true), Zone::Temperate);
        // Cairo: BW.
        let t = [
            14.0, 15.0, 18.0, 21.5, 25.0, 27.5, 28.5, 28.5, 26.5, 24.0, 19.5, 15.5,
        ];
        let p = [5.0, 4.0, 4.0, 1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0, 4.0, 6.0];
        assert_eq!(classify(&t, &p, true), Zone::Desert);
        // Singapore: A.
        let t = [
            26.5, 27.0, 27.5, 28.0, 28.5, 28.5, 28.0, 28.0, 27.5, 27.5, 27.0, 26.5,
        ];
        let p = [
            220.0, 105.0, 150.0, 160.0, 170.0, 130.0, 160.0, 175.0, 170.0, 190.0, 250.0, 290.0,
        ];
        assert_eq!(classify(&t, &p, true), Zone::Tropical);
        // Moscow: D.
        let t = [
            -6.5, -6.0, -1.0, 6.5, 13.0, 17.0, 19.5, 17.5, 11.5, 5.5, -0.5, -4.5,
        ];
        let p = [
            53.0, 44.0, 39.0, 37.0, 61.0, 78.0, 84.0, 78.0, 66.0, 70.0, 52.0, 51.0,
        ];
        assert_eq!(classify(&t, &p, true), Zone::Continental);
        // Svalbard airport: E.
        let t = [
            -12.0, -12.5, -12.0, -8.5, -2.5, 3.0, 7.0, 5.5, 1.5, -4.0, -7.5, -10.5,
        ];
        let p = [20.0; 12];
        assert_eq!(classify(&t, &p, true), Zone::Polar);
        // Alice Springs (southern hemisphere, summer rain): BW.
        let t = [
            29.0, 28.0, 25.0, 20.5, 15.5, 12.5, 12.0, 14.5, 19.0, 23.0, 26.0, 28.0,
        ];
        let p = [
            40.0, 45.0, 30.0, 15.0, 15.0, 13.0, 10.0, 8.0, 8.0, 20.0, 25.0, 40.0,
        ];
        assert_eq!(classify(&t, &p, false), Zone::Desert);
    }
}
