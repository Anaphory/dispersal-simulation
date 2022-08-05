use std::ops::{Div, Index, Mul};

/**
OneYearResources is a number with units.
*/
#[derive(Default, Clone, Copy, Serialize, Deserialize)]
pub struct OneYearResources {
    c: f64,
}
/**
We need comparability, so we cheat this implementation.
 */
impl PartialEq for OneYearResources {
    fn eq(&self, other: &Self) -> bool {
        self.c == other.c
    }
}
impl Eq for OneYearResources {}
impl PartialOrd for OneYearResources {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.c.partial_cmp(&other.c) {
            None => Some(std::cmp::Ordering::Equal),
            Some(c) => Some(c),
        }
    }
}
impl Ord for OneYearResources {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.c.partial_cmp(&other.c) {
            None => std::cmp::Ordering::Equal,
            Some(c) => c,
        }
    }
}

/**
Attach the unit to a number
*/
impl From<f64> for OneYearResources {
    fn from(c: f64) -> Self {
        OneYearResources { c: c }
    }
}
/**
Numbers with units can be multiplied/divided by scalars
 */
impl Mul<f64> for OneYearResources {
    type Output = Self;
    fn mul(self, f: f64) -> Self {
        OneYearResources { c: self.c * f }
    }
}
impl Div for OneYearResources {
    type Output = f64;
    fn div(self, f: Self) -> f64 {
        self.c / f.c
    }
}
/**
Numbers with unit can be added, and subtracted
 */
impl std::ops::Sub for OneYearResources {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        OneYearResources {
            c: self.c - other.c,
        }
    }
}
impl std::ops::AddAssign for OneYearResources {
    fn add_assign(&mut self, other: Self) {
        self.c += other.c
    }
}
impl std::iter::Sum for OneYearResources {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        OneYearResources {
            c: iter.map(|o| o.c).sum(),
        }
    }
}
impl<'a> std::iter::Sum<&'a OneYearResources> for OneYearResources {
    fn sum<I: Iterator<Item = &'a Self>>(iter: I) -> Self {
        OneYearResources {
            c: iter.map(|o| o.c).sum(),
        }
    }
}
impl std::ops::SubAssign for OneYearResources {
    fn sub_assign(&mut self, other: Self) {
        self.c -= other.c
    }
}
impl std::ops::Add for OneYearResources {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        OneYearResources {
            c: self.c + other.c,
        }
    }
}
/**
Dividing out the unit gives back a scalar
 */
impl Div<f64> for OneYearResources {
    type Output = Self;
    fn div(self, f: f64) -> Self {
        OneYearResources { c: self.c / f }
    }
}
/**
To get a feeling for the number of year-resources, transform it into kCal (at 2000 kCal a day, 365 days a year).
 */
impl std::fmt::Debug for OneYearResources {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:.2} ({:.2} kCal)", self.c, self.c * 2000. * 365.)
    }
}

use std::num::ParseFloatError;
use std::str::FromStr;
impl FromStr for OneYearResources {
    type Err = ParseFloatError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(OneYearResources {
            c: f64::from_str(s)?,
        })
    }
}

mod array_serde;
use array_serde::BigArray;
pub const ATTESTED_ECOREGIONS: usize = 303;

use serde_derive::{Deserialize, Serialize};
#[derive(Clone, Copy, Serialize, Deserialize)]
pub struct Ecovector {
    #[serde(with = "BigArray")]
    pub entries: [f64; ATTESTED_ECOREGIONS + 1],
    pub minimum: f64,
}

impl Mul<f64> for Ecovector {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        #[allow(clippy::suspicious_arithmetic_impl)]
        let mut result = [0.; ATTESTED_ECOREGIONS + 1];
        for (i, item) in self.entries.iter().enumerate() {
            result[i] = f64::max(self.minimum, item * rhs); //TODO: Move this part of the logic elsewhere
        }
        Ecovector {
            entries: result,
            minimum: self.minimum,
        }
    }
}

use std::slice::Iter;
impl<'a> IntoIterator for &'a Ecovector {
    type Item = &'a f64;
    type IntoIter = Iter<'a, f64>;
    fn into_iter(self) -> Self::IntoIter {
        self.entries.iter()
    }
}

impl Index<usize> for Ecovector {
    type Output = f64;

    fn index(&self, index: usize) -> &f64 {
        &self.entries[index]
    }
}

impl Default for Ecovector {
    fn default() -> Self {
        Ecovector {
            entries: [0.5; ATTESTED_ECOREGIONS + 1],
            minimum: 0.5,
        }
    }
}

impl From<f64> for Ecovector {
    fn from(min: f64) -> Self {
        Ecovector {
            entries: [min; ATTESTED_ECOREGIONS + 1],
            minimum: min,
        }
    }
}
