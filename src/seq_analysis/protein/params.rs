// use crate::seq_analysis::protein::AminoAcidCount;
// use bio::utils::TextSlice;
use lazy_static::lazy_static;
use std::collections::BTreeMap;

pub mod isoelectric_point {
    use super::*;
    const N_TERM_PKA_DEFAULT: f32 = 7.5;
    const C_TERM_PKA_DEFAULT: f32 = 3.55;
    pub enum Charge {
        Positive,
        Negative,
    }
    lazy_static! {
        pub static ref pKa_table: BTreeMap<u8, (f32, Charge)> = {
            let mut m = BTreeMap::new();
            m.insert(b'K', (10.0, Charge::Positive));
            m.insert(b'R', (12.0, Charge::Positive));
            m.insert(b'H', (5.98, Charge::Positive));
            m.insert(b'D', (4.05, Charge::Negative));
            m.insert(b'E', (4.45, Charge::Negative));
            m.insert(b'C', (9.00, Charge::Negative));
            m.insert(b'Y', (10.0, Charge::Negative));
            m
        };
        pub static ref n_terminal_pKa_table: BTreeMap<u8, f32> = {
            let mut m = BTreeMap::new();
            m.insert(b'A', 7.59);
            m.insert(b'M', 7.00);
            m.insert(b'S', 6.93);
            m.insert(b'P', 8.36);
            m.insert(b'T', 6.82);
            m.insert(b'V', 7.44);
            m.insert(b'E', 7.70);
            m
        };
        pub static ref c_terminal_pKa_table: BTreeMap<u8, f32> = {
            let mut m = BTreeMap::new();
            m.insert(b'D', 4.55);
            m.insert(b'E', 4.75);
            m
        };
    }
}
