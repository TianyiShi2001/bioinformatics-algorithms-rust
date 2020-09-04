//! Alignment with affine gap penalty in linear space, by combining Gotoh's (1982) and
//! Hirschberg's (1975) ideas, which was first implemented in C (Myers & Miller 1988).
//!
//! # Time Complexity
//!
//! O(n * m) for strings of length m and n.
//!
//! # Space Complexity
//!
//! The space usage depends on the `cost_only` method of [GotohSpaceEfficientAligner](struct.GotohSpaceEfficientAligner),
//! which uses 6 scalars and 2 vectors of length (n + 1), where n is the length of the shorter sequence.
//! [See also](struct.GotohSpaceEfficientAligner.html#space-complexity)
//!
//! # Example
//!
//! ```
//! use bioinformatics_algorithms::alignment::pairwise::gotoh_space_efficient::GotohSpaceEfficientAligner;
//! use bioinformatics_algorithms::alignment::Scoring;
//! let x = b"ATGATGATG";
//! let y = b"ATGAATG";
//! let mut aligner = GotohSpaceEfficientAligner::new(Scoring::from_scores(-5, -1, 2, -1));
//! let res = aligner.global(x, y);
//! ```
//!
//! # References
//!
//! - [Eugene W. Myers and Webb Miller (1988) Optimal alignments in linear space. _Bioinformatics_ **4**: 11-17.](https://doi.org/10.1093/bioinformatics/4.1.11)
//! - [Hirschberg, D. S. (1975) A linear space algorithm for computing maximal common subsequences. _Commun. Assoc. Comput. Mach._ **18**: 341-343.](https://doi.org/10.1145/360825.360861)
//! - [Gotoh, O. (1982) An improved algorithm for matching biological sequences. _J. Molec. Biol._ **162**: 705-708.](https://doi.org/10.1016/0022-2836(82)90398-9)

use crate::alignment::*;
use std::cmp::max;

pub struct GotohSpaceEfficientAligner<F: MatchFunc> {
    scoring: Scoring<F>,
}

impl<F: MatchFunc> GotohSpaceEfficientAligner<F> {
    pub fn new(scoring: Scoring<F>) -> Self {
        GotohSpaceEfficientAligner { scoring }
    }
    pub fn global<'a>(&mut self, x: &'a Seq, y: &'a Seq) -> AlignmentResult<'a> {
        let alignment = self.nw_recursive(
            x,
            y,
            x.len(),
            y.len(),
            self.scoring.gap_open,
            self.scoring.gap_open,
        );
        let score = self.cost_only(x, y, false, self.scoring.gap_open).0[y.len()];
        return AlignmentResult {
            alignment,
            score,
            x,
            y,
            i: 0,
            j: 0,
        };
    }
    /// Recursively compute alignments of sub-sequences and concatenating them
    fn nw_recursive(
        &self,
        x: &Seq,
        y: &Seq,
        m: usize,
        n: usize,
        tb: Score,
        te: Score,
    ) -> Vec<Direction> {
        // * m = x.len(); n = y.len()
        if n == 0 {
            return vec![Direction::Del; m];
        }
        if m == 0 {
            return vec![Direction::Ins; n];
        }
        if m == 1 {
            return self.nw_onerow(x[0], y, n, tb, te);
        }
        let (imid, jmid, join_by_deletion) = self.find_mid(x, y, m, n, tb, te);
        return if join_by_deletion {
            [
                self.nw_recursive(&x[..imid - 1], &y[..jmid], imid - 1, jmid, tb, 0),
                vec![Direction::Del; 2],
                self.nw_recursive(&x[imid + 1..], &y[jmid..], m - imid - 1, n - jmid, 0, te),
            ]
            .concat()
        } else {
            [
                self.nw_recursive(
                    &x[..imid],
                    &y[..jmid],
                    imid,
                    jmid,
                    tb,
                    self.scoring.gap_open,
                ),
                self.nw_recursive(
                    &x[imid..],
                    &y[jmid..],
                    m - imid,
                    n - jmid,
                    0,
                    self.scoring.gap_open,
                ),
            ]
            .concat()
        };
    }

    fn find_mid(
        &self,
        x: &Seq,
        y: &Seq,
        m: usize,
        n: usize,
        tb: Score,
        te: Score,
    ) -> (usize, usize, bool) {
        let imid = m / 2;
        let (cc_upper, dd_upper) = self.cost_only(&x[..imid], y, false, tb);
        let (cc_lower, dd_lower) = self.cost_only(&x[imid..], y, true, te);
        let mut max = Score::MIN;
        let mut jmid = 0;
        let mut join_by_deletion = false;
        for j in 0..=n {
            let c = cc_upper[j] + cc_lower[n - j];
            if c > max {
                max = c;
                jmid = j;
                join_by_deletion = false;
            }
            let d = dd_upper[j] + dd_lower[n - j] - self.scoring.gap_open; // subtract duplicating open!
            if d > max {
                max = d;
                jmid = j;
                join_by_deletion = true;
            }
        }
        (imid, jmid, join_by_deletion)
    }

    /// Cost-only (score-only) Gotoh's algorithm in linear space
    /// # Space Complexity
    /// Use six scalars and two vectors of length (N + 1), where N is the length
    /// of the shorter sequence.
    fn cost_only(&self, x: &Seq, y: &Seq, rev: bool, tx: Score) -> (Vec<Score>, Vec<Score>) {
        let m = x.len() + 1;
        let n = y.len() + 1;
        let mut cc: Vec<Score> = vec![0; n]; // match/mismatch
        let mut dd: Vec<Score> = vec![0; n]; // deletion
        let mut e: Score; // I(i, j-1)
        let mut c: Score; // C(i, j-1)
        let mut s: Score; // C(i-1, j-1)
        let mut t: Score;
        t = self.scoring.gap_open; // originally self.scoring.gap_open;
        for j in 1..n {
            t += self.scoring.gap_extend;
            cc[j] = t;
            dd[j] = Score::MIN;
        }
        t = tx; // originally self.scoring.gap_open;
        for i in 1..m {
            s = cc[0];
            t += self.scoring.gap_extend;
            c = t;
            cc[0] = c;
            e = Score::MIN;
            for j in 1..n {
                e = max(e, c + self.scoring.gap_open) + self.scoring.gap_extend; // update e to I[i,j]
                dd[j] = max(dd[j], cc[j] + self.scoring.gap_open) + self.scoring.gap_extend; // cc[j] = C[i-1, j]
                c = if rev {
                    max(
                        max(dd[j], e),
                        s + self.scoring.match_fn.score(x[m - i - 1], y[n - j - 1]),
                    )
                } else {
                    max(
                        max(dd[j], e),
                        s + self.scoring.match_fn.score(x[i - 1], y[j - 1]),
                    )
                };
                s = cc[j];
                cc[j] = c;
            }
        }
        (cc, dd)
    }
    fn nw_onerow(&self, x: u8, y: &Seq, n: usize, tb: Score, te: Score) -> Vec<Direction> {
        let score_by_indels_only =
            max(tb, te) + self.scoring.gap_extend * (n as Score + 1) + self.scoring.gap_open;
        let mut max = score_by_indels_only;
        let score_with_one_substitution_base =
            (n as Score - 1) * self.scoring.gap_extend + self.scoring.gap_open; // plus substitution score and possibly one more gap_open
        let mut maxj_ = 0usize;
        for j_ in 0..n {
            // index of sequence instead of matrix; y[j] instead of j[j-1] is the jth character
            let score = score_with_one_substitution_base
                + self.scoring.match_fn.score(x, y[j_])
                + if j_ == 0 || j_ == n - 1 {
                    0
                } else {
                    self.scoring.gap_open
                };
            if score > max {
                max = score;
                maxj_ = j_;
            }
        }
        return if max == score_by_indels_only {
            let mut res = Vec::with_capacity(n + 1);
            res.push(Direction::Del);
            for _j in 0..n {
                res.push(Direction::Ins)
            }
            res
        } else {
            let mut res = Vec::with_capacity(n);
            for _j in 0..maxj_ {
                res.push(Direction::Ins)
            }
            res.push(Direction::Sub);
            for _j in 0..(n - maxj_ - 1) {
                res.push(Direction::Ins)
            }
            res
        };
    }
}
