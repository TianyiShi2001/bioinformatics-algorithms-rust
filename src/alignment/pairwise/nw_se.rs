//! Space-efficient version of Needleman-Wunsch's algorithm

use crate::alignment::*;
use std::cmp::max;
use std::mem;

pub struct NwSpaceEfficientAligner<F: MatchFunc> {
    scoring: Scoring<F>,
}

impl<F: MatchFunc> NwSpaceEfficientAligner<F> {
    pub fn new(scoring: Scoring<F>) -> Self {
        NwSpaceEfficientAligner { scoring }
    }
    pub fn global<'a>(&self, x: &'a Seq, y: &'a Seq) -> AlignmentResult<'a> {
        let alignment = self.compute_recursive(x, y, x.len(), y.len());
        let score = self.cost_only(x, y, false)[y.len()];
        return AlignmentResult {
            alignment,
            score,
            x,
            y,
            xstart: 0,
            ystart: 0,
            xend: x.len(),
            yend: y.len(),
        };
    }

    fn compute_recursive(&self, x: &Seq, y: &Seq, m: usize, n: usize) -> Vec<AlignmentOperation> {
        if n == 0 {
            return vec![AlignmentOperation::Del; m];
        }
        if m == 1 {
            return self.nw_onerow(x[0], y, n);
        }
        let (imid, jmid) = self.find_mid(x, y, m, n);

        return [
            self.compute_recursive(&x[..imid], &y[..jmid], imid, jmid),
            self.compute_recursive(&x[imid..], &y[jmid..], m - imid, n - jmid),
        ]
        .concat();
    }
    fn find_mid(&self, x: &Seq, y: &Seq, m: usize, n: usize) -> (usize, usize) {
        let imid = m / 2;
        let upper = self.cost_only(&x[..imid], y, false);
        let lower = self.cost_only(&x[imid..], y, true);
        let mut max = Score::MIN;
        let mut jmid = 0;
        for j in 0..=n {
            let x = upper[j] + lower[n - j];
            if x > max {
                max = x;
                jmid = j;
            }
        }
        (imid, jmid)
    }

    fn cost_only_0(&self, x: &Seq, y: &Seq, rev: bool) -> Vec<Score> {
        let m = x.len() + 1;
        let n = y.len() + 1;
        let mut prev: Vec<Score> = vec![0; n];
        let mut curr: Vec<Score> = vec![0; n];
        for j in 1..n {
            curr[j] = curr[j - 1] + self.scoring.gap_extend; // 0th row
        }
        for i in 1..m {
            mem::swap(&mut prev, &mut curr);
            curr[0] = prev[0] + self.scoring.gap_extend;
            for j in 1..n {
                let (xi, yj) = if rev {
                    (x[m - i - 1], y[n - j - 1])
                } else {
                    (x[i - 1], y[j - 1])
                };
                curr[j] = self
                    .scoring
                    .max_score(prev[j], curr[j - 1], prev[j - 1], xi, yj)
            }
        }
        curr
    }
    fn cost_only_1(&self, x: &Seq, y: &Seq, rev: bool) -> Vec<Score> {
        let m = x.len() + 1;
        let n = y.len() + 1;
        let mut v0: Vec<Score> = vec![0; n];
        let mut p = &mut v0 as *mut Vec<Score>;
        let mut v1: Vec<Score> = vec![0; n];
        let mut c = &mut v1 as *mut Vec<Score>;
        let mut tmp: *mut Vec<Score>; // for swapping
        for j in 1..n {
            v1[j] = v1[j - 1] + self.scoring.gap_extend; // 0th row
        }
        for i in 1..m {
            unsafe {
                tmp = p;
                p = c;
                c = tmp;
                (*c)[0] = (*p)[0] + self.scoring.gap_extend;
                for j in 1..n {
                    let (xi, yj) = if rev {
                        (x[m - i - 1], y[n - j - 1])
                    } else {
                        (x[i - 1], y[j - 1])
                    };
                    (*c)[j] = self
                        .scoring
                        .max_score((*p)[j], (*c)[j - 1], (*p)[j - 1], xi, yj)
                }
            }
        }
        if m % 2 == 0 {
            v0
        } else {
            v1
        }
    }
    /// Cost-only NW with only one vector and one scalar
    fn cost_only(&self, x: &Seq, y: &Seq, rev: bool) -> Vec<Score> {
        let m = x.len() + 1;
        let n = y.len() + 1;
        let mut cc: Vec<Score> = vec![0; n];
        let mut s: Score;
        for j in 1..n {
            cc[j] = cc[j - 1] + self.scoring.gap_extend; // 0th row
        }
        for i in 1..m {
            s = cc[0];
            cc[0] = s + self.scoring.gap_extend;
            for j in 1..n {
                // ---old cc[j]---
                let up = cc[j] + self.scoring.gap_extend; // update UP
                let left = cc[j - 1] + self.scoring.gap_extend; // update LEFT
                let diag = if rev {
                    // compute DIAG
                    s + self.scoring.match_fn.score(x[m - i - 1], y[n - j - 1])
                } else {
                    s + self.scoring.match_fn.score(x[i - 1], y[j - 1])
                };
                s = cc[j]; // old cc[j] becomes old cc[j-1], i.e. `s`, for the next round
                cc[j] = max(up, max(left, diag)); // ---new cc[j]---
            }
        }
        cc
    }
    fn nw_onerow(&self, x: u8, y: &Seq, n: usize) -> Vec<AlignmentOperation> {
        let mut S = Vec::<Score>::with_capacity(n + 1);
        let mut T = Vec::<AlignmentOperation>::with_capacity(n + 1);
        S.push(self.scoring.gap_extend);
        for j in 0..n {
            let up = (j as Score + 2) * self.scoring.gap_extend;
            let left = S[j] + self.scoring.gap_extend;
            let (mut diag, diag_operation) = self.scoring.match_fn.score_with_operation(x, y[j]);
            diag += (j as Score) * self.scoring.gap_extend;
            let (max, dir) = max_score_and_operation_precomputed(up, left, diag, diag_operation);
            S.push(max);
            T.push(dir);
        }
        let mut tmp = AlignmentOperation::None;
        match T.iter().rposition(|&x| {
            tmp = x;
            x == AlignmentOperation::Match || x == AlignmentOperation::Subst
        }) {
            Some(n) => {
                for j in 0..n {
                    T[j] = AlignmentOperation::Ins;
                }
                T[n] = tmp
            }
            None => {
                for j in 0..n {
                    T[j] = AlignmentOperation::Ins
                }
                T.push(AlignmentOperation::Del);
            }
        }
        T
    }
}
