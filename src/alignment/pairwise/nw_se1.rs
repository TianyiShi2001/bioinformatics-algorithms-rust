//! Space-efficient version of Needleman-Wunsch's algorithm

use crate::alignment::*;
use std::cmp::max;
use std::mem;
pub struct NwSpaceEfficientAligner<F: MatchFunc> {
    scoring: Scoring<F>,
}
use std::ptr;

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
            println!("{}", "n=0!");
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
    fn cost_only(&self, x: &Seq, y: &Seq, rev: bool) -> Vec<Score> {
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
                let up = prev[j] + self.scoring.gap_extend;
                let left = curr[j - 1] + self.scoring.gap_extend;
                let diag = if rev {
                    prev[j - 1] + self.scoring.match_fn.score(x[m - i - 1], y[n - j - 1])
                } else {
                    prev[j - 1] + self.scoring.match_fn.score(x[i - 1], y[j - 1])
                };
                curr[j] = max(up, max(left, diag));
            }
        }
        curr
    }
    fn nw_onerow(&self, x: u8, y: &Seq, n: usize) -> Vec<AlignmentOperation> {
        let mut S = Vec::<Score>::with_capacity(n + 1);
        let mut T = Vec::<AlignmentOperation>::with_capacity(n + 1);
        S.push(self.scoring.gap_extend);
        for j in 0..n {
            let up = (j as Score + 2) * self.scoring.gap_extend;
            let left = S[j] + self.scoring.gap_extend;
            let diag =
                (j as Score) * self.scoring.gap_extend + self.scoring.match_fn.score(x, y[j]);
            let (max, dir) = compute_max_score_and_operation(up, left, diag);
            S.push(max);
            T.push(dir);
        }
        match T.iter().rposition(|&x| x == AlignmentOperation::Sub) {
            Some(n) => {
                for j in 0..n {
                    T[j] = AlignmentOperation::Ins;
                }
                T[n] = AlignmentOperation::Sub
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
