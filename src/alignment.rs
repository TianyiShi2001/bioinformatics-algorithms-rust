pub mod pairwise;

#[derive(Debug, Default)]
pub struct AlignmentResult<'a> {
    pub alignment: Vec<Direction>,
    pub score: i32,
    pub x: &'a [u8],
    pub y: &'a [u8],
    pub i: usize,
    pub j: usize,
}

impl<'a> AlignmentResult<'a> {
    pub fn as_strings(&self, gap_char: char) -> (String, String) {
        let gap_char = gap_char as u8;
        let mut x: Vec<u8> = Vec::with_capacity(self.alignment.len());
        let mut y: Vec<u8> = Vec::with_capacity(self.alignment.len());
        let mut i = self.i;
        let mut j = self.j;
        for dir in &self.alignment {
            match dir {
                Direction::Del => {
                    x.push(self.x[i]);
                    y.push(gap_char);
                    i += 1;
                }
                Direction::Ins => {
                    x.push(gap_char);
                    y.push(self.y[j]);
                    j += 1;
                }
                Direction::Sub => {
                    x.push(self.x[i]);
                    y.push(self.y[j]);
                    i += 1;
                    j += 1;
                }
                Direction::None => {}
            }
        }
        unsafe {
            let x = String::from_utf8_unchecked(x);
            let y = String::from_utf8_unchecked(y);
            (x, y)
        }
    }
}

/// Trait required to instantiate a Scoring instance
pub trait MatchFunc {
    fn score(&self, a: u8, b: u8) -> i32;
}

/// A concrete data structure which implements trait MatchFunc with constant
/// match and mismatch scores
#[derive(Debug, Clone)]
pub struct MatchParams {
    pub match_score: i32,
    pub mismatch_score: i32,
}

impl MatchParams {
    /// Create new MatchParams instance with given match and mismatch scores
    ///
    /// # Arguments
    ///
    /// * `match_score` - the score for a match (should not be negative)
    /// * `mismatch_score` - the score for a mismatch (should not be positive)
    pub fn new(match_score: i32, mismatch_score: i32) -> Self {
        assert!(match_score >= 0, "match_score can't be negative");
        assert!(mismatch_score <= 0, "mismatch_score can't be positive");
        MatchParams {
            match_score,
            mismatch_score,
        }
    }
}

impl MatchFunc for MatchParams {
    #[inline]
    fn score(&self, a: u8, b: u8) -> i32 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// The trait Matchfunc is also implemented for Fn(u8, u8) -> i32 so that Scoring
/// can be instantiated using closures and custom user defined functions
impl<F> MatchFunc for F
where
    F: Fn(u8, u8) -> i32,
{
    fn score(&self, a: u8, b: u8) -> i32 {
        (self)(a, b)
    }
}

/// Details of scoring are encapsulated in this structure.
/// An affine gap score model is used so that the gap score for a length 'k' is:
/// GapScore(k) = gap_open + gap_extend * k
#[derive(Debug, Clone)]
pub struct Scoring<F: MatchFunc> {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub match_fn: F,
    pub match_scores: Option<(i32, i32)>,
}

impl Scoring<MatchParams> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// match and mismatch scores.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_score` - the score for a match
    /// * `mismatch_score` - the score for a mismatch
    pub fn from_scores(
        gap_open: i32,
        gap_extend: i32,
        match_score: i32,
        mismatch_score: i32,
    ) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn: MatchParams::new(match_score, mismatch_score),
            match_scores: Some((match_score, mismatch_score)),
        }
    }
}

impl<F: MatchFunc> Scoring<F> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// and the score function. The clip penalties are set to MIN_SCORE by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions (also see bio::scores)
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn,
            match_scores: None,
        }
    }
}

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Direction {
    Del, // up
    Ins, // left
    Sub, // diagonal
    None,
}

pub type Score = i32;
pub type Matrix<T> = Vec<Vec<T>>;
pub type ScoreMatrix = Matrix<Score>;
pub type TracebackMatrix = Matrix<Direction>;
pub type Seq = [u8];
pub type Coords = (usize, usize);

pub fn compute_max_score_and_direction(up: Score, left: Score, diag: Score) -> (Score, Direction) {
    let mut max = up;
    let mut dir = Direction::Del;
    if left > max {
        max = left;
        dir = Direction::Ins;
    }
    if diag > max {
        max = diag;
        dir = Direction::Sub;
    }
    (max, dir)
}
