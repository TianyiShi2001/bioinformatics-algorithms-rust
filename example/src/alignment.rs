use bioinformatics_algorithms::alignment::pairwise::gotoh_space_efficient::GotohSpaceEfficientAligner;
use bioinformatics_algorithms::alignment::pairwise::nw_se::NwSpaceEfficientAligner;
use bioinformatics_algorithms::alignment::AlignmentResult;
use bioinformatics_algorithms::alignment::MatchParams;
use bioinformatics_algorithms::alignment::Scoring;

use lazy_static::lazy_static;

lazy_static! {
    static ref ALIGNER_GOTOH_SPACE_EFFICIENT: GotohSpaceEfficientAligner<MatchParams> =
        GotohSpaceEfficientAligner::new(Scoring::from_scores(-5, -1, 2, -1));
    static ref ALIGNER_NW_SPACE_EFFICIENT: NwSpaceEfficientAligner<MatchParams> =
        NwSpaceEfficientAligner::new(Scoring::from_scores(0, -1, 2, -1));
}

static X1: &[u8] = b"ATGATGATGATGATGATGATGCG"; // ATGATGATGATGATGATGATGCG
static Y1: &[u8] = b"ATGAATGCG"; //               ATGA--------------ATGCG (14G, 9M)
static X2: &[u8] = b"AAAAAAAGGGTTTCCCCCCCCCC"; // AAAAAAAGGGTTTCCCCCCCCCC
static Y2: &[u8] = b"AAAAGGGTTT"; //              ---AAAAGGGTTT----------

pub fn show_alignment_result(res: AlignmentResult) {
    println!("Score, {}; {:?}", res.score, res.as_strings('-'));
}

pub fn test1() {
    show_alignment_result(ALIGNER_GOTOH_SPACE_EFFICIENT.global(X1, Y1));
    show_alignment_result(ALIGNER_NW_SPACE_EFFICIENT.global(X1, Y1));
}

pub fn test2() {
    show_alignment_result(ALIGNER_GOTOH_SPACE_EFFICIENT.global(X2, Y2));
    show_alignment_result(ALIGNER_NW_SPACE_EFFICIENT.global(X2, Y2));
}
