use bioinformatics_algorithms::alignment::pairwise::gotoh_space_efficient::GotohSpaceEfficientAligner;
use bioinformatics_algorithms::alignment::Scoring;

pub fn test1() {
    let x = b"ATGATGATG";
    let y = b"ATGAATG";
    let mut aligner = GotohSpaceEfficientAligner::new(Scoring::from_scores(-5, -1, 2, -1));
    let res = aligner.global(x, y);
    println!("{:?}", res.as_strings('-'));
}
