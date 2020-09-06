use bioinformatics_algorithms::seqanalysis::protein::count_aa;
use bioinformatics_algorithms::seqanalysis::protein::isoelectric_point;
fn main() {
    let res = isoelectric_point(b"MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV");
    println!("{:?}", res);
}
