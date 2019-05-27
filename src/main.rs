mod integrals;

fn main() {
    let za = 1.8;
    let zb = 2.8;
    let ra = [0.0, 0.0, 0.0];
    let rb = [0.5, 0.8, -0.2];
    println!("{}", integrals::get_overlap(za, zb, ra, rb, [2, 1, 0, 1, 1, 0]));
}
