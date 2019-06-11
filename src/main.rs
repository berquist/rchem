use chemfiles;
use cpython::Python;
use ndarray as nd;

mod basis;
mod integrals;
mod shell;

fn main() {
    // http://www.patorjk.com/software/taag/#p=display&f=3D%20Diagonal&t=rchem
    let logo = r#"
                   ,---,                       ____
                 ,--.' |                     ,'  , `.
  __  ,-.        |  |  :                  ,-+-,.' _ |
,' ,'/ /|        :  :  :               ,-+-. ;   , ||
'  | |' | ,---.  :  |  |,--.   ,---.  ,--.'|'   |  ||
|  |   ,'/     \ |  :  '   |  /     \|   |  ,', |  |,
'  :  / /    / ' |  |   /' : /    /  |   | /  | |--'
|  | ' .    ' /  '  :  | | |.    ' / |   : |  | ,
;  : | '   ; :__ |  |  ' | :'   ;   /|   : |  |/
|  , ; '   | '.'||  :  :_:,''   |  / |   | |`-'
 ---'  |   :    :|  | ,'    |   :    |   ;/
        \   \  / `--''       \   \  /'---'
         `----'               `----'
"#;
    println!("{}", logo);

    // let za = 1.8;
    // let zb = 2.8;
    // let ra = [0.0, 0.0, 0.0];
    // let rb = [0.5, 0.8, -0.2];
    // println!(
    //     "{}",
    //     integrals::get_overlap(za, zb, ra, rb, [2, 1, 0, 1, 1, 0])
    // );

    let mut trajectory = chemfiles::Trajectory::open("water.xyz", 'r').unwrap();
    let mut frame = chemfiles::Frame::new();
    trajectory.read(&mut frame).unwrap();
    let natom = frame.size();
    let atomcoords = frame.positions();
    let atomnos: Vec<_> = (0..natom).map(|i| frame.atom(i).atomic_number()).collect();
    println!("{:?}", atomcoords);
    println!("{:?}", atomnos);

    let a = basis::PGTO::new(
        [atomcoords[0][0], atomcoords[0][1], atomcoords[0][2]],
        1.0,
        [0, 0, 0],
    );
    // println!("{:?}", a);
    // println!("{}", basis::S(&a, &a));

    let gil = Python::acquire_gil();
    let bseresult = basis::get_bse_json(gil.python(), "STO-3G", atomnos);
    println!("{:#?}", bseresult);
}
