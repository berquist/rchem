use chemfiles;

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

    let mut trajectory = chemfiles::Trajectory::open("water.xyz", 'r').unwrap();
    let mut frame = chemfiles::Frame::new();
    trajectory.read(&mut frame).unwrap();
    let natom = frame.size();
    let atomcoords = frame.positions();
    let atomnos: Vec<_> = (0..natom).map(|i| frame.atom(i).atomic_number()).collect();
    let basis_set = basis::Basis::new(&atomnos, &atomcoords, "STO-3G");
}
