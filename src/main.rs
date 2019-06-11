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

    let za = 1.8;
    let zb = 2.8;
    let ra = [0.0, 0.0, 0.0];
    let rb = [0.5, 0.8, -0.2];
    // println!(
    //     "{}",
    //     integrals::get_overlap(za, zb, ra, rb, [2, 1, 0, 1, 1, 0])
    // );

    let coords = nd::arr2(&[
        [0.000000000000, -0.143225816552, 0.000000000000],
        [1.638036840407, 1.136548822547, 0.000000000000],
        [-1.638036840407, 1.136548822547, 0.000000000000],
    ]);
    // println!("{}", coords);

    let a = basis::PGTO::new(
        [coords[[0, 0]], coords[[0, 1]], coords[[0, 2]]],
        1.0,
        [0, 0, 0],
    );
    // println!("{:?}", a);
    // println!("{}", basis::S(&a, &a));

    let gil = Python::acquire_gil();
    let elements = vec![8, 1, 1];
    let bseresult = basis::get_bse_json(gil.python(), "STO-3G", elements);
    // println!("{:#?}", bseresult);
}
