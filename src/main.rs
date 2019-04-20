use std::f64::consts::PI;

fn main() {
    println!("Hello, world!");
}

enum X2kind {
    S,  // overlap
    T,  // kinetic energy
    V,  // nuclear-electron attraction
    M,  // multipole moment
    L,  // angular momentum (not implemented)
    E,  // electric field (not implemented)
    J,  // spin-orbit (not implemented)
    FC, // Fermi contact (not implemented)
}

struct X2 {
    scale: f64,
    prefactors: Vec<usize>,
    q: [u8; 6],
    kind: X2kind,
    operator: [u8; 3],
    d: u8,
    order: u8,
}

struct X4 {
    scale: f64,
    prefactors: Vec<f64>,
    q: [u8; 12],
    order: u8,
}

fn get_coulomb(
    za: f64,
    zb: f64,
    zc: f64,
    zd: f64,
    ra: [f64; 3],
    rb: [f64; 3],
    rc: [f64; 3],
    rd: [f64; 3],
    c: [u8; 12],
) -> f64 {
    return 0.0;
}

fn get_bi_center(z1: f64, z2: f64, r1: [f64; 3], r2: [f64; 3]) -> [f64; 3] {
    let z = z1 + z2;
    let rx = (z1 * r1[0] + z2 * r2[0]) / z;
    let ry = (z1 * r1[1] + z2 * r2[1]) / z;
    let rz = (z1 * r1[2] + z2 * r2[2]) / z;
    return [rx, ry, rz];
}

fn get_r12_squared(r1: [f64; 3], r2: [f64; 3]) -> f64 {
    return (r1[0] - r2[0]).powi(2) + (r1[1] - r2[1]).powi(2) + (r1[2] - r2[2]).powi(2);
}

fn find_fun_to_lower(q: Vec<u8>) -> usize {
    // Determine the total angular momentum on each center.
    let n = q.len() / 3;
    let mut l = vec!();
    for i in 0..n {
        l.push(q[i * 3] + q[i * 3 + 1] + q[i * 3 + 2])
    }

    // find function to lower
    // start with lowest angular momentum above s
    let mut fun = 0;
    let mut kmax = 1 + *l.iter().max().expect("There are no centers?");

    for i in 0..n {
        let k = l[i];
        // If we're larger than a s-function...
        if k > 0 {
            if k < kmax {
                kmax = k;
                fun = i;
            }
        }
    }
    return fun
}

fn apply_os2(mut x: X2, kind: X2kind) -> Vec<X2> {
    let orders = vec![x.q[0], x.q[1], x.q[2], x.q[3], x.q[4], x.q[5], x.operator[0], x.operator[1], x.operator[2]];

    // base case
    let order_sum: u8 = orders.iter().sum();
    if order_sum == 0 {
        x.kind = kind;
        return vec!(x);
    }

    // Determine which basis function and component to lower.
    // The component is one of (x, y, z).
    let fun = find_fun_to_lower(orders);

    return vec!();
}

fn get_overlap(za: f64, zb: f64, ra: [f64; 3], rb: [f64; 3], c: [u8; 6]) -> f64 {
    let z = za + zb;
    let e = za * zb / (za + zb);
    let rp = get_bi_center(za, zb, ra, rb);
    let ab = get_r12_squared(ra, rb);
    let aux = (-e * ab).exp() * (PI / z).powf(1.5);

    let prefac = vec![
        rp[0] - ra[0],
        rp[0] - rb[0],
        rp[1] - ra[1],
        rp[1] - rb[1],
        rp[2] - ra[2],
        rp[2] - rb[2],
        0.5 / z,
        0.5 / z,
    ];

    let fun = X2 {
        scale: 1.0,
        prefactors: vec!(),
        q: c,
        kind: X2kind::S,
        operator: [0, 0, 0],
        d: 0,
        order: 0,
    };
    let expansion = apply_os2(fun, X2kind::S);
    let mut integral = 0.0;
    for f in expansion.iter() {
        let mut g = 1.0;
        for k in f.prefactors.iter() {
            g *= prefac[*k];
        }
        integral += f.scale * aux * g;
    }
    return integral;
}

fn get_kinetic(za: f64, zb: f64, ra: [f64; 3], rb: [f64; 3], c: [u8; 6]) -> f64 {
    return 0.0;
}

fn get_nuclear(za: f64, zb: f64, ra: [f64; 3], rb: [f64; 3], rc: [f64; 3], c: [u8; 6]) -> f64 {
    return 0.0;
}

fn get_moment(
    za: f64,
    zb: f64,
    ra: [f64; 3],
    rb: [f64; 3],
    rc: [f64; 3],
    c: [u8; 6],
    order: [u8; 3],
) -> f64 {
    return 0.0;
}

#[cfg(test)]
mod tests {
    use super::find_fun_to_lower;
    use super::get_coulomb;
    use super::get_kinetic;
    use super::get_moment;
    use super::get_nuclear;
    use super::get_overlap;

    #[test]
    fn test_find_fun_to_lower() {
        assert_eq!(find_fun_to_lower([1, 0, 0, 0, 0, 0].to_vec()), 0);
        assert_eq!(find_fun_to_lower([0, 1, 0, 0, 0, 0].to_vec()), 0);
        assert_eq!(find_fun_to_lower([0, 0, 1, 0, 0, 0].to_vec()), 0);
        assert_eq!(find_fun_to_lower([0, 0, 0, 1, 0, 0].to_vec()), 1);
        assert_eq!(find_fun_to_lower([0, 0, 0, 0, 1, 0].to_vec()), 1);
        assert_eq!(find_fun_to_lower([0, 0, 0, 0, 0, 1].to_vec()), 1);
        assert_eq!(find_fun_to_lower([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0].to_vec()), 1);
        assert_eq!(find_fun_to_lower([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1].to_vec()), 3);
        assert_eq!(find_fun_to_lower([1, 0, 0, 0, 0, 0, 0, 0, 1].to_vec()), 0);
        assert_eq!(find_fun_to_lower([0, 0, 0, 0, 0, 0, 0, 0, 1].to_vec()), 2);
        assert_eq!(find_fun_to_lower([0, 0, 0, 0, 1, 0, 0, 0, 1].to_vec()), 1);
        assert_eq!(find_fun_to_lower([0, 0, 0, 0, 2, 0, 0, 0, 1].to_vec()), 2);
    }

    #[test]
    fn test_get_coulomb() {
        let za = 1.1;
        let zb = 1.2;
        let zc = 1.3;
        let zd = 1.4;

        let ra = [1.0, 0.0, 1.0];
        let rb = [0.0, 1.0, 2.0];
        let rc = [0.0, 0.0, 3.0];
        let rd = [0.0, 0.0, 4.0];

        let thresh = 1.0e-16;

        let reference = 1.71817807954e-05;

        let integral = get_coulomb(
            za,
            zb,
            zc,
            zd,
            ra,
            rb,
            rc,
            rd,
            [2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0],
        );

        assert!((integral - reference).abs() < thresh);
    }

    #[test]
    fn test_get_overlap() {
        let za = 1.8;
        let zb = 2.8;
        let ra = [0.0, 0.0, 0.0];
        let rb = [0.5, 0.8, -0.2];

        let thresh = 1.0e-16;

        let integral = get_overlap(za, zb, ra, rb, [0, 0, 0, 0, 0, 0]);
        assert!((integral - 0.20373275913014607).abs() < thresh);

        let integral = get_overlap(za, zb, ra, rb, [1, 0, 0, 0, 0, 0]);
        assert!((integral - 0.062005622343957505).abs() < thresh);

        let integral = get_overlap(za, zb, ra, rb, [1, 1, 0, 1, 1, 0]);
        assert!((integral - -0.00043801221837779696).abs() < thresh);

        let integral = get_overlap(za, zb, ra, rb, [2, 1, 0, 1, 1, 0]);
        assert!((integral - -0.0002385994651113168).abs() < thresh);
    }

    #[test]
    fn test_get_kinetic() {
        let za = 1.8;
        let zb = 2.0;
        let ra = [0.0, 0.0, 0.0];
        let rb = [0.5, 0.8, -0.2];

        let thresh = 1.0e-16;

        let integral = get_kinetic(za, zb, ra, rb, [0, 0, 0, 0, 0, 0]);
        assert!((integral - 0.3652714583525358).abs() < thresh);

        let integral = get_kinetic(za, zb, ra, rb, [1, 0, 0, 0, 0, 0]);
        assert!((integral - 0.2514265587836556).abs() < thresh);

        let integral = get_kinetic(za, zb, ra, rb, [2, 2, 2, 2, 2, 2]);
        assert!((integral - -7.40057384314e-05).abs() < thresh);
    }

    #[test]
    fn test_get_nuclear() {
        let za = 1.8;
        let zb = 2.0;
        let ra = [0.0, 0.0, 0.0];
        let rb = [0.5, 0.8, -0.2];
        let rc = [0.5, 0.8, 0.2];

        // TODO why is this thresh lower now?
        let thresh = 1.0e-15;

        let integral = get_nuclear(za, zb, ra, rb, rc, [0, 0, 0, 0, 0, 0]);
        assert!((integral - -0.49742209545104593).abs() < thresh);

        let integral = get_nuclear(za, zb, ra, rb, rc, [1, 0, 0, 0, 0, 0]);
        assert!((integral - -0.15987439458254471).abs() < thresh);

        let integral = get_nuclear(za, zb, ra, rb, rc, [2, 2, 2, 0, 0, 0]);
        assert!((integral - -0.003801373531942607).abs() < thresh);

        let integral = get_nuclear(za, zb, ra, rb, rc, [1, 1, 1, 1, 1, 1]);
        assert!((integral - 8.8415484347060993e-5).abs() < thresh);
    }

    #[test]
    fn test_get_moment() {
        let za = 1.8;
        let zb = 2.0;
        let ra = [0.0, 0.0, 0.0];
        let rb = [0.5, 0.8, -0.2];
        let rc = [0.0, 0.0, 0.0];

        let thresh = 1.0e-16;

        let integral = get_moment(za, zb, ra, rb, rc, [0, 0, 2, 0, 0, 0], [0, 0, 1]);
        assert!((integral - -0.01330515491323708).abs() < thresh);
    }
}
