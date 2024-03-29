use std::convert::TryInto;
use std::f64::consts::PI;

extern crate rgsl;

use arrayvec::ArrayVec;

use boys::micb25::boys;

fn remove_item<T: PartialEq>(v: &mut Vec<T>, item: &T) {
    v.iter().position(|x| x == item).map(|x| v.remove(x));
}

#[derive(Clone, PartialEq)]
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

#[derive(Clone)]
struct X2 {
    scale: f64,
    prefactors: Vec<usize>,
    q: [i8; 6],
    kind: X2kind,
    operator: [i8; 3],
    d: u8,
    order: u8,
}

#[derive(Clone)]
struct X4 {
    scale: f64,
    prefactors: Vec<usize>,
    q: [i8; 12],
    order: u8,
}

fn apply_os4(x: X4) -> Vec<X4> {
    assert!(x.q.iter().all(|&i| i >= 0));

    let q_sum: i8 = x.q.iter().sum();
    if q_sum == 0 {
        return vec![x];
    }

    let fun = find_fun_to_lower(&x.q.to_vec(), 4);
    if fun == Err(false) {
        return vec![x];
    }

    let fun = fun.unwrap();
    // let possible_components_to_lower = [&orders[fun*3], &orders[fun*3 + 1], &orders[fun*3 + 2]];
    let mut possible_components_to_lower = [0, 0, 0];
    let fun_orders = &x.q[fun * 3..fun * 3 + 3];
    possible_components_to_lower.clone_from_slice(fun_orders);
    let component = find_component_to_lower(&possible_components_to_lower);

    let (i1, i2): (u8, u8) = match component {
        Ok(0) => ([0, 1, 2, 3][fun], [4, 4, 5, 5][fun]),
        Ok(1) => ([6, 7, 8, 9][fun], [10, 10, 11, 11][fun]),
        Ok(2) => ([12, 13, 14, 15][fun], [16, 16, 17, 17][fun]),
        _ => unreachable!(),
    };

    let mut bra = ArrayVec::<usize, 2>::new();
    bra.push(0);
    bra.push(1);
    let mut ket = ArrayVec::<usize, 2>::new();
    ket.push(2);
    ket.push(3);

    let pre = if bra.contains(&fun) {
        [i1, i2, 18, 20, 18, 20, 22, 22]
    } else {
        [i1, i2, 19, 21, 19, 21, 22, 22]
    };

    let a = fun;
    let (b, c, d) = if bra.contains(&a) {
        // remove_item(&mut bra, &a);
        bra.iter().position(|&x| x == a).map(|x| bra.remove(x));
        (bra.pop().unwrap(), ket.pop().unwrap(), ket.pop().unwrap())
    } else {
        // remove_item(&mut ket, &a);
        ket.iter().position(|&x| x == a).map(|x| ket.remove(x));
        (ket.pop().unwrap(), bra.pop().unwrap(), bra.pop().unwrap())
    };

    let component = component.unwrap();
    let num_terms = 8;
    let mut x_copy = Vec::new();
    for term in 0..num_terms {
        x_copy.push(x.clone());
        x_copy[term].q[fun * 3 + component] -= 1;
    }

    let terms = [1, 3, 5, 6, 7];
    for term in terms.iter() {
        x_copy[*term].order += 1;
    }

    x_copy[2].q[a * 3 + component] -= 1;
    x_copy[3].q[a * 3 + component] -= 1;
    x_copy[4].q[b * 3 + component] -= 1;
    x_copy[5].q[b * 3 + component] -= 1;
    x_copy[6].q[c * 3 + component] -= 1;
    x_copy[7].q[d * 3 + component] -= 1;

    let n = [
        1,
        1,
        x.q[a * 3 + component] - 1,
        x.q[a * 3 + component] - 1,
        x.q[b * 3 + component],
        x.q[b * 3 + component],
        x.q[c * 3 + component],
        x.q[d * 3 + component],
    ];

    let mut x_list = Vec::new();
    for term in 0..num_terms {
        if n[term] > 0 {
            if x_copy[term].q.iter().all(|&order| order >= 0) {
                if n[term] > 1 {
                    x_copy[term].scale *= n[term] as f64;
                }
                x_copy[term].prefactors.push(pre[term].try_into().unwrap());
                x_list.push(x_copy[term].clone());
            }
        }
    }

    x_list
        .iter()
        .flat_map(|x4| {
            if x4.q.iter().all(|&order| order == 0) {
                vec![x4.clone()]
            } else {
                apply_os4(x4.clone())
            }
        })
        .collect()
}

fn get_k(z1: f64, z2: f64, r1: &[f64; 3], r2: &[f64; 3]) -> f64 {
    let r12 = get_r12_squared(r1, r2);
    let f0 = z1 + z2;
    let f2 = if r12 > 0.0 {
        (-z1 * z2 * r12 / f0).exp()
    } else {
        1.0
    };
    (2.0_f64).sqrt() * f2 * PI.powf(5.0 / 4.0) / f0
}

fn get_aux(
    za: f64,
    zb: f64,
    zc: f64,
    zd: f64,
    ra: &[f64; 3],
    rb: &[f64; 3],
    rc: &[f64; 3],
    rd: &[f64; 3],
) -> f64 {
    let k1 = get_k(za, zb, ra, rb);
    let k2 = get_k(zc, zd, rc, rd);
    k1 * k2 / (za + zb + zc + zd).sqrt()
}

pub fn get_coulomb(
    za: f64,
    zb: f64,
    zc: f64,
    zd: f64,
    ra: &[f64; 3],
    rb: &[f64; 3],
    rc: &[f64; 3],
    rd: &[f64; 3],
    c: &[usize; 12],
) -> f64 {
    let rp = get_bi_center(za, zb, ra, rb);
    let rq = get_bi_center(zc, zd, rc, rd);
    let rw = get_bi_center(za + zb, zc + zd, &rp, &rq);
    let z = za + zb;
    let n = zc + zd;
    let rho = z * n / (z + n);
    let t = rho * get_r12_squared(&rp, &rq);
    let s = get_aux(za, zb, zc, zd, ra, rb, rc, rd);

    let prefac = [
        rp[0] - ra[0],
        rp[0] - rb[0],
        rq[0] - rc[0],
        rq[0] - rd[0],
        rw[0] - rp[0],
        rw[0] - rq[0],
        rp[1] - ra[1],
        rp[1] - rb[1],
        rq[1] - rc[1],
        rq[1] - rd[1],
        rw[1] - rp[1],
        rw[1] - rq[1],
        rp[2] - ra[2],
        rp[2] - rb[2],
        rq[2] - rc[2],
        rq[2] - rd[2],
        rw[2] - rp[2],
        rw[2] - rq[2],
        0.5 / z,
        0.5 / n,
        -0.5 * rho / (z * z),
        -0.5 * rho / (n * n),
        0.5 / (z + n),
    ];

    let q: [i8; 12] = [
        c[0] as i8,
        c[1] as i8,
        c[2] as i8,
        c[3] as i8,
        c[4] as i8,
        c[5] as i8,
        c[6] as i8,
        c[7] as i8,
        c[8] as i8,
        c[9] as i8,
        c[10] as i8,
        c[11] as i8,
    ];

    let fun = X4 {
        scale: 1.0,
        prefactors: vec![],
        q: q,
        order: 0,
    };
    let expansion = apply_os4(fun);
    let mut integral = 0.0;
    let c_sum: usize = c.iter().sum();
    for i in 0..(c_sum + 1) {
        let b = boys(i as u64, t) * s;
        for f in expansion.iter() {
            if f.order as usize == i {
                let mut g = 1.0;
                for k in f.prefactors.iter() {
                    g *= prefac[*k];
                }
                integral += f.scale * b * g;
            }
        }
    }
    integral
}

fn get_bi_center(z1: f64, z2: f64, r1: &[f64; 3], r2: &[f64; 3]) -> [f64; 3] {
    let z = z1 + z2;
    let rx = (z1 * r1[0] + z2 * r2[0]) / z;
    let ry = (z1 * r1[1] + z2 * r2[1]) / z;
    let rz = (z1 * r1[2] + z2 * r2[2]) / z;
    [rx, ry, rz]
}

fn get_r12_squared(r1: &[f64; 3], r2: &[f64; 3]) -> f64 {
    (r1[0] - r2[0]).powi(2) + (r1[1] - r2[1]).powi(2) + (r1[2] - r2[2]).powi(2)
}

fn find_fun_to_lower(q: &Vec<i8>, n: usize) -> Result<usize, bool> {
    // Determine the total angular momentum on each center.
    let mut l = ArrayVec::<i8, 12>::new();
    for i in 0..n {
        l.push(q[i * 3] + q[i * 3 + 1] + q[i * 3 + 2])
    }

    // find function to lower
    // start with lowest angular momentum above s
    let mut fun: isize = -1;
    let mut kmax = 1 + *l.iter().max().expect("There are no centers?");

    for i in 0..n {
        let k = l[i];
        // If we're larger than a s-function...
        if k > 0 {
            if k < kmax {
                kmax = k;
                fun = i as isize;
            }
        }
    }

    if fun > -1 {
        return Ok(fun as usize);
    }
    return Err(false);
}

fn find_component_to_lower(fun: &[i8; 3]) -> Result<usize, bool> {
    for (i, c) in fun.iter().enumerate() {
        if *c > 0 {
            return Ok(i);
        }
    }
    return Err(false);
}

impl X2 {
    fn orders(&self) -> [i8; 9] {
        [
            self.q[0],
            self.q[1],
            self.q[2],
            self.q[3],
            self.q[4],
            self.q[5],
            self.operator[0],
            self.operator[1],
            self.operator[2],
        ]
    }
}

fn apply_os2(mut x: X2, kind: X2kind) -> Vec<X2> {
    let orders = x.orders();
    assert!(orders.iter().all(|&i| i >= 0));

    // base case
    let order_sum: i8 = orders.iter().sum();
    if order_sum == 0 {
        x.kind = kind;
        return vec![x];
    }

    // Determine which basis function and component to lower.
    // The component is one of (x, y, z).
    let mut fun = find_fun_to_lower(&orders.to_vec(), 3);
    // Make sure to not choose the operator vrr until q is exhausted.
    let q_sum: i8 = x.q.iter().sum();
    if fun == Ok(2) && q_sum > 0 {
        fun = find_fun_to_lower(&x.q.to_vec(), 2)
    }
    if fun == Err(false) {
        x.kind = kind;
        return vec![x];
    }
    let fun = fun.unwrap();
    // let possible_components_to_lower = [&orders[fun*3], &orders[fun*3 + 1], &orders[fun*3 + 2]];
    let mut possible_components_to_lower = [0, 0, 0];
    let fun_orders = &orders[fun * 3..fun * 3 + 3];
    possible_components_to_lower.clone_from_slice(fun_orders);
    let component = find_component_to_lower(&possible_components_to_lower);

    // Determine the index of q to descend on.
    // where q = [xa, xb, ya, yb, za, zb].
    // Note: I thought the order of q == [xa, ya, za, xb, yb, zb]?
    let i1: u8 = match component {
        Ok(0) => [0, 1, 6][fun],
        Ok(1) => [2, 3, 7][fun],
        Ok(2) => [4, 5, 8][fun],
        _ => unreachable!(),
    };

    let mut pre = ArrayVec::<u8, 12>::new();
    pre.push(i1);

    if kind == X2kind::S {
        // The vrr for overlap integrals consists of three "terms".
        pre.push(6);
        pre.push(7);
    }

    if kind == X2kind::T {
        pre.push(6);
        pre.push(7);
        pre.push(8);
        pre.push([9, 10][fun]);
    }

    if kind == X2kind::V {
        let i2: u8 = match component {
            Ok(0) => 6,
            Ok(1) => 7,
            Ok(2) => 8,
            _ => unreachable!(),
        };
        pre.push(i2);
        pre.push(9);
        pre.push(10);
        pre.push(9);
        pre.push(10);
    }

    if kind == X2kind::M {
        pre.push(9);
        pre.push(10);
        pre.push(11);
    }

    if kind == X2kind::L {
        pre.push(6);
        pre.push(7);
        pre.push(8 + x.d);
        pre.push(11);
        pre.push(12);
        pre.push(13);
    }

    // Determine which of the basis functions is ("a", "b").
    let (a, b) = match fun {
        2 => (2, 2),
        _ => {
            let mut l = ArrayVec::<usize, 2>::new();
            l.push(0);
            l.push(1);
            // remove_item(&mut l, &fun);
            l.iter().position(|&x| x == fun).map(|x| l.remove(x));
            (fun, l[0])
        }
    };

    // These are the number of integrals that appear in the main
    // recursion equations for each kind.
    let num_terms = match kind {
        X2kind::S => 3,
        X2kind::T => 5,
        X2kind::V => 6,
        X2kind::M => match fun {
            2 => 2,
            _ => 4,
        },
        X2kind::L => 7,
        X2kind::E => 0,
        X2kind::J => 0,
        X2kind::FC => 0,
    };

    // Make copies of the current integral to manipulate later,
    // one fo reach term in the recursion expression.
    let mut x_copy = Vec::new();
    for _ in 0..num_terms {
        x_copy.push(x.clone());
    }

    // These are the terms in [A19] with (m + 1).
    if kind == X2kind::V {
        let terms: [usize; 3] = [1, 3, 5];
        for term in terms.iter() {
            x_copy[*term].order += 1;
        }
    }

    // Look at the last line of [A12].
    if kind == X2kind::T {
        x_copy[3].kind = X2kind::S;
        x_copy[4].kind = X2kind::S;
    }

    // Look at the last two lines of [A31].
    if kind == X2kind::L {
        x_copy[3].kind = X2kind::S;
        x_copy[4].kind = X2kind::S;
        x_copy[5].kind = X2kind::S;
        x_copy[6].kind = X2kind::S;
    }

    let component = component.unwrap();
    // 1. Lower the target component for all three terms.
    // 2. Lower again on center "a".
    // 3. Lower again on center "b".
    if kind == X2kind::S {
        x_copy[0].q[fun * 3 + component] -= 1;
        x_copy[1].q[fun * 3 + component] -= 1;
        x_copy[2].q[fun * 3 + component] -= 1;
        x_copy[1].q[a * 3 + component] -= 1;
        x_copy[2].q[b * 3 + component] -= 1;
    }

    if kind == X2kind::T {
        x_copy[0].q[fun * 3 + component] -= 1;
        x_copy[1].q[fun * 3 + component] -= 1;
        x_copy[2].q[fun * 3 + component] -= 1;
        // term 4 (x_copy[3]) has the same components but becomes an
        // overlap integral
        x_copy[4].q[fun * 3 + component] -= 2;
        x_copy[1].q[a * 3 + component] -= 1;
        x_copy[2].q[b * 3 + component] -= 1;
    }

    if kind == X2kind::V {
        x_copy[0].q[fun * 3 + component] -= 1;
        x_copy[2].q[fun * 3 + component] -= 1;
        x_copy[4].q[fun * 3 + component] -= 1;
        x_copy[2].q[a * 3 + component] -= 1;
        x_copy[4].q[b * 3 + component] -= 1;
        x_copy[1].q = x_copy[0].q;
        x_copy[3].q = x_copy[2].q;
        x_copy[5].q = x_copy[4].q;
    }

    if kind == X2kind::M {
        // For the case of the moment operator over s functions.
        if fun == 2 {
            x_copy[0].operator[component] -= 1;
            x_copy[1].operator[component] -= 2;
        } else {
            x_copy[0].q[fun * 3 + component] -= 1;
            x_copy[1].q[fun * 3 + component] -= 1;
            x_copy[2].q[fun * 3 + component] -= 1;
            x_copy[3].q[fun * 3 + component] -= 1;
            x_copy[1].q[a * 3 + component] -= 1;
            x_copy[2].q[b * 3 + component] -= 1;
            x_copy[3].operator[component] -= 1;
        }
    }

    if kind == X2kind::L {
        x_copy[0].q[fun * 3 + component] -= 1;
        x_copy[1].q[fun * 3 + component] -= 1;
        x_copy[2].q[fun * 3 + component] -= 1;
        x_copy[3].q[fun * 3 + component] -= 1;
        x_copy[4].q[fun * 3 + component] -= 1;
        x_copy[5].q[fun * 3 + component] -= 1;
        x_copy[6].q[fun * 3 + component] -= 1;
        x_copy[1].q[a * 3 + component] -= 1;
        x_copy[2].q[b * 3 + component] -= 1;
        x_copy[4].q[b * 3 + 0] -= 1;
        x_copy[5].q[b * 3 + 1] -= 1;
        x_copy[6].q[b * 3 + 2] -= 1;
    }

    // Now that the descending part of the vrr has been performed, keep track
    // of which new terms are going to be zero/non-zero.
    let mut n: Vec<i8> = Vec::new();
    // The first term is always going to be non-zero, otherwise we wouldn't
    // even be in the vrr routine.
    n.push(1);

    if kind == X2kind::S {
        n.push(x.q[a * 3 + component] - 1);
        n.push(x.q[b * 3 + component]);
    }

    if kind == X2kind::T {
        n.push(x.q[a * 3 + component] - 1);
        n.push(x.q[b * 3 + component]);
        n.push(1);
        n.push(x.q[a * 3 + component] - 1);
    }

    if kind == X2kind::V {
        n.push(1);
        n.push(x.q[a * 3 + component] - 1);
        n.push(x.q[a * 3 + component] - 1);
        n.push(x.q[b * 3 + component]);
        n.push(x.q[b * 3 + component]);
    }

    if kind == X2kind::M {
        if fun == 2 {
            n.push(x.operator[component] - 1);
        } else {
            n.push(x.q[a * 3 + component] - 1);
            n.push(x.q[b * 3 + component]);
            n.push(x.operator[component]);
        }
    }

    if kind == X2kind::L {
        n.push(x.q[a * 3 + component] - 1);
        n.push(x.q[b * 3 + component]);
        n.push(
            (0..3)
                .map(|d| if d == component { 1 } else { 0 })
                .collect::<Vec<_>>()[x.d as usize],
        );
        n.push(
            (0..3)
                .map(|d| if d == component && d == 0 { 1 } else { 0 })
                .collect::<Vec<_>>()[x.d as usize],
        );
        n.push(
            (0..3)
                .map(|d| if d == component && d == 1 { 1 } else { 0 })
                .collect::<Vec<_>>()[x.d as usize],
        );
        n.push(
            (0..3)
                .map(|d| if d == component && d == 2 { 1 } else { 0 })
                .collect::<Vec<_>>()[x.d as usize],
        );
    }

    // TODO remainder of X2kind

    // Generate a list of all non-zero terms for an expression.
    let mut x_list: Vec<X2> = Vec::new();
    for term in 0..num_terms {
        if n[term] > 0 {
            if x_copy[term].orders().iter().all(|&order| order >= 0) {
                if n[term] > 1 {
                    x_copy[term].scale *= n[term] as f64;
                }
                x_copy[term].prefactors.push(pre[term].try_into().unwrap());
                x_list.push(x_copy[term].clone());
            }
        }
    }

    // If we've hit the base case where all the components are zero,
    // terminate, otherwise recurse through each term.
    x_list
        .iter()
        .flat_map(|x2| {
            if x2.orders().iter().all(|&order| order == 0) {
                vec![x2.clone()]
            } else {
                apply_os2(x2.clone(), x2.kind.clone())
            }
        })
        .collect()
}

pub fn get_overlap(za: f64, zb: f64, ra: &[f64; 3], rb: &[f64; 3], c: &[usize; 6]) -> f64 {
    let z = za + zb;
    let e = za * zb / (za + zb);
    let rp = get_bi_center(za, zb, ra, rb);
    let ab = get_r12_squared(ra, rb);
    let aux = (-e * ab).exp() * (PI / z).powf(1.5);

    let prefac = [
        rp[0] - ra[0],
        rp[0] - rb[0],
        rp[1] - ra[1],
        rp[1] - rb[1],
        rp[2] - ra[2],
        rp[2] - rb[2],
        0.5 / z,
        0.5 / z,
    ];

    let q: [i8; 6] = [
        c[0] as i8, c[1] as i8, c[2] as i8, c[3] as i8, c[4] as i8, c[5] as i8,
    ];

    let fun = X2 {
        scale: 1.0,
        prefactors: vec![],
        q: q,
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
    integral
}

pub fn get_kinetic(za: f64, zb: f64, ra: &[f64; 3], rb: &[f64; 3], c: &[usize; 6]) -> f64 {
    let z = za + zb;
    let e = za * zb / (za + zb);
    let rp = get_bi_center(za, zb, ra, rb);
    let ab = get_r12_squared(ra, rb);
    let aux = (-e * ab).exp() * (PI / z).powf(1.5);

    let prefac = [
        rp[0] - ra[0],
        rp[0] - rb[0],
        rp[1] - ra[1],
        rp[1] - rb[1],
        rp[2] - ra[2],
        rp[2] - rb[2],
        0.5 / z,
        0.5 / z,
        2.0 * e,
        -e / za,
        -e / zb,
    ];

    let q: [i8; 6] = [
        c[0] as i8, c[1] as i8, c[2] as i8, c[3] as i8, c[4] as i8, c[5] as i8,
    ];

    let fun = X2 {
        scale: 1.0,
        prefactors: vec![],
        q: q,
        kind: X2kind::T,
        operator: [0, 0, 0],
        d: 0,
        order: 0,
    };
    let expansion = apply_os2(fun, X2kind::T);
    let mut integral = 0.0;
    for f in expansion.iter() {
        let mut g = match f.kind {
            X2kind::T => e * (3.0 - 2.0 * e * ab),
            X2kind::S => 1.0,
            _ => unreachable!(),
        };
        for k in f.prefactors.iter() {
            g *= prefac[*k];
        }
        integral += f.scale * aux * g;
    }
    integral
}

pub fn get_nuclear(
    za: f64,
    zb: f64,
    ra: &[f64; 3],
    rb: &[f64; 3],
    rc: &[f64; 3],
    c: &[usize; 6],
) -> f64 {
    let z = za + zb;
    let rp = get_bi_center(za, zb, ra, rb);
    let pc = get_r12_squared(&rp, rc);
    let u = z * pc;
    let aux = -2.0 * (z / PI).powf(0.5) * get_overlap(za, zb, ra, rb, &[0, 0, 0, 0, 0, 0]);

    let prefac = [
        rp[0] - ra[0],
        rp[0] - rb[0],
        rp[1] - ra[1],
        rp[1] - rb[1],
        rp[2] - ra[2],
        rp[2] - rb[2],
        -rp[0] + rc[0],
        -rp[1] + rc[1],
        -rp[2] + rc[2],
        0.5 / z,
        -0.5 / z,
    ];

    let q: [i8; 6] = [
        c[0] as i8, c[1] as i8, c[2] as i8, c[3] as i8, c[4] as i8, c[5] as i8,
    ];

    let fun = X2 {
        scale: 1.0,
        prefactors: vec![],
        q: q,
        kind: X2kind::V,
        operator: [0, 0, 0],
        d: 0,
        order: 0,
    };
    let expansion = apply_os2(fun, X2kind::V);
    let mut integral = 0.0;
    for f in expansion.iter() {
        let mut g = 1.0;
        for k in f.prefactors.iter() {
            g *= prefac[*k];
        }
        integral += f.scale * aux * g * boys(f.order.into(), u);
    }
    integral
}

fn get_moment(
    za: f64,
    zb: f64,
    ra: &[f64; 3],
    rb: &[f64; 3],
    rc: &[f64; 3],
    c: &[u8; 6],
    order: &[u8; 3],
) -> f64 {
    let z = za + zb;
    let rp = get_bi_center(za, zb, ra, rb);
    let aux = get_overlap(za, zb, ra, rb, &[0, 0, 0, 0, 0, 0]);

    let prefac = [
        rp[0] - ra[0],
        rp[0] - rb[0],
        rp[1] - ra[1],
        rp[1] - rb[1],
        rp[2] - ra[2],
        rp[2] - rb[2],
        rp[0] - rc[0],
        rp[1] - rc[1],
        rp[2] - rc[2],
        0.5 / z,
        0.5 / z,
        0.5 / z,
    ];

    let q: [i8; 6] = [
        c[0] as i8, c[1] as i8, c[2] as i8, c[3] as i8, c[4] as i8, c[5] as i8,
    ];
    let operator: [i8; 3] = [order[0] as i8, order[1] as i8, order[2] as i8];

    let fun = X2 {
        scale: 1.0,
        prefactors: vec![],
        q: q,
        kind: X2kind::M,
        operator: operator,
        d: 0,
        order: 0,
    };

    let expansion = apply_os2(fun, X2kind::M);
    let mut integral = 0.0;
    for f in expansion.iter() {
        let mut g = 1.0;
        for k in f.prefactors.iter() {
            g *= prefac[*k];
        }
        integral += f.scale * aux * g;
    }
    integral
}

#[cfg(test)]
mod tests {
    use super::find_component_to_lower;
    use super::find_fun_to_lower;
    use super::get_coulomb;
    use super::get_kinetic;
    use super::get_moment;
    use super::get_nuclear;
    use super::get_overlap;

    #[test]
    fn test_find_fun_to_lower() {
        assert_eq!(find_fun_to_lower(&[1, 0, 0, 0, 0, 0].to_vec(), 2), Ok(0));
        assert_eq!(find_fun_to_lower(&[0, 1, 0, 0, 0, 0].to_vec(), 2), Ok(0));
        assert_eq!(find_fun_to_lower(&[0, 0, 1, 0, 0, 0].to_vec(), 2), Ok(0));
        assert_eq!(find_fun_to_lower(&[0, 0, 0, 1, 0, 0].to_vec(), 2), Ok(1));
        assert_eq!(find_fun_to_lower(&[0, 0, 0, 0, 1, 0].to_vec(), 2), Ok(1));
        assert_eq!(find_fun_to_lower(&[0, 0, 0, 0, 0, 1].to_vec(), 2), Ok(1));
        assert_eq!(
            find_fun_to_lower(&[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0].to_vec(), 4),
            Ok(1)
        );
        assert_eq!(
            find_fun_to_lower(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1].to_vec(), 4),
            Ok(3)
        );
        assert_eq!(
            find_fun_to_lower(&[1, 0, 0, 0, 0, 0, 0, 0, 1].to_vec(), 3),
            Ok(0)
        );
        assert_eq!(
            find_fun_to_lower(&[0, 0, 0, 0, 0, 0, 0, 0, 1].to_vec(), 3),
            Ok(2)
        );
        assert_eq!(
            find_fun_to_lower(&[0, 0, 0, 0, 1, 0, 0, 0, 1].to_vec(), 3),
            Ok(1)
        );
        assert_eq!(
            find_fun_to_lower(&[0, 0, 0, 0, 2, 0, 0, 0, 1].to_vec(), 3),
            Ok(2)
        );
    }

    #[test]
    fn test_find_component_to_lower() {
        assert_eq!(find_component_to_lower(&[0, 0, 1]), Ok(2));
        assert_eq!(find_component_to_lower(&[0, 1, 1]), Ok(1));
        assert_eq!(find_component_to_lower(&[1, 0, 1]), Ok(0));
        assert_eq!(find_component_to_lower(&[0, 0, 0]), Err(false));
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

        let thresh = 1.0e-15;

        let reference = 0.08608517834596989;
        let integral = get_coulomb(
            za,
            zb,
            zc,
            zd,
            &ra,
            &rb,
            &rc,
            &rd,
            &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        );
        assert!((integral - reference).abs() < thresh);

        let reference = -0.046888060557776973;
        let integral = get_coulomb(
            za,
            zb,
            zc,
            zd,
            &ra,
            &rb,
            &rc,
            &rd,
            &[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        );
        assert!((integral - reference).abs() < thresh);

        let reference = 1.71817807954e-05;
        let integral = get_coulomb(
            za,
            zb,
            zc,
            zd,
            &ra,
            &rb,
            &rc,
            &rd,
            &[2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0],
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

        let integral = get_overlap(za, zb, &ra, &rb, &[0, 0, 0, 0, 0, 0]);
        assert!((integral - 0.20373275913014607).abs() < thresh);

        let integral = get_overlap(za, zb, &ra, &rb, &[1, 0, 0, 0, 0, 0]);
        assert!((integral - 0.062005622343957505).abs() < thresh);

        let integral = get_overlap(za, zb, &ra, &rb, &[1, 1, 0, 1, 1, 0]);
        assert!((integral - -0.00043801221837779696).abs() < thresh);

        let integral = get_overlap(za, zb, &ra, &rb, &[2, 1, 0, 1, 1, 0]);
        assert!((integral - -0.0002385994651113168).abs() < thresh);
    }

    #[test]
    fn test_get_kinetic() {
        let za = 1.8;
        let zb = 2.0;
        let ra = [0.0, 0.0, 0.0];
        let rb = [0.5, 0.8, -0.2];

        let thresh = 1.0e-16;

        let integral = get_kinetic(za, zb, &ra, &rb, &[0, 0, 0, 0, 0, 0]);
        assert!((integral - 0.3652714583525358).abs() < thresh);

        let integral = get_kinetic(za, zb, &ra, &rb, &[1, 0, 0, 0, 0, 0]);
        assert!((integral - 0.2514265587836556).abs() < thresh);

        let integral = get_kinetic(za, zb, &ra, &rb, &[2, 2, 2, 2, 2, 2]);
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

        let integral = get_nuclear(za, zb, &ra, &rb, &rc, &[0, 0, 0, 0, 0, 0]);
        assert!((integral - -0.49742209545104593).abs() < thresh);

        let integral = get_nuclear(za, zb, &ra, &rb, &rc, &[1, 0, 0, 0, 0, 0]);
        assert!((integral - -0.15987439458254471).abs() < thresh);

        let integral = get_nuclear(za, zb, &ra, &rb, &rc, &[2, 2, 2, 0, 0, 0]);
        assert!((integral - -0.003801373531942607).abs() < thresh);

        let integral = get_nuclear(za, zb, &ra, &rb, &rc, &[1, 1, 1, 1, 1, 1]);
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

        let integral = get_moment(za, zb, &ra, &rb, &rc, &[0, 0, 2, 0, 0, 0], &[0, 0, 1]);
        assert!((integral - -0.01330515491323708).abs() < thresh);
    }
}
