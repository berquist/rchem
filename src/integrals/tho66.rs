#![allow(clippy::too_many_arguments)]
#![allow(clippy::many_single_char_names)]
use factorial::SignedFactorial;
use num_integer;
use std::convert::TryInto;
use std::f64;

fn fact_ratio2(a: i64, b: i64) -> i64 {
    a.factorial() / b.factorial() / (a - 2 * b).factorial()
}

fn bfunc(i: i64, r: i64, g: f64) -> f64 {
    (fact_ratio2(i, r) as f64) * (4.0 * g).powf((r - i) as f64)
}

fn binomial_prefactor(s: i64, ia: i64, ib: i64, xpa: f64, xpb: f64) -> f64 {
    let mut result = 0.0;
    for t in 0..=s {
        if ((s - ia) <= t) && (t <= ib) {
            result += (num_integer::binomial(ia, s - t) * num_integer::binomial(ib, t)) as f64
                * xpa.powi((ia - s + t) as i32)
                * xpb.powi((ib - t) as i32);
        }
    }
    result
}

const SMALL: f64 = 1.0e-8;

fn fgamma(m: f64, x: f64) -> f64 {
    let actual_x = if x.abs() < SMALL { SMALL } else { x };
    0.5 * actual_x.powf(-m - 0.5)
        * rgsl::gamma_beta::incomplete_gamma::gamma_inc_P(m + 0.5, actual_x)
}

fn fb(i: i64, l1: i64, l2: i64, px: f64, ax: f64, bx: f64, r: i64, g: f64) -> f64 {
    binomial_prefactor(i, l1, l2, px - ax, px - bx) * bfunc(i, r, g)
}

fn dist2(xa: f64, ya: f64, za: f64, xb: f64, yb: f64, zb: f64) -> f64 {
    (xa - xb).powf(2.0) + (ya - yb).powf(2.0) + (za - zb).powf(2.0)
}

fn product_center_1d(za: f64, xa: f64, zb: f64, xb: f64) -> f64 {
    ((za * xa) + (zb * xb)) / (za + zb)
}

fn b_term(
    i1: i64,
    i2: i64,
    r1: i64,
    r2: i64,
    u: i64,
    l1: i64,
    l2: i64,
    l3: i64,
    l4: i64,
    px: f64,
    ax: f64,
    bx: f64,
    qx: f64,
    cx: f64,
    dx: f64,
    gamma1: f64,
    gamma2: f64,
    delta: f64,
) -> f64 {
    fb(i1, l1, l2, px, ax, bx, r1, gamma1)
        * ((-1_i64).pow(i2.try_into().unwrap()) as f64)
        * fb(i2, l3, l4, qx, cx, dx, r2, gamma2)
        * ((-1_i64).pow(u.try_into().unwrap()) as f64)
        * (fact_ratio2(i1 + i2 - 2 * (r1 + r2), u) as f64)
        * (qx - px).powf((i1 + i2 - 2 * (r1 + r2) - 2 * u) as f64)
        / delta.powi((i1 + i2 - 2 * (r1 + r2) - u) as i32)
}

fn b_array(
    l1: i64,
    l2: i64,
    l3: i64,
    l4: i64,
    p: f64,
    a: f64,
    b: f64,
    q: f64,
    c: f64,
    d: f64,
    g1: f64,
    g2: f64,
    delta: f64,
) -> Vec<f64> {
    let mut res = vec![0.0; (l1 + l2 + l3 + l4 + 1) as usize];
    for i1 in 0..l1 + l2 + 1 {
        for i2 in 0..l3 + l4 + 1 {
            for r1 in 0..i1 / 2 + 1 {
                for r2 in 0..i2 / 2 + 1 {
                    for u in 0..(i1 + i2) / 2 - r1 - r2 + 1 {
                        let i = (i1 + i2 - 2 * (r1 + r2) - u) as usize;
                        res[i] += b_term(
                            i1, i2, r1, r2, u, l1, l2, l3, l4, p, a, b, q, c, d, g1, g2, delta,
                        );
                    }
                }
            }
        }
    }
    res
}

// fn coulomb_repulsion(
//     xa: f64,
//     ya: f64,
//     za: f64,
//     norma: f64,
//     la: i64,
//     ma: i64,
//     na: i64,
//     alphaa: f64,
//     xb: f64,
//     yb: f64,
//     zb: f64,
//     normb: f64,
//     lb: i64,
//     mb: i64,
//     nb: i64,
//     alphab: f64,
//     xc: f64,
//     yc: f64,
//     zc: f64,
//     normc: f64,
//     lc: i64,
//     mc: i64,
//     nc: i64,
//     alphac: f64,
//     xd: f64,
//     yd: f64,
//     zd: f64,
//     normd: f64,
//     ld: i64,
//     md: i64,
//     nd: i64,
//     alphad: f64,
// ) -> f64 {
//     let rab2 = dist2(xa, ya, za, xb, yb, zb);
//     let rcd2 = dist2(xc, yc, zc, xd, yd, zd);
//     let xp = product_center_1d(alphaa, xa, alphab, xb);
//     let yp = product_center_1d(alphaa, ya, alphab, yb);
//     let zp = product_center_1d(alphaa, za, alphab, zb);
//     let xq = product_center_1d(alphac, xc, alphad, xd);
//     let yq = product_center_1d(alphac, yc, alphad, yd);
//     let zq = product_center_1d(alphac, zc, alphad, zd);
//     let rpq2 = dist2(xp, yp, zp, xq, yq, zq);
//     let gamma1 = alphaa + alphab;
//     let gamma2 = alphac + alphad;
//     let delta = (1. / gamma1 + 1. / gamma2) / 4.;
//     let bx = b_array(
//         la, lb, lc, ld, xp, xa, xb, xq, xc, xd, gamma1, gamma2, delta,
//     );
//     let by = b_array(
//         ma, mb, mc, md, yp, ya, yb, yq, yc, yd, gamma1, gamma2, delta,
//     );
//     let bz = b_array(
//         na, nb, nc, nd, zp, za, zb, zq, zc, zd, gamma1, gamma2, delta,
//     );
//     let mut sum = 0.0;
//     let imax = (la + lb + lc + ld + 1) as usize;
//     let jmax = (ma + mb + mc + md + 1) as usize;
//     let kmax = (na + nb + nc + nd + 1) as usize;
//     for i in 0..imax {
//         for j in 0..jmax {
//             for k in 0..kmax {
//                 sum += bx[i] * by[j] * bz[k] * fgamma((i + j + k) as f64, 0.25 * rpq2 / delta);
//             }
//         }
//     }
//     2.0 * f64::consts::PI.powf(2.5)
//         * (-alphaa * alphab * rab2 / gamma1).exp()
//         * (-alphac * alphad * rcd2 / gamma2).exp()
//         * sum
//         * norma
//         * normb
//         * normc
//         * normd
//         / (gamma1 * gamma2 * (gamma1 + gamma2).sqrt())
// }

#[link(name = "pyquante2", kind = "static")]
extern "C" {
    fn coulomb_repulsion(
        xa: f64,
        ya: f64,
        za: f64,
        norma: f64,
        la: i64,
        ma: i64,
        na: i64,
        alphaa: f64,
        xb: f64,
        yb: f64,
        zb: f64,
        normb: f64,
        lb: i64,
        mb: i64,
        nb: i64,
        alphab: f64,
        xc: f64,
        yc: f64,
        zc: f64,
        normc: f64,
        lc: i64,
        mc: i64,
        nc: i64,
        alphac: f64,
        xd: f64,
        yd: f64,
        zd: f64,
        normd: f64,
        ld: i64,
        md: i64,
        nd: i64,
        alphad: f64,
    ) -> f64;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dist2() {
        assert_abs_diff_eq!(dist2(0.5, 0.6, 0.7, 0.8, 0.9, 1.0), 0.27);
    }

    #[test]
    fn test_product_center_1d() {
        assert_eq!(product_center_1d(2.0, 3.0, 4.0, 5.0), 26.0 / 6.0);
    }

    #[test]
    fn test_binomial_prefactor() {
        assert_abs_diff_eq!(binomial_prefactor(1, 1, 1, 0.1, 0.2), 0.3);
        assert_abs_diff_eq!(binomial_prefactor(1, 1, 1, 0.3, 0.4), 0.7);
        assert_abs_diff_eq!(binomial_prefactor(1, 3, 1, 0.1, 0.2), 0.007);
        assert_abs_diff_eq!(binomial_prefactor(2, 3, 1, 0.1, 0.2), 0.09);
    }

    #[test]
    fn test_coulomb_repulsion() {
        let za = 1.1;
        let zb = 1.2;
        let zc = 1.3;
        let zd = 1.4;

        let ra = [1.0, 0.0, 1.0];
        let rb = [0.0, 1.0, 2.0];
        let rc = [0.0, 0.0, 3.0];
        let rd = [0.0, 0.0, 4.0];

        let thresh = 1.0e-12;

        let reference = 0.08608517834596989;
        let integral = unsafe {
            coulomb_repulsion(
                ra[0], ra[1], ra[2], 1.0, 0, 0, 0, za, rb[0], rb[1], rb[2], 1.0, 0, 0, 0, zb,
                rc[0], rc[1], rc[2], 1.0, 0, 0, 0, zc, rd[0], rd[1], rd[2], 1.0, 0, 0, 0, zd,
            )
        };
        println!("{}", integral);
        println!("{}", reference);
        assert!((integral - reference).abs() < thresh);
    }
}
