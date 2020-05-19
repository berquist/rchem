//! Two-electron integrals based on the Obara--Saika algorith, provided by [Simint](http://www.bennyp.org/research/simint/).

#![allow(clippy::unreadable_literal)]
#![allow(dead_code)]
#![allow(non_snake_case)]
include!(concat!(env!("OUT_DIR"), "/bindings_simint.rs"));

impl simint_shell {
    //! Create a new **uninitialized** internal `simint_shell` with
    //! mutable null pointers.
    fn new() -> simint_shell {
        simint_shell {
            am: 0,
            nprim: 0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
            alpha: std::ptr::null_mut(),
            coef: std::ptr::null_mut(),
            memsize: 0,
            ptr: std::ptr::null_mut(),
        }
    }
}

#[derive(Clone, Debug)]
struct SimintShell {
    am: i32,
    nprim: i32,
    x: f64,
    y: f64,
    z: f64,
    alpha: Vec<f64>,
    coef: Vec<f64>,
}

impl SimintShell {
    // fn new() -> SimintShell {
    //     SimintShell {
    //         am: 0,
    //         nprim: 0,
    //         x: 0.0,
    //         y: 0.0,
    //         z: 0.0,
    //         alpha: Vec::new(),
    //         coef: Vec::new(),
    //     }
    // }

    fn from(s: &simint_shell) -> SimintShell {
        if s.am < 0 {
            panic!();
        }
        if s.nprim < 0 {
            panic!();
        }
        SimintShell {
            am: s.am,
            nprim: s.nprim,
            x: s.x,
            y: s.y,
            z: s.z,
            alpha: unsafe { std::slice::from_raw_parts(s.alpha, s.nprim as usize).to_vec() },
            coef: unsafe { std::slice::from_raw_parts(s.coef, s.nprim as usize).to_vec() },
        }
        // TODO Since I'm not doing anything with `s` outside this
        // function, I'm probably leaking memory.
    }

    fn to(self) -> simint_shell {
        unsafe {
            simint_init();
        };
        // Strategy: initialize an empty simint_shell, allocate its
        // memory, then copy over this SimintShell's contents.
        let mut internal_shell = simint_shell::new();
        let p_internal_shell = &mut internal_shell as *mut simint_shell;
        unsafe {
            simint_allocate_shell(self.nprim, p_internal_shell);
        };
        if internal_shell.alpha.is_null() {
            panic!("internal_shell.alpha is null");
        }
        if internal_shell.coef.is_null() {
            panic!("internal_shell.coef is null");
        }
        internal_shell.am = self.am;
        internal_shell.nprim = self.nprim;
        internal_shell.x = self.x;
        internal_shell.y = self.y;
        internal_shell.z = self.z;
        unsafe {
            self.alpha
                .as_ptr()
                .copy_to(internal_shell.alpha, self.alpha.len());
            self.coef
                .as_ptr()
                .copy_to(internal_shell.coef, self.coef.len());
        }
        internal_shell
    }
}

impl PartialEq for SimintShell {
    fn eq(&self, other: &Self) -> bool {
        self.am == other.am
            && self.nprim == other.nprim
            && self.x == other.x
            && self.y == other.y
            && self.z == other.z
            && self.alpha == other.alpha
            && self.coef == other.coef
    }
}

fn normalize_shells(shells: Vec<SimintShell>) -> Vec<SimintShell> {
    let mut internal_shells: Vec<_> = shells.iter().map(|shell| shell.clone().to()).collect();
    unsafe {
        simint_normalize_shells(internal_shells.len() as i32, internal_shells.as_mut_ptr());
    };
    internal_shells
        .iter()
        .map(|shell| SimintShell::from(&shell.clone()))
        .collect()
}

impl simint_multi_shellpair {
    //! Create a new **uninitialized** internal `simint_multi_shellpair` with
    //! mutable null pointers.
    fn new() -> simint_multi_shellpair {
        simint_multi_shellpair {
            am1: 0,
            am2: 0,
            nprim: 0,
            nshell12: 0,
            nshell12_clip: 0,
            nprim12: std::ptr::null_mut(),
            AB_x: std::ptr::null_mut(),
            AB_y: std::ptr::null_mut(),
            AB_z: std::ptr::null_mut(),
            x: std::ptr::null_mut(),
            y: std::ptr::null_mut(),
            z: std::ptr::null_mut(),
            PA_x: std::ptr::null_mut(),
            PA_y: std::ptr::null_mut(),
            PA_z: std::ptr::null_mut(),
            PB_x: std::ptr::null_mut(),
            PB_y: std::ptr::null_mut(),
            PB_z: std::ptr::null_mut(),
            alpha: std::ptr::null_mut(),
            prefac: std::ptr::null_mut(),
            screen: std::ptr::null_mut(),
            screen_max: 0.0,
            memsize: 0,
            ptr: std::ptr::null_mut(),
        }
    }
}

struct SimintMultiShellpair {
    am1: i32,
    am2: i32,
    nprim: i32,
    nshell12: i32,
    nshell12_clip: i32,
    nprim12: Vec<i32>,
    dist_ab_x: Vec<f64>,
    dist_ab_y: Vec<f64>,
    dist_ab_z: Vec<f64>,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    pa_x: Vec<f64>,
    pa_y: Vec<f64>,
    pa_z: Vec<f64>,
    pb_x: Vec<f64>,
    pb_y: Vec<f64>,
    pb_z: Vec<f64>,
    alpha: Vec<f64>,
    // alpha2: Vec<f64>,
    // beta2: Vec<f64>,
    prefac: Vec<f64>,
    screen_method: SimintScreenMethod,
    screen: Vec<f64>,
    screen_max: f64,
}

#[derive(Clone, Copy)]
enum SimintScreenMethod {
    None,
    Schwarz,
    FastSchwarz,
}

impl SimintScreenMethod {
    fn from(sm: u32) -> SimintScreenMethod {
        match sm {
            SIMINT_SCREEN_NONE => SimintScreenMethod::None,
            SIMINT_SCREEN_SCHWARZ => SimintScreenMethod::Schwarz,
            SIMINT_SCREEN_FASTSCHWARZ => SimintScreenMethod::FastSchwarz,
            _ => panic!(),
        }
    }

    fn to(sm: SimintScreenMethod) -> u32 {
        match sm {
            SimintScreenMethod::None => SIMINT_SCREEN_NONE,
            SimintScreenMethod::Schwarz => SIMINT_SCREEN_SCHWARZ,
            SimintScreenMethod::FastSchwarz => SIMINT_SCREEN_FASTSCHWARZ,
        }
    }
}

impl SimintMultiShellpair {
    fn from(s: &simint_multi_shellpair, screen_method: SimintScreenMethod) -> SimintMultiShellpair {
        // TODO check null pointers
        let nprim = s.nprim as usize;
        let npair = s.nshell12 as usize;
        SimintMultiShellpair {
            am1: s.am1,
            am2: s.am2,
            nprim: s.nprim,
            nshell12: s.nshell12,
            nshell12_clip: s.nshell12_clip,
            nprim12: unsafe { std::slice::from_raw_parts(s.nprim12, npair).to_vec() },
            dist_ab_x: unsafe { std::slice::from_raw_parts(s.AB_x, npair).to_vec() },
            dist_ab_y: unsafe { std::slice::from_raw_parts(s.AB_y, npair).to_vec() },
            dist_ab_z: unsafe { std::slice::from_raw_parts(s.AB_z, npair).to_vec() },
            x: unsafe { std::slice::from_raw_parts(s.x, nprim).to_vec() },
            y: unsafe { std::slice::from_raw_parts(s.y, nprim).to_vec() },
            z: unsafe { std::slice::from_raw_parts(s.z, nprim).to_vec() },
            pa_x: unsafe { std::slice::from_raw_parts(s.PA_x, nprim).to_vec() },
            pa_y: unsafe { std::slice::from_raw_parts(s.PA_y, nprim).to_vec() },
            pa_z: unsafe { std::slice::from_raw_parts(s.PA_z, nprim).to_vec() },
            pb_x: unsafe { std::slice::from_raw_parts(s.PB_x, nprim).to_vec() },
            pb_y: unsafe { std::slice::from_raw_parts(s.PB_y, nprim).to_vec() },
            pb_z: unsafe { std::slice::from_raw_parts(s.PB_z, nprim).to_vec() },
            alpha: unsafe { std::slice::from_raw_parts(s.alpha, nprim).to_vec() },
            // TODO We make the choice to hard-code the assumption that SIMINT
            // cannot handle derivatives.
            // alpha2: Vec::new(),
            // beta2: Vec::new(),
            prefac: unsafe { std::slice::from_raw_parts(s.prefac, nprim).to_vec() },
            screen_method: screen_method,
            screen: match screen_method {
                SimintScreenMethod::None => Vec::new(),
                _ => unsafe { std::slice::from_raw_parts(s.screen, nprim).to_vec() },
            },
            screen_max: s.screen_max,
        }
    }

    // fn to(self) -> simint_multi_shellpair {
    //     unsafe {
    //         simint_init();
    //     };
    //     // Strategy: initialize an empty simint_multi_shellpair, allocate its
    //     // memory, then copy over this SimintMultiShellpair's contents.
    //     let mut internal_multi_shellpair = simint_multi_shellpair::new();
    //     let p_internal_multi_shellpair =
    //         &mut internal_multi_shellpair as *mut simint_multi_shellpair;
    //     unsafe {
    //         simint_initialize_multi_shellpair(p_internal_multi_shellpair);
    //     };
    //     internal_multi_shellpair.am1 = self.am1;
    //     internal_multi_shellpair.am2 = self.am2;
    //     internal_multi_shellpair
    // }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    /// A reimplementation of [simint/examples/example1.c](https://github.com/simint-chem/simint-generator/blob/c589bd70e53bbdde1753df093bea6db2eb63c971/skel/examples/example1.c)
    fn example1() {
        let s_shells = vec![
            SimintShell {
                am: 0,
                nprim: 3,
                x: 0.00000,
                y: -0.14322,
                z: 0.0,
                alpha: vec![130.7093200, 23.8088610, 6.4436083],
                coef: vec![0.15432897, 0.53532814, 0.44463454],
            },
            SimintShell {
                am: 0,
                nprim: 3,
                x: 0.00000,
                y: -0.14322,
                z: 0.0,
                alpha: vec![5.0331513, 1.1695961, 0.3803890],
                coef: vec![-0.09996723, 0.39951283, 0.70011547],
            },
            SimintShell {
                am: 0,
                nprim: 3,
                x: 1.63804,
                y: 1.13654,
                z: 0.0,
                alpha: vec![3.42525091, 0.62391373, 0.16885540],
                coef: vec![0.15432897, 0.53532814, 0.44463454],
            },
            SimintShell {
                am: 0,
                nprim: 3,
                x: -1.63804,
                y: 1.13654,
                z: 0.0,
                alpha: vec![3.42525091, 0.62391373, 0.16885540],
                coef: vec![0.15432897, 0.53532814, 0.44463454],
            },
        ];
        let p_shells = vec![SimintShell {
            am: 1,
            nprim: 3,
            x: 0.0,
            y: -0.14322,
            z: 0.0,
            alpha: vec![5.0331513, 1.1695961, 0.3803890],
            coef: vec![0.15591627, 0.60768372, 0.39195739],
        }];
        let s_shells_simint: Vec<_> = s_shells.iter().map(|shell| shell.clone().to()).collect();
        let p_shells_simint: Vec<_> = p_shells.iter().map(|shell| shell.clone().to()).collect();
        let s_shells_roundtrip: Vec<_> = s_shells_simint
            .iter()
            .map(|shell| SimintShell::from(&shell.clone()))
            .collect();
        let p_shells_roundtrip: Vec<_> = p_shells_simint
            .iter()
            .map(|shell| SimintShell::from(&shell.clone()))
            .collect();
        assert_eq!(s_shells, s_shells_roundtrip);
        assert_eq!(p_shells, p_shells_roundtrip);
        let s_shells_normalized = normalize_shells(s_shells_roundtrip);
        let p_shells_normalized = normalize_shells(p_shells_roundtrip);
        unsafe {
            simint_finalize();
        };
    }
}
