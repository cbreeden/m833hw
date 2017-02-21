#![allow(unused_imports, unused_mut, dead_code, unused_variables, non_snake_case)]
extern crate rand;
use rand::Rng;
use rand::distributions::{IndependentSample, Range};

#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate itertools;

use itertools::Itertools;

mod exercise1 {
    use super::*;

    fn generate_xy(n: usize) -> (Vec<f64>, Vec<f64>) {
        let space = Range::new(0f64, 1f64);
        let mut rng = rand::thread_rng();

        let mut xs = Vec::with_capacity(n);
        let mut ys = Vec::with_capacity(n);

        for _ in 0..n {
            let t1 = space.ind_sample(&mut rng);
            let t2 = 2. * space.ind_sample(&mut rng);
            let t3 = space.ind_sample(&mut rng);
            let t4 = space.ind_sample(&mut rng);
            let t5 = space.ind_sample(&mut rng);
            let t6 = space.ind_sample(&mut rng);

            let x = (t1 + t4)
                .min(t1 + t3 + t5 + t6)
                .min(t2 + t5 + t6)
                .min(t2 + t3 + t4);

            xs.push(x);
            ys.push(t1 + t4);
        }

        (xs, ys)
    }

    pub fn anithetic1(n: usize) -> (Vec<f64>, Vec<f64>) {
        let space = Range::new(0f64, 1f64);
        let mut rng = rand::thread_rng();

        let mut xs = Vec::with_capacity(n);
        let mut ys = Vec::with_capacity(n);

        for _ in 0..n {
            let t1 = space.ind_sample(&mut rng);
            let t2 = space.ind_sample(&mut rng);
            let t3 = space.ind_sample(&mut rng);
            let t4 = space.ind_sample(&mut rng);
            let t5 = space.ind_sample(&mut rng);
            let t6 = space.ind_sample(&mut rng);

            let s1 = 1. - t2;
            let s2 = 2. * (1. - t1);
            let s3 = 1. - t3;
            let s4 = 1. - t4;
            let s5 = 1. - t5;
            let s6 = 1. - t6;

            let t2 = 2. * t2;

            let x = (t1 + t4)
                .min(t1 + t3 + t5 + t6)
                .min(t2 + t5 + t6)
                .min(t2 + t3 + t4);

            let y = (s1 + s4)
                .min(s1 + s3 + s5 + s6)
                .min(s2 + s5 + s6)
                .min(s2 + s3 + s4);

            xs.push(x);
            ys.push(y);
        }

        (xs, ys)
    }

    pub fn anithetic2(n: usize) -> (Vec<f64>, Vec<f64>) {
        let space = Range::new(0f64, 1f64);
        let mut rng = rand::thread_rng();

        let mut xs = Vec::with_capacity(n);
        let mut ys = Vec::with_capacity(n);

        for _ in 0..n {
            let t1 = space.ind_sample(&mut rng);
            let t2 = 2. * space.ind_sample(&mut rng);
            let t3 = space.ind_sample(&mut rng);
            let t4 = space.ind_sample(&mut rng);
            let t5 = space.ind_sample(&mut rng);
            let t6 = space.ind_sample(&mut rng);

            let s1 = 1. - t1;
            let s2 = 2. - t2;
            let s3 = 1. - t3;
            let s4 = 1. - t4;
            let s5 = 1. - t5;
            let s6 = 1. - t6;

            let x = (t1 + t4)
                .min(t1 + t3 + t5 + t6)
                .min(t2 + t5 + t6)
                .min(t2 + t3 + t4);

            let y = (s1 + s4)
                .min(s1 + s3 + s5 + s6)
                .min(s2 + s5 + s6)
                .min(s2 + s3 + s4);

            xs.push(x);
            ys.push(y);
        }

        (xs, ys)
    }

    pub fn run() {
        println!("# Exercise 1 #");
        let n = 100000;
        let (xs, ys) = generate_xy(n);

        println!("## Monte Carlo Method ##");
        let mc = xs.iter().fold(0.0, |acc, &x| acc + x) / (n as f64);
        let var = xs.iter()
            .map(|x| (x - mc).powi(2) )
            .sum::<f64>() / (n as f64 - 1.0);

        println!("Estimate: {}", mc);
        println!("Variance: {}", var);


        println!("");
        println!("## Using control variate ##");

        let alpha = 6.0 * xs.iter().zip(ys.iter())
            .map(|(x,y)| x*y - mc)
            .sum::<f64>() / (n as f64);

        let est_1 = xs.iter().zip(ys.iter())
            .map(|(x,y)| x - (y - 1.0))
            .sum::<f64>() / (n as f64);

        let est_a = xs.iter().zip(ys.iter())
            .map(|(x,y)| x - alpha * (y - 1.0))
            .sum::<f64>() / (n as f64);

        println!("Alpha estimated: {}", alpha);
        println!("Estimator with alpha: {}", est_a);
        println!("Estimator with 1: {}", est_1);

        let var_1 = xs.iter().zip(ys.iter())
            .map(|(x,y)| (x - (y - 1.0) - est_1)*(x - (y - 1.0) - est_1) )
            .sum::<f64>() / (n as f64 - 1.0);

        let var_a = xs.iter().zip(ys.iter())
            .map(|(x,y)| (x - alpha*(y - 1.0) - est_a)*(x - alpha*(y - 1.0) - est_a) )
            .sum::<f64>() / (n as f64 - 1.0);

        println!("Variance with 1: {}", var_1);
        println!("Variance with alpha {}", var_a);

        println!("");
        println!("## Antithetic ##");
        let (xs, ys) = anithetic1(n);
        let (xs2, ys2) = anithetic2(n);

        let est1 = xs.iter().chain(ys.iter())
            .sum::<f64>() / (2. * (n as f64));

        let est2 = xs2.iter().chain(ys2.iter())
            .sum::<f64>() / (2. * (n as f64));

        println!("Estimator1: {}", est1);
        println!("Estimator2: {}", est2);

        let var1 = xs.iter().zip(ys.iter())
            .map(|(x,y)| ((x + y)/ 2. - est1).powi(2))
            .sum::<f64>() / (n as f64 - 1.);

        let var2 = xs2.iter().zip(ys2.iter())
            .map(|(x,y)| ((x + y)/ 2. - est2).powi(2))
            .sum::<f64>() / (n as f64 - 1.);

        println!("Variance1: {}", var1);
        println!("Variance2: {}", var2);

        stratta(n);
    }

    fn stratta(n: usize) {
        let n = 3usize.pow(6)*200;
        let m = [0.0, 1.0/3.0, 2.0/3.0];
        let indices = iproduct!(m.iter(), m.iter(), m.iter(), m.iter(), m.iter(), m.iter());
        let j = n / 3usize.pow(6);

        let space = Range::new(0f64, 1f64);
        let mut rng = rand::thread_rng();

        let mut sum = 0f64;
        let mut var = 0f64;

        for (s1, s2, s3, s4, s5, s6) in indices {
            let mut sqr = 0f64;
            let mut sum_tmp = 0f64;
            for _ in 0..j {
                // generate new variables given strata indices
                let t1 = space.ind_sample(&mut rng) / 3.0 + s1;
                let t2 = 2.0 * space.ind_sample(&mut rng) / 3.0 + 2.0 * s2;
                let t3 = space.ind_sample(&mut rng) / 3.0 + s3;
                let t4 = space.ind_sample(&mut rng) / 3.0 + s4;
                let t5 = space.ind_sample(&mut rng) / 3.0 + s5;
                let t6 = space.ind_sample(&mut rng) / 3.0 + s6;

                let x = (t1 + t4)
                    .min(t1 + t3 + t5 + t6)
                    .min(t2 + t5 + t6)
                    .min(t2 + t3 + t4);

                sum_tmp += x;
                sqr += x * x;
            }

            sum += sum_tmp;
            var += (sqr - sum_tmp.powi(2) / (j as f64) ) / ((j - 1) as f64);
        }

        sum = sum * 1.0 / (3usize.pow(6) as f64) * 1.0 / (j as f64);
        println!("Strat estimator: {}", sum);
        println!("Strat Variance: {}", var);
    }


}


fn main() {
    exercise1::run();
}