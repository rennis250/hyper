#![feature(cfg_target_feature)]
extern crate simd;
extern crate lodepng;
extern crate nalgebra;
extern crate csv;
extern crate byteorder;

use std::path::Path;
use std::fs::File;
use std::io::{BufReader, BufRead, Read, Write, Cursor};
use std::iter::FromIterator;

use byteorder::{BigEndian, LittleEndian, ReadBytesExt, WriteBytesExt};

use nalgebra::Matrix3;

#[cfg(target_feature = "sse2")]
use simd::x86::sse2::*;

mod hyper {
    use super::*;

    #[derive(Debug)]
    pub struct Hyper {
        data: Vec<u16>,
        cal_data: Vec<f64x2>,
        pub xyz: Vec<f64>,
        pub rgb: Vec<u16>,
        rad_corr: Vec<f64x2>,
        pub lines: usize,
        pub samples: usize,
        bands: usize,
        adsline: usize,
        specb: usize,
        spatb: usize,
        tint: f64,
        wln: Vec<f64>,
        rads: f64,
        erads: f64,
        pub fname: String,
    }

    impl Hyper {
        pub fn new() -> Hyper {
            Hyper {
                data: vec![0u16; 0],
                cal_data: vec![f64x2::splat(0f64); 0],
                xyz: vec![0f64; 0],
                rgb: vec![0u16; 0],
                rad_corr: vec![f64x2::splat(0f64); 0],
                lines: 0usize,
                samples: 0usize,
                bands: 0usize,
                adsline: 0usize,
                specb: 0usize,
                spatb: 0usize,
                tint: 0f64,
                wln: vec![0f64; 0],
                rads: 1000f64,
                erads: 100000f64,
                fname: "".to_string(),
            }
        }

        pub fn parse_header(&mut self) {
            let mut wlnc = 0;

            let mut processing_wavelengths = false;
            let file = File::open(&self.fname).unwrap();
            for line in BufReader::new(file).lines() {
                let l = line.unwrap();
                if !processing_wavelengths {
                    if l.starts_with("samples = ") {
                        self.samples =
                            l.split_whitespace().nth(2).unwrap().parse::<usize>().unwrap();
                    }
                    if l.starts_with("lines   = ") {
                        self.lines = l.split_whitespace().nth(2).unwrap().parse::<usize>().unwrap();
                    }
                    if l.starts_with("bands = ") {
                        self.bands = l.split_whitespace().nth(2).unwrap().parse::<usize>().unwrap();
                        self.wln = vec![0f64; self.bands];
                    }
                    if l.starts_with("autodarkstartline   = ") {
                        self.adsline =
                            l.split_whitespace().nth(2).unwrap().parse::<usize>().unwrap();
                    }
                    if l.starts_with("binning = {") {
                        self.specb = l.split(',')
                            .nth(0)
                            .unwrap()
                            .split('{')
                            .nth(1)
                            .unwrap()
                            .trim()
                            .parse::<usize>()
                            .unwrap();
                        self.spatb = l.split(',')
                            .nth(1)
                            .unwrap()
                            .split('}')
                            .nth(0)
                            .unwrap()
                            .trim()
                            .parse::<usize>()
                            .unwrap();
                    }
                    if l.starts_with("tint =  ") {
                        self.tint = l.split_whitespace().nth(2).unwrap().parse::<f64>().unwrap();
                    }
                    if l.starts_with("Wavelength = { ") {
                        processing_wavelengths = true;
                    }
                } else {
                    if l.starts_with("}") {
                        break;
                    }
                    self.wln[wlnc] = l.split(',').nth(0).unwrap().trim().parse::<f64>().unwrap();
                    wlnc += 1;
                }
            }
        }

        pub fn read_cal(&mut self) {
            let mut samples: usize = 0;
            let mut bands: usize = 0;

            let file =
                File::open(&Path::new("/Users/rje/my_docs/projects/hyper/Radiometric_calibration_1x1_Scanner.hdr")).unwrap();
            for line in BufReader::new(file).lines() {
                let l = line.unwrap();
                if l.starts_with("samples = ") {
                    samples = l.split_whitespace().nth(2).unwrap().parse::<usize>().unwrap();
                }
                if l.starts_with("bands = ") {
                    bands = l.split_whitespace().nth(2).unwrap().parse::<usize>().unwrap();
                }
            }

            let mut cald = vec![0f64; samples * bands];

            let mut file =
                File::open(&Path::new("/Users/rje/my_docs/projects/hyper/Radiometric_calibration_1x1_Scanner.cal"))
                    .unwrap();
            let file_len = file.metadata().unwrap().len();
            let mut bytes: Vec<u8> = Vec::with_capacity(file_len as usize + 1);
            file.read_to_end(&mut bytes).unwrap();
            let mut rdr = Cursor::new(bytes);
            for x in cald.iter_mut() {
                *x = rdr.read_f32::<LittleEndian>().unwrap() as f64;
            }

            cald[0] = cald[samples];
            cald[1] = cald[samples + 1];

            if samples == 1600 && bands == 800 && self.specb != 1 && self.spatb != 1 {
                let bbands = bands / self.specb;
                let bsamples = samples / self.spatb;

                let mut binned_cal_file = vec![0f64; bsamples * bbands];

                let mut spatc = 0;
                let mut specc = 0;
                for wln in (0..bands).filter(|x| (x % self.specb == 0)) {
                    for y in (0..samples).filter(|x| (x % self.spatb == 0)) {
                        let mut accum = 0.0;
                        for wlns in 0..self.specb {
                            for ys in 0..self.spatb {
                                accum += cald[(wln + wlns) * samples + y + ys];
                            }
                        }
                        binned_cal_file[specc * bsamples + spatc] =
                            accum / (self.spatb * self.specb) as f64;
                        spatc += 1;
                    }
                    spatc = 0;
                    specc += 1;
                }

                // self.rad_corr = binned_cal_file;

                let mut c = 0;
                self.rad_corr = vec![f64x2::splat(0.0f64); bsamples * bbands/2];
                for x in (0..(bsamples * bbands)).filter(|x| (x % 2 == 0)) {
                    self.rad_corr[c] = f64x2::new(binned_cal_file[x], binned_cal_file[x + 1]);
                    c += 1;
                }

                return;
            }

            // self.rad_corr = cald;

            let mut c = 0;
            self.rad_corr = vec![f64x2::splat(0.0f64); samples * bands/2];
            for x in (0..(samples * bands)).filter(|x| (x % 2 == 0)) {
                self.rad_corr[c] = f64x2::new(cald[x], cald[x + 1]);
                c += 1;
            }
        }

        pub fn read_raw(&mut self) {
            let nelems = self.samples * self.bands / 2 * self.adsline;
            self.cal_data = vec![f64x2::splat(0f64); nelems];

            let nelems = self.samples * self.bands * self.lines;
            self.data = vec![0u16; nelems];

            let mut file = File::open(self.fname[0..(self.fname.len() - 4)].to_string() + ".raw")
                .unwrap();
            let file_len = file.metadata().unwrap().len();
            let mut bytes: Vec<u8> = Vec::with_capacity(file_len as usize + 1);
            file.read_to_end(&mut bytes).unwrap();
            let mut rdr = Cursor::new(bytes);
            for x in self.data.iter_mut() {
                *x = rdr.read_u16::<LittleEndian>().unwrap();
            }
            self.data[0] = 0;
        }

        pub fn calibrate_raw(&mut self) {
            let corr_factor =
                f64x2::splat(1.0e-5 / (self.tint / (self.rads / (self.spatb * self.specb) as f64)));

            self.rad_corr = self.rad_corr.iter().map(|&x| x * corr_factor).collect();

            let mut idx = 0;
            let mut accum = f64x2::splat(0.0f64);
            let den = f64x2::splat((self.lines - self.adsline - 1) as f64);
            let mut y2 = 0;
            for y in (0..self.bands).filter(|x| (x % 2 == 0)) {
                for z in 0..self.samples {
                    accum = f64x2::splat(0.0f64);
                    for x in (0..self.lines).rev() {
                        idx = x * self.bands * self.samples + y * self.samples + z;
                        let v = f64x2::new(self.data[idx] as f64, self.data[idx + 1] as f64);

                        if x >= self.adsline + 1 {
                            accum = accum + v;
                        } else if x == self.adsline {
                            continue;
                        } else {
                            idx = x * self.bands / 2 * self.samples + y2 * self.samples + z;
                            self.cal_data[idx] = self.rad_corr[y2 * self.samples + z] *
                                                 (v - (accum / den));
                        }
                    }
                }
                y2 += 1;
            }

            self.cal_data = self.cal_data
                .iter()
                .map(|&x| if x.lt(f64x2::splat(0.0f64)).any() {
                    x.max(f64x2::splat(0.0))
                } else {
                    x
                })
                .collect();

            self.lines = self.adsline;

            return;
        }

        // pub fn write_dat(&self) {
        //     let mut file = File::create(self.fname[0..(self.fname.len() - 4)].to_string() + ".dat")
        //         .unwrap();
        //     let mut wtr = vec![];
        //     for x in &self.cal_data {
        //         wtr.write_u16::<LittleEndian>((x * self.erads) as u16).unwrap();
        //     }
        //     file.write_all(&mut wtr).unwrap();
        //     file.flush().unwrap();
        // }

        // pub fn xyz_to_rgb(&mut self) {
        //     let mut monxyYx: Vec<f64> = Vec::new();
        //     let mut monxyYy: Vec<f64> = Vec::new();
        //     let mut monxyYY: Vec<f64> = Vec::new();

        //     let mut rdr =
        //         csv::Reader::from_file(&Path::new("/Users/rje/my_docs/projects/hyper/OLEDxyY.csv"))
        //             .unwrap()
        //             .has_headers(false);
        //     for record in rdr.decode() {
        //         let (x, y, Y): (f64, f64, f64) = record.unwrap();
        //         monxyYx.push(x);
        //         monxyYy.push(y);
        //         monxyYY.push(Y);
        //     }

        //     let mut monxyzx = vec![0f64; 3];
        //     let mut monxyzy = vec![0f64; 3];
        //     let mut monxyzz = vec![0f64; 3];

        //     for i in 0..3 {
        //         monxyzx[i] = (monxyYY[i] / monxyYy[i]) * monxyYx[i];
        //         monxyzy[i] = monxyYY[i];
        //         monxyzz[i] = (monxyYY[i] / monxyYy[i]) * (1.0 - monxyYx[i] - monxyYy[i]);
        //     }

        //     let matrix = Matrix3::new(monxyzx[0],
        //                               monxyzx[1],
        //                               monxyzx[2],
        //                               monxyzy[0],
        //                               monxyzy[1],
        //                               monxyzy[2],
        //                               monxyzz[0],
        //                               monxyzz[1],
        //                               monxyzz[2]);
        //     let monXYZ2RGB = matrix.try_inverse().unwrap();

        //     let mut rgb = vec![f64x2::splat(0f64); self.lines * (self.samples / 2) * 3];

        //     let mut idx = 0;
        //     let A = vec![1.0267, 0.9983, 1.0337];
        //     let B = vec![1.0/2.2546, 1.0/2.2473, 1.0/2.2121];
        //     let mut c = f64x2::splat(0f64);
        //     let mut h = vec![0f64; 2];
        //     let mut y2 = 0;
        //     for z in 0..3 {
        //         y2 = 0;
        //         for y in (0..self.samples).filter(|x| (x % 2 == 0)) {
        //             for x in 0..self.lines {
        //                 idx = 3 * self.lines * y + 3 * x;
        //                 c = monXYZ2RGB[(z, 0)] * self.xyz[idx] +
        //                     monXYZ2RGB[(z, 1)] * self.xyz[idx + 1] +
        //                     monXYZ2RGB[(z, 2)] * self.xyz[idx + 2];

        //                 c.store(&mut h, 0);
        //                 h.iter_mut().map(|x| *x.powf(B[z]));

        //                 idx = 3 * self.lines * y2 + 3 * x;
        //                 rgb[idx + z] = h / A[z];
        //             }
        //             y2 += 1;
        //         }
        //     }

        //     self.rgb = rgb.iter()
        //         .map(|&x| if x.lt(f64x2::splat(0.0f64)).any() {
        //             x.max(f64x2::splat(0.0))
        //         } else {
        //             x
        //         })
        //         .collect();
        //         .map(|&x| if x < 0.0 {
        //             0.0 as u16
        //         } else {
        //             (x * 65535f64) as u16
        //         })
        //         .collect();
        // }

        pub fn raw_to_rgb(&mut self) {
            let mut ciex: Vec<f64x2> = Vec::new();
            let mut ciey: Vec<f64x2> = Vec::new();
            let mut ciez: Vec<f64x2> = Vec::new();

            let mut rdr =
                csv::Reader::from_file(&Path::new("/Users/rje/my_docs/projects/hyper/corrCMF.csv"))
                    .unwrap()
                    .has_headers(false);

            let x = 0;
            for record in rdr.decode() {
                let (x, y, z): (f64, f64, f64) = record.unwrap();
                ciex.push(f64x2::new(x));
                ciey.push(y);
                ciez.push(z);
            }

            ciex = ciex.iter()
                .map(|&x| if x < 0.0 { 0.0 } else { x })
                .collect();

            ciey = ciey.iter()
                .map(|&x| if x < 0.0 { 0.0 } else { x })
                .collect();

            ciez = ciez.iter()
                .map(|&x| if x < 0.0 { 0.0 } else { x })
                .collect();

            let mut xyz = vec![f64x2::splat(0f64); self.lines * self.samples * 2];

            let mut idx = 0;
            let mut ax1 = 0f64;
            let mut ay1 = 0f64;
            let mut az1 = 0f64;
            let mut ax2 = 0f64;
            let mut ay2 = 0f64;
            let mut az2 = 0f64;
            for x in 0..self.lines {
                for z in 0..self.samples {
                    ax1 = 0f64;
                    ay1 = 0f64;
                    az1 = 0f64;
                    ax2 = 0f64;
                    ay2 = 0f64;
                    az2 = 0f64;
                    for (y, v) in ciex.iter().enumerate() {
                        idx = x * (self.bands / 2) * self.samples + (y + 4) * self.samples + z;
                        ax1 += *v * self.cal_data[idx];
                        ay1 += ciey[y] * self.cal_data[idx];
                        az1 += ciez[y] * self.cal_data[idx];

                        idx = x * (self.bands / 2) * self.samples + (y + 4) * self.samples + z + 1;
                        ax2 += *v * self.cal_data[idx];
                        ay2 += ciey[y] * self.cal_data[idx];
                        az2 += ciez[y] * self.cal_data[idx];
                    }
                    idx = 3 * self.lines * z + 3 * x;
                    xyz[idx] = accumx;
                    xyz[idx + 1] = accumy;
                    xyz[idx + 2] = accumz;
                }
            }

            self.xyz = xyz.iter()
                .map(|&x| if x.lt(f64x2::splat(0.0f64)).any() {
                    x.max(f64x2::splat(0.0))
                } else {
                    x
                })
                .collect();
            // self.xyz_to_rgb()
        }
    }
}

fn main() {
    let mut h = hyper::Hyper::new();

    h.fname = std::env::args().nth(1).unwrap();

    h.parse_header();
    h.read_cal();
    h.read_raw();
    h.calibrate_raw();
    // h.write_dat();
    h.raw_to_rgb();

    // let mut rgb = vec![];
    // for x in &h.rgb {
    //     rgb.write_u16::<BigEndian>(*x).unwrap();
    // }
    // if let Err(e) = lodepng::encode_file(&Path::new(&(h.fname[0..(h.fname.len() - 4)].to_string() +
    //                                                   ".png")),
    //                                      &rgb,
    //                                      h.lines,
    //                                      h.samples,
    //                                      lodepng::LCT_RGB,
    //                                      16) {
    //     panic!("failed to write png: {:?}", e);
    // }

    // let xyz2: Vec<u16> = h.xyz
    //     .iter()
    //     .map(|&x| (x / f64x2::splat(100.0) * f64x2::splat(65535f64)) as u16)
    //     .collect();

    // let mut xyz = vec![];
    // for x in &xyz2 {
    //     xyz.write_u16::<BigEndian>(*x).unwrap();
    // }
    // if let Err(e) = lodepng::encode_file(&Path::new(&(h.fname[0..(h.fname.len() - 4)].to_string() +
    //                                                   "XYZ.png")),
    //                                      &xyz,
    //                                      h.lines,
    //                                      h.samples,
    //                                      lodepng::LCT_RGB,
    //                                      16) {
    //     panic!("failed to write png: {:?}", e);
    // }
}
