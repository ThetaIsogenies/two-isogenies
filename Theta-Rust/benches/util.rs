#[cfg(target_arch = "x86")]
pub fn core_cycles() -> u64 {
    use core::arch::x86::{_mm_lfence, _rdtsc};
    unsafe {
        _mm_lfence();
        _rdtsc()
    }
}

#[cfg(target_arch = "x86_64")]
pub fn core_cycles() -> u64 {
    use core::arch::x86_64::{_mm_lfence, _rdtsc};
    unsafe {
        _mm_lfence();
        _rdtsc()
    }
}

#[cfg(target_arch = "aarch64")]
pub fn core_cycles() -> u64 {
    use core::arch::asm;
    let mut x: u64;
    unsafe {
        asm!("dsb sy", "mrs {}, pmccntr_el0", out(reg) x);
    }
    x
}

use rand_core::{CryptoRng, RngCore};
use sha2::{Digest, Sha512};

// Fake RNG for benchmarks only. NOT ACTUALLY SECURE! DO NOT USE!
pub struct DRNG {
    buf: [u8; 64],
    ptr: usize,
}

impl DRNG {
    pub fn new() -> Self {
        Self::from_seed(&u64::MAX.to_le_bytes())
    }

    pub fn from_seed(seed: &[u8]) -> Self {
        let mut d = Self {
            buf: [0u8; 64],
            ptr: 0,
        };
        let mut sh = Sha512::new();
        sh.update(seed);
        d.buf[..].copy_from_slice(&sh.finalize());
        d
    }

    pub fn reseed(&mut self, seed: &[u8]) {
        let mut sh = Sha512::new();
        sh.update(&self.buf[32..]);
        sh.update(seed);
        self.buf[..].copy_from_slice(&sh.finalize());
        self.ptr = 0;
    }
}

impl RngCore for DRNG {
    fn next_u32(&mut self) -> u32 {
        let mut buf = [0u8; 4];
        self.fill_bytes(&mut buf);
        u32::from_le_bytes(buf)
    }

    fn next_u64(&mut self) -> u64 {
        let mut buf = [0u8; 8];
        self.fill_bytes(&mut buf);
        u64::from_le_bytes(buf)
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        let len = dest.len();
        let mut off = 0;
        while off < len {
            let mut clen = 32 - self.ptr;
            if clen > (len - off) {
                clen = len - off;
            }
            dest[off..off + clen].copy_from_slice(&self.buf[self.ptr..self.ptr + clen]);
            self.ptr += clen;
            off += clen;
            if self.ptr == 32 {
                let mut sh = Sha512::new();
                sh.update(&self.buf);
                self.buf[..].copy_from_slice(&sh.finalize());
                self.ptr = 0;
            }
        }
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

impl CryptoRng for DRNG {}
