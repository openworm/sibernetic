//! SPH Kernel Functions
//!
//! Smoothing kernels for SPH simulation. These functions determine how
//! particle properties are interpolated over space.
//!
//! ## References
//!
//! - Müller et al. (2003): Particle-based fluid simulation for interactive applications
//! - Monaghan (1992): Smoothed Particle Hydrodynamics

use std::f32::consts::PI;

/// Wpoly6 kernel - used for density calculation
///
/// W(r, h) = 315 / (64 * π * h^9) * (h² - r²)³  for 0 ≤ r ≤ h
///
/// Properties:
/// - Smooth at r = 0
/// - Zero at r = h
/// - Non-negative everywhere
#[inline]
pub fn wpoly6(r_sq: f32, h: f32, h_sq: f32, coeff: f32) -> f32 {
    if r_sq >= h_sq {
        return 0.0;
    }
    let diff = h_sq - r_sq;
    coeff * diff * diff * diff
}

/// Gradient of Wpoly6 kernel
///
/// ∇W(r, h) = -945 / (32 * π * h^9) * r * (h² - r²)²
#[inline]
pub fn grad_wpoly6(r: [f32; 3], r_len: f32, h: f32, h_sq: f32) -> [f32; 3] {
    let r_sq = r_len * r_len;
    if r_sq >= h_sq || r_len < 1e-10 {
        return [0.0, 0.0, 0.0];
    }

    let diff = h_sq - r_sq;
    let coeff = -945.0 / (32.0 * PI * h.powi(9));
    let factor = coeff * diff * diff;

    [factor * r[0], factor * r[1], factor * r[2]]
}

/// Wspiky kernel - used for pressure force calculation
///
/// W(r, h) = 15 / (π * h^6) * (h - r)³  for 0 ≤ r ≤ h
///
/// Properties:
/// - Has a sharp peak at r = 0
/// - Better for pressure forces (avoids particle clumping)
#[inline]
pub fn wspiky(r: f32, h: f32) -> f32 {
    if r >= h {
        return 0.0;
    }
    let coeff = 15.0 / (PI * h.powi(6));
    let diff = h - r;
    coeff * diff * diff * diff
}

/// Gradient of Wspiky kernel
///
/// ∇W(r, h) = -45 / (π * h^6) * (h - r)² * r̂
///
/// This gradient points from the neighbor toward the particle.
#[inline]
pub fn grad_wspiky(r: [f32; 3], r_len: f32, h: f32, coeff: f32) -> [f32; 3] {
    if r_len >= h || r_len < 1e-10 {
        return [0.0, 0.0, 0.0];
    }

    let diff = h - r_len;
    let factor = coeff * diff * diff / r_len;

    [factor * r[0], factor * r[1], factor * r[2]]
}

/// Precomputed gradient Wspiky coefficient
///
/// -45 / (π * h^6)
#[inline]
pub fn grad_wspiky_coeff(h: f32) -> f32 {
    -45.0 / (PI * h.powi(6))
}

/// Wviscosity kernel - used for viscosity force calculation
///
/// W(r, h) = 15 / (2 * π * h³) * (-r³/(2h³) + r²/h² + h/(2r) - 1)
///
/// Properties:
/// - Laplacian is non-negative everywhere (important for viscosity)
#[inline]
pub fn wviscosity(r: f32, h: f32) -> f32 {
    if r >= h || r < 1e-10 {
        return 0.0;
    }

    let h3 = h * h * h;
    let coeff = 15.0 / (2.0 * PI * h3);
    let r2 = r * r;
    let r3 = r2 * r;
    let h2 = h * h;

    coeff * (-r3 / (2.0 * h3) + r2 / h2 + h / (2.0 * r) - 1.0)
}

/// Laplacian of Wviscosity kernel
///
/// ∇²W(r, h) = 45 / (π * h^6) * (h - r)
///
/// This is always non-negative for r < h, which is important for stable viscosity.
#[inline]
pub fn lap_wviscosity(r: f32, h: f32, coeff: f32) -> f32 {
    if r >= h {
        return 0.0;
    }
    coeff * (h - r)
}

/// Precomputed Laplacian Wviscosity coefficient
///
/// 45 / (π * h^6)
#[inline]
pub fn lap_wviscosity_coeff(h: f32) -> f32 {
    45.0 / (PI * h.powi(6))
}

/// Cubic spline kernel (alternative to Wpoly6)
///
/// More commonly used in astrophysical SPH
#[inline]
pub fn cubic_spline(r: f32, h: f32) -> f32 {
    let q = r / h;
    let sigma = 1.0 / (PI * h * h * h);

    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let diff = 2.0 - q;
        sigma * 0.25 * diff * diff * diff
    } else {
        0.0
    }
}

/// Surface tension kernel (cohesion)
///
/// Used for modeling liquid surface tension effects
#[inline]
pub fn surface_tension_kernel(r: f32, h: f32) -> f32 {
    if r >= h {
        return 0.0;
    }

    let h2 = h * h;
    let r2 = r * r;
    let coeff = 32.0 / (PI * h.powi(9));

    if r < h / 2.0 {
        let inner = h.powi(6) / 64.0;
        coeff * 2.0 * (h2 - r2).powi(3) - inner
    } else {
        coeff * (h2 - r2).powi(3)
    }
}

/// Kernel for elastic force calculation
///
/// Simple spring-like force based on distance from rest length
#[inline]
pub fn elastic_force(
    r: f32,
    rest_length: f32,
    stiffness: f32,
) -> f32 {
    stiffness * (r - rest_length)
}

/// Boundary kernel for wall repulsion
///
/// Lennard-Jones style repulsion to keep particles inside boundaries
#[inline]
pub fn boundary_force(distance: f32, cutoff: f32, strength: f32) -> f32 {
    if distance >= cutoff || distance < 1e-10 {
        return 0.0;
    }

    let ratio = cutoff / distance;
    let ratio_6 = ratio.powi(6);
    let ratio_12 = ratio_6 * ratio_6;

    strength * (ratio_12 - ratio_6) / distance
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wpoly6_at_zero() {
        let h = 1.0;
        let h_sq = h * h;
        let coeff = 315.0 / (64.0 * PI * h.powi(9));

        let w = wpoly6(0.0, h, h_sq, coeff);
        assert!(w > 0.0, "Wpoly6 should be positive at r=0");
    }

    #[test]
    fn test_wpoly6_at_h() {
        let h = 1.0;
        let h_sq = h * h;
        let coeff = 315.0 / (64.0 * PI * h.powi(9));

        let w = wpoly6(h_sq, h, h_sq, coeff);
        assert_eq!(w, 0.0, "Wpoly6 should be zero at r=h");
    }

    #[test]
    fn test_wspiky_decreasing() {
        let h = 1.0;

        let w1 = wspiky(0.1, h);
        let w2 = wspiky(0.5, h);
        let w3 = wspiky(0.9, h);

        assert!(w1 > w2, "Wspiky should decrease with distance");
        assert!(w2 > w3, "Wspiky should decrease with distance");
    }

    #[test]
    fn test_lap_wviscosity_positive() {
        let h = 1.0;
        let coeff = lap_wviscosity_coeff(h);

        for r in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let lap = lap_wviscosity(r, h, coeff);
            assert!(lap >= 0.0, "Laplacian should be non-negative for r < h");
        }
    }
}
