#pragma once
// Exact Coulomb (gravitational) scattering wavefunction, eq. (154):
//
//   psi_k = e^{i k.r} Gamma(1 - i/(k r_b)) e^{pi/(2 k r_b)}
//                     M(i/(k r_b), 1, i k r (1 - khat.rhat))
//
// Units: k in 1/r_b, r in r_b (so k r_b -> k and r/r_b -> r below).
// Write xi_k = 1/(k r_b) = 1/k  (this is per-momentum, NOT the global CONST_XI)
// and u = k r (1 - khat.rhat) >= 0, so
//
//   psi = C(xi_k) e^{i k r cos} M(i xi_k, 1, i u),   C = Gamma(1-i xi_k) e^{pi xi_k/2}.
//
// Normalisation convention is the same as the old eq.(158) routine: an overall
// 1/sqrt(k r_b) is factored OUT, i.e. this returns sqrt(k r_b) * psi_exact, so
//   |C| = sqrt(2 pi / (1 - e^{-2 pi xi_k}))  ->  sqrt(2 pi)   (focused, xi_k >> 1)
//                                            ->  sqrt(k r_b) (unfocused, xi_k << 1)
// reproducing eq.(158) in the first case and a bare plane wave in the second.
//
// The only large parameter is |a z| = xi_k * u = (r/r_b)(1 - khat.rhat) <= 2 r/r_b,
// which is bounded by the integration volume, never by k. That is what makes the
// Kummer series usable: a large |z| always comes with a small |a| and vice versa.
// Validated against mpmath to <= 1e-7 relative for (r/r_b)(1-cos) <= 70.

#include <cmath>
#include <complex>

namespace coulomb {

using cplx = std::complex<double>;

// ---------------------------------------------------------------- log Gamma
// Lanczos (g=7, n=9). Always used inside exp(), so 2*pi*i branch offsets of the
// principal logs are harmless.
inline cplx log_gamma(cplx z) {
    static const double c[9] = {
        0.99999999999980993,    676.5203681218851,     -1259.1392167224028,
        771.32342877765313,     -176.61502916214059,   12.507343278686905,
        -0.13857109526572012,   9.9843695780195716e-6, 1.5056327351493116e-7};
    if (z.real() < 0.5)  // reflection formula
        return std::log(M_PI / std::sin(M_PI * z)) - log_gamma(1.0 - z);
    z -= 1.0;
    cplx x = c[0];
    for (int i = 1; i < 9; ++i) x += c[i] / (z + double(i));
    const cplx t = z + 7.5;
    return 0.5 * std::log(2 * M_PI) + (z + 0.5) * std::log(t) - t + std::log(x);
}

// Quantities that depend on a = i xi_k only; build once per momentum.
struct Context {
    cplx a{0.0, 1.0};
    cplx inv_gamma_a{1.0, 0.0};    // 1/Gamma(a)
    cplx inv_gamma_1ma{1.0, 0.0};  // 1/Gamma(1-a)
    cplx exp_ipi_a{1.0, 0.0};      // e^{+i pi a}
    cplx exp_mipi_a{1.0, 0.0};     // e^{-i pi a}

    Context() = default;
    explicit Context(cplx a_in) : a(a_in) {
        inv_gamma_a = std::exp(-log_gamma(a));
        inv_gamma_1ma = std::exp(-log_gamma(1.0 - a));
        exp_ipi_a = std::exp(cplx(0.0, M_PI) * a);
        exp_mipi_a = 1.0 / exp_ipi_a;
    }
};

// ------------------------------------------------- 1F1(a;1;z), Kummer series
// M(a,1,z) = sum_n (a)_n z^n / (n!)^2, convergent for every z. The relative
// accuracy is limited by cancellation, ~ max|term| / |sum|, reported through
// `cancellation`.
inline cplx M_b1_series(cplx a, cplx z, double* cancellation = nullptr) {
    cplx term(1.0, 0.0), sum(1.0, 0.0);
    double max_norm = 1.0;
    const double abs_z = std::abs(z), abs_az = std::abs(a * z);
    // terms grow until roughly n_peak; never test convergence before that
    const double n_peak = 0.5 * (abs_z + std::sqrt(abs_z * abs_z + 4.0 * abs_az));
    const int NMAX = 4000;
    const double EPS2 = 1e-34;  // (1e-17)^2, compare squared moduli
    for (int n = 0; n < NMAX; ++n) {
        const double inv = 1.0 / (double(n + 1) * double(n + 1));
        term *= (a + double(n)) * z * inv;
        sum += term;
        const double tn = std::norm(term);
        if (tn > max_norm) max_norm = tn;
        if (double(n) > n_peak + 2.0) {
            if (tn <= EPS2 * std::norm(sum)) break;
        } else if (tn > 1e300) {  // pathological |a z|: bail out, flagged below
            max_norm = 1e300;
            break;
        }
    }
    if (cancellation) *cancellation = std::sqrt(max_norm / std::norm(sum));
    return sum;
}

// ------------------------------------------ 1F1(a;1;z), large-|z| asymptotic
// DLMF 13.7.2 with b = 1:
//   M(a,1,z) ~ e^z z^{a-1}/Gamma(a)   sum_s [(1-a)_s]^2 / (s! z^s)
//            + e^{+-i pi a} z^{-a}/Gamma(1-a) sum_s [(a)_s]^2 / (s! (-z)^s),
// upper/lower sign for Im z >/< 0. `achieved` returns the size of the smallest
// term reached, i.e. the best relative accuracy available at these arguments.
inline cplx M_b1_asymptotic(cplx a, cplx z, const Context& ctx,
                            double* achieved = nullptr) {
    auto series = [](cplx p, cplx w, double* floor_out) {
        // sum_s [(p)_s]^2 / (s! w^s), truncated at the smallest term
        cplx term(1.0, 0.0), sum(1.0, 0.0);
        double prev = 1.0, smallest = 1.0;
        for (int s = 1; s <= 60; ++s) {
            const cplx pk = p + double(s - 1);
            const cplx next = term * pk * pk / (double(s) * w);
            const double at = std::abs(next);
            if (at > prev) break;  // asymptotic series turned around
            term = next;
            sum += term;
            prev = smallest = at;
            if (at < 1e-17) break;
        }
        if (floor_out) *floor_out = smallest;
        return sum;
    };
    double f1 = 1.0, f2 = 1.0;
    const cplx log_z = std::log(z);
    const cplx t1 = std::exp(z + (a - 1.0) * log_z) * ctx.inv_gamma_a *
                    series(1.0 - a, z, &f1);
    const cplx t2 = (z.imag() >= 0.0 ? ctx.exp_ipi_a : ctx.exp_mipi_a) *
                    std::exp(-a * log_z) * ctx.inv_gamma_1ma * series(a, -z, &f2);
    if (achieved) {
        const double d = std::abs(t1 + t2);
        *achieved = (std::abs(t1) * f1 + std::abs(t2) * f2) / (d > 0.0 ? d : 1.0);
    }
    return t1 + t2;
}

// ------------------------------------------------------------- dispatcher
inline cplx M_b1(const Context& ctx, cplx z, double* rel_err = nullptr) {
    const double abs_z = std::abs(z);
    if (abs_z < 15.0) {  // series is at machine precision here
        double cancel = 1.0;
        const cplx v = M_b1_series(ctx.a, z, &cancel);
        if (rel_err) *rel_err = cancel * 1e-16;
        return v;
    }
    if (abs_z > 60.0) {  // series cancellation hopeless, asymptotic is accurate
        return M_b1_asymptotic(ctx.a, z, ctx, rel_err);
    }
    // crossover band: both are usable, keep the one with the better estimate
    double cancel = 1.0, achieved = 1.0;
    const cplx vs = M_b1_series(ctx.a, z, &cancel);
    const cplx va = M_b1_asymptotic(ctx.a, z, ctx, &achieved);
    if (cancel * 1e-16 <= achieved) {
        if (rel_err) *rel_err = cancel * 1e-16;
        return vs;
    }
    if (rel_err) *rel_err = achieved;
    return va;
}

// --------------------------------------------------------- the wavefunction
// One of these per scattering momentum; build it when the momentum is drawn,
// then call psi() for every sample point.
struct ScatteringState {
    double k = 1.0;     // |k| in 1/r_b
    double xi_k = 1.0;  // 1/(k r_b), the focusing parameter of this momentum
    cplx prefactor{1.0, 0.0};
    Context ctx;

    ScatteringState() = default;
    explicit ScatteringState(double k_in) : k(k_in), xi_k(1.0 / k_in) {
        ctx = Context(cplx(0.0, xi_k));
        // |sqrt(k r_b) Gamma(1-i xi) e^{pi xi/2}| = sqrt(2 pi/(1-e^{-2 pi xi}))
        const double mag = std::sqrt(2.0 * M_PI / -std::expm1(-2.0 * M_PI * xi_k));
        const double phase = -std::imag(log_gamma(cplx(1.0, xi_k)));  // arg G(1-i xi)
        prefactor = mag * std::exp(cplx(0.0, phase));
    }

    // r = |r|/r_b, cos_kr = khat.rhat
    cplx psi(double r, double cos_kr, double* rel_err = nullptr) const {
        if (!(r > 0.0)) return prefactor;  // r == 0 (cos_kr is then undefined)
        const double kr = k * r;
        const cplx M = M_b1(ctx, cplx(0.0, kr * (1.0 - cos_kr)), rel_err);
        return prefactor * std::exp(cplx(0.0, kr * cos_kr)) * M;
    }
};

}  // namespace coulomb
