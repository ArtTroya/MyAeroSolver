#include "Prandtl_meyer.h"
#include "Earth_atmosphere.h"

#include <cmath>

namespace Prandtl_meyer {

    using Gas_dynamics::GasState;
    using Gas_dynamics::PI;

    double OmegaAngleDeg(double M) {
        if (M < 1.0) M = 1.0;
        const double k = Earth_atmosphere::GAMMA;
        double s1 = std::sqrt((k + 1.0) / (k - 1.0));
        double s2 = std::sqrt((k - 1.0) / (k + 1.0));
        double term1 = s1 * std::atan(s2 * std::sqrt(M * M - 1.0));
        double term2 = std::atan(std::sqrt(M * M - 1.0));
        double omega_rad = term1 - term2;
        return omega_rad * 180.0 / PI;
    }

    double SolveMachFromOmega(double omega_deg) {
        // Бисекция по M в диапазоне [1, 50]
        double target = omega_deg;
        double lo = 1.0, hi = 50.0;

        // Если omega меньше нуля — поток не расширяется
        if (target <= 0.0) return 1.0;

        for (int i = 0; i < 80; ++i) {
            double mid = 0.5 * (lo + hi);
            double w = OmegaAngleDeg(mid);
            if (w < target) lo = mid;
            else            hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    GasState ExpandFlow(const GasState& G, double delta_beta_deg) {
        const double k = Earth_atmosphere::GAMMA;

        double omega_in = OmegaAngleDeg(G.M);
        double omega_out = omega_in + delta_beta_deg;
        double M_out = SolveMachFromOmega(omega_out);

        GasState out = G;
        out.M = M_out;

        // Изоэнтропические соотношения через параметры торможения
        double denom = 1.0 + 0.5 * (k - 1.0) * M_out * M_out;
        out.p = G.p0 / std::pow(denom, k / (k - 1.0));
        out.rho = G.rho0 / std::pow(denom, 1.0 / (k - 1.0));
        out.T = G.T0 / denom;
        out.a = std::sqrt(k * Earth_atmosphere::R_AIR * out.T);
        out.v = out.M * out.a;

        Gas_dynamics::FillTransport(out);
        Gas_dynamics::FillStagnation(out);

        return out;
    }
}