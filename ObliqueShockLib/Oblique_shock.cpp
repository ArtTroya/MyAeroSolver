#include "Oblique_shock.h"
#include "Earth_atmosphere.h"

#include <cmath>
#include <optional>

namespace Oblique_shock {

    using Gas_dynamics::GasState;
    using Gas_dynamics::PI;

    // Кубическое x^3 + a x^2 + b x + c = 0.
    // Возвращает наименьший положительный корень.
    static std::optional<double> SmallestPositiveRoot(double a, double b, double c) {
        double p = -a * a / 3.0 + b;
        double q = 2.0 * (a / 3.0) * (a / 3.0) * (a / 3.0) - a * b / 3.0 + c;
        double Q = (p / 3.0) * (p / 3.0) * (p / 3.0) + (q / 2.0) * (q / 2.0);

        if (Q < 0.0) {
            double arg = q / (2.0 * std::sqrt(-(p / 3.0) * (p / 3.0) * (p / 3.0)));
            if (arg > 1.0) arg = 1.0;
            if (arg < -1.0) arg = -1.0;
            double alpha = std::acos(arg);
            double x = 2.0 * std::sqrt(-p / 3.0)
                * std::cos(alpha / 3.0 + PI / 3.0)
                - a / 3.0;
            return x;
        }
        else if (std::abs(Q) < 1e-12) {
            double x = std::cbrt(-q / 2.0) - a / 3.0;
            return x;
        }
        return std::nullopt;
    }

    ShockSolution SolveShockAngle(double M, double beta_deg) {
        ShockSolution result{ 0.0, SolutionType::Detached, false };

        if (M <= 1.0 || beta_deg <= 0.0) return result;

        const double k = Earth_atmosphere::GAMMA;
        double beta = beta_deg * PI / 180.0;
        double tg_b = std::tan(beta);

        double a1 = ((k - 1.0) / 2.0 * M * M + 1.0) * tg_b;
        double b1 = 1.0 - M * M;
        double c1 = tg_b * (1.0 + (k + 1.0) / 2.0 * M * M);
        double d1 = 1.0;

        double a = b1 / a1;
        double b = c1 / a1;
        double c = d1 / a1;

        auto root = SmallestPositiveRoot(a, b, c);
        if (!root || *root <= 0.0) return result;

        result.theta_c_deg = std::atan(*root) * 180.0 / PI;
        result.type = SolutionType::AttachedWeak;
        result.ok = true;
        return result;
    }

    double MaxDeflectionAngleDeg(double M) {
        // beta_max достигается при Q = 0; ищем бисекцией
        double lo = 0.0, hi = 45.0;
        for (int i = 0; i < 60; ++i) {
            double mid = 0.5 * (lo + hi);
            auto s = SolveShockAngle(M, mid);
            if (s.ok) lo = mid; else hi = mid;
        }
        return lo;
    }

    GasState JumpRelations(const GasState& G, double beta_deg, double* theta_c_out) {
        const double k = Earth_atmosphere::GAMMA;

        ShockSolution sol = SolveShockAngle(G.M, beta_deg);
        if (!sol.ok) return G;                 // вернём как есть, если решения нет

        double theta_c = sol.theta_c_deg * PI / 180.0;
        double beta = beta_deg * PI / 180.0;

        double sin_t = std::sin(theta_c);
        double Mn = G.M * sin_t;
        double Mn2 = Mn * Mn;

        GasState out = G;
        out.theta_shock = sol.theta_c_deg;

        double p_ratio = 2.0 * k / (k + 1.0) * Mn2 - (k - 1.0) / (k + 1.0);
        double rho_ratio = (k + 1.0) * Mn2 / (2.0 + (k - 1.0) * Mn2);
        double T_ratio = p_ratio / rho_ratio;

        out.p = G.p * p_ratio;
        out.rho = G.rho * rho_ratio;
        out.T = G.T * T_ratio;
        out.a = std::sqrt(k * Earth_atmosphere::R_AIR * out.T);
        out.v = G.v * std::cos(theta_c) / std::cos(theta_c - beta);
        out.M = out.v / out.a;

        Gas_dynamics::FillTransport(out);
        Gas_dynamics::FillStagnation(out);

        if (theta_c_out) *theta_c_out = sol.theta_c_deg;
        return out;
    }
}