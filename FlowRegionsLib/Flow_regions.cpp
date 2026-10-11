#include "Flow_regions.h"
#include "Earth_atmosphere.h"
#include "Oblique_shock.h"
#include "Prandtl_meyer.h"

#include <cmath>
#include <algorithm>

namespace Flow_regions {

    using Gas_dynamics::GasState;
    using Gas_dynamics::PI;

    DeltaSolution SolveDelta(const GasState& G3, const GasState& G4, double beta_deg) {
        const double k = Earth_atmosphere::GAMMA;

        DeltaSolution result;
        result.delta_deg = 0.0;

        const double tol = 1e-4;    // относительная разница p5 и p6
        const double step = 0.001;   // градусы
        const int    max_iter = 20000;

        for (int it = 0; it < max_iter; ++it) {
            double b5 = beta_deg + result.delta_deg;
            double b6 = beta_deg - result.delta_deg;

            auto S5 = Oblique_shock::SolveShockAngle(G3.M, b5);
            auto S6 = Oblique_shock::SolveShockAngle(G4.M, b6);
            if (!S5.ok || !S6.ok) { result.delta_deg += step; continue; }

            double tc5 = S5.theta_c_deg * PI / 180.0;
            double tc6 = S6.theta_c_deg * PI / 180.0;
            double Mn5 = G3.M * std::sin(tc5);
            double Mn6 = G4.M * std::sin(tc6);

            double p5 = G3.p * (2.0 * k / (k + 1.0) * Mn5 * Mn5 - (k - 1.0) / (k + 1.0));
            double p6 = G4.p * (2.0 * k / (k + 1.0) * Mn6 * Mn6 - (k - 1.0) / (k + 1.0));

            double rel = std::abs(p5 - p6) / std::max(std::abs(p5), 1.0);
            if (rel < tol) break;

            result.delta_deg += step;
        }

        result.G5 = Oblique_shock::JumpRelations(G3, beta_deg + result.delta_deg);
        result.G6 = Oblique_shock::JumpRelations(G4, beta_deg - result.delta_deg);
        return result;
    }

    FlowField ComputeAll(double M_inf, double p_inf, double T_inf,
        double alpha_deg, double beta_deg) {
        FlowField F;
        F.beta_deg = beta_deg;
        F.alpha_deg = alpha_deg;

        GasState Ginf = Gas_dynamics::FromMach(p_inf, T_inf, M_inf);

        double beta1 = beta_deg - alpha_deg;
        double beta2 = beta_deg + alpha_deg;

        // Области 1, 2 — за косыми скачками
        F.regions[0] = Oblique_shock::JumpRelations(Ginf, beta1);
        F.regions[1] = Oblique_shock::JumpRelations(Ginf, beta2);
        F.theta_c_deg[0] = F.regions[0].theta_shock;
        F.theta_c_deg[1] = F.regions[1].theta_shock;

        // Области 3, 4 — за веером Прандтля–Майера на угол 2β
        double turn = 2.0 * beta_deg;
        F.regions[2] = Prandtl_meyer::ExpandFlow(F.regions[0], turn);
        F.regions[3] = Prandtl_meyer::ExpandFlow(F.regions[1], turn);

        // Области 5, 6 — за системой скачков в донном следе
        DeltaSolution D = SolveDelta(F.regions[2], F.regions[3], beta_deg);
        F.regions[4] = D.G5;
        F.regions[5] = D.G6;
        F.delta_deg = D.delta_deg;

        F.beta_flow_deg[0] = beta1;
        F.beta_flow_deg[1] = beta2;
        F.beta_flow_deg[2] = 0.0;
        F.beta_flow_deg[3] = 0.0;
        F.beta_flow_deg[4] = beta_deg + D.delta_deg;
        F.beta_flow_deg[5] = beta_deg - D.delta_deg;

        F.theta_c_deg[2] = 0.0;
        F.theta_c_deg[3] = 0.0;
        F.theta_c_deg[4] = F.regions[4].theta_shock;
        F.theta_c_deg[5] = F.regions[5].theta_shock;

        return F;
    }
}