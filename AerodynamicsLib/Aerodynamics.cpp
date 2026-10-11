#include "Aerodynamics.h"
#include "Earth_atmosphere.h"

#include <cmath>

namespace Aerodynamics {

    using Gas_dynamics::PI;

    AeroCoefficients Compute(const Flow_regions::FlowField& F,
        double p_inf, double rho_inf, double V_inf,
        double c, double b, double alpha_deg) {

        AeroCoefficients R{};

        // p1..p4 — давления в областях на гранях
        double dp1 = F.regions[0].p - p_inf;
        double dp2 = F.regions[1].p - p_inf;
        double dp3 = F.regions[2].p - p_inf;
        double dp4 = F.regions[3].p - p_inf;

        double X1 = dp1 * c / 2.0;
        double X2 = dp2 * c / 2.0;
        double X3 = dp3 * c / 2.0;
        double X4 = dp4 * c / 2.0;

        double Y1 = dp1 * b / 2.0;
        double Y2 = dp2 * b / 2.0;
        double Y3 = dp3 * b / 2.0;
        double Y4 = dp4 * b / 2.0;

        double Mz1 = Y1 * b / 4.0 + X1 * c / 4.0;
        double Mz2 = -Y2 * b / 4.0 - X2 * c / 4.0;
        double Mz3 = Y3 * 3.0 * b / 4.0 - X3 * c / 4.0;
        double Mz4 = -Y4 * 3.0 * b / 4.0 + X4 * c / 4.0;

        R.X = X1 + X2 - X3 - X4;
        R.Y = -Y1 + Y2 - Y3 + Y4;
        R.Mz = Mz1 + Mz2 + Mz3 + Mz4;

        R.q_inf = 0.5 * rho_inf * V_inf * V_inf;

        R.Cx = R.X / (R.q_inf * b);
        R.Cy = R.Y / (R.q_inf * b);
        R.mz = R.Mz / (R.q_inf * b * b);

        double alpha = alpha_deg * PI / 180.0;
        double sa = std::sin(alpha), ca = std::cos(alpha);

        R.Cxa = R.Cy * sa + R.Cx * ca;
        R.Cya = R.Cy * ca - R.Cx * sa;
        R.K = R.Cya / R.Cxa;
        R.Cd = -R.mz / R.Cy;

        return R;
    }
}