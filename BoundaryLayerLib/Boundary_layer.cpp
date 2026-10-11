#include "Boundary_layer.h"
#include "Earth_atmosphere.h"

#include <cmath>
#include <algorithm>

namespace Boundary_layer {

    using Gas_dynamics::PI;

    ReferenceState ComputeReference(double T_wall, double T_e, double M_e,
        double Cp_e, double mu_e, double lambda_e,
        bool turbulent) {
        ReferenceState R{};
        R.T_e = T_e;
        R.T_wall = T_wall;

        double r = turbulent ? 0.88 : 0.83;
        double T_star_old = 1e7;
        double T_star_new = 0.0;

        for (int it = 0; it < 200; ++it) {
            R.Tr = T_e * (1.0 + r * 0.2 * M_e * M_e);
            T_star_new = 0.5 * (T_wall + T_e) + 0.22 * (R.Tr - T_e);

            R.Cp_star = Cp_e * std::pow(T_star_new / 288.15, 0.1);
            R.mu_star = mu_e * std::pow(T_star_new / 288.15, 0.76);
            R.lambda_star = lambda_e * std::pow(T_star_new / 261.0, 0.85);
            R.Pr_star = R.Cp_star * R.mu_star / R.lambda_star;

            r = turbulent ? std::cbrt(R.Pr_star) : std::sqrt(R.Pr_star);

            if (std::abs(T_star_new - T_star_old) < 0.01 * T_star_new) {
                T_star_old = T_star_new;
                break;
            }
            T_star_old = T_star_new;
        }

        R.T_star = T_star_old;
        R.r = r;
        return R;
    }

    double CorrectedCriticalRe(double Re_cr_base, double T_wall, double Tr,
        double M_e) {
        double x = (T_wall / Tr - 1.0) / (M_e * M_e);
        double y = 1.0 - 16.0 * x - 412.5 * x * x
            - 35000.0 * x * x * x - 375000.0 * x * x * x * x;
        return y * Re_cr_base;
    }

    BLPoint ComputeLaminarPoint(double x, double V_e, double rho_e, double mu_e,
        const ReferenceState& ref) {
        BLPoint P{};
        P.x = x;
        P.x_fictive = 0.0;
        P.type = 0;

        if (x <= 0.0) {
            P.Re_x = 0.0;
            return P;
        }

        P.Re_x = V_e * rho_e * x / mu_e;
        double ratio = ref.T_star / ref.T_e;
        double corr_delta = std::pow(ratio, (N_MU_LAYER + 1.0) / 2.0);
        double corr_tau = std::pow(ratio, (N_MU_LAYER - 1.0) / 2.0);

        double sqRe = std::sqrt(P.Re_x);
        P.delta = 4.64 * x / sqRe * corr_delta;
        P.delta_star = 1.74 * x / sqRe * corr_delta;
        P.delta_star_star = 0.646 * x / sqRe * corr_delta;
        P.tau = 0.323 * rho_e * V_e * V_e / sqRe * corr_tau;
        P.Cfx = 0.646 / sqRe * corr_tau;
        P.Cf_avg = 1.3 / sqRe * corr_tau;
        return P;
    }

    BLPoint ComputeTurbulentPoint(double x, double x_fictive,
        double V_e, double rho_e, double mu_e,
        const ReferenceState& ref) {
        BLPoint P{};
        P.x = x;
        P.x_fictive = x_fictive;
        P.type = 1;

        if (x_fictive <= 0.0) {
            P.Re_x = 0.0;
            return P;
        }

        P.Re_x = V_e * rho_e * x_fictive / mu_e;
        double ratio = ref.T_star / ref.T_e;
        double corr_delta = std::pow(ratio, (N_MU_LAYER + 1.0) / 5.0);
        double corr_tau = std::pow(ratio, (N_MU_LAYER - 4.0) / 5.0);

        double Re15 = std::pow(P.Re_x, 0.2);
        P.delta = 0.37 * x_fictive / Re15 * corr_delta;
        P.delta_star = P.delta / (K_BAR + 1.0);
        P.delta_star_star = K_BAR * P.delta / ((K_BAR + 1.0) * (K_BAR + 2.0));
        P.tau = 0.0289 * rho_e * V_e * V_e / Re15 * corr_tau;
        P.Cfx = 0.0578 / Re15 * corr_tau;
        P.Cf_avg = 0.074 / Re15 * corr_tau;
        return P;
    }

    FaceResult SolveMixedFace(double L,
        double V_e, double rho_e, double mu_e,
        double T_e, double M_e,
        double Cp_e, double lambda_e,
        double T_wall, double Re_cr_base) {
        FaceResult F{};
        F.L = L;
        F.has_laminar = true;

        F.ref_lam = ComputeReference(T_wall, T_e, M_e, Cp_e, mu_e, lambda_e, false);
        F.ref_turb = ComputeReference(T_wall, T_e, M_e, Cp_e, mu_e, lambda_e, true);

        F.Re_cr = CorrectedCriticalRe(Re_cr_base, T_wall, F.ref_lam.Tr, M_e);
        F.x_cr = F.Re_cr * mu_e / (V_e * rho_e);
        if (F.x_cr > L) F.x_cr = L;

        BLPoint lam_cr = ComputeLaminarPoint(F.x_cr, V_e, rho_e, mu_e, F.ref_lam);

        double ratio_lam = F.ref_lam.T_star / T_e;
        double ratio_turb = F.ref_turb.T_star / T_e;

        F.delta_x = std::pow(
            std::pow(lam_cr.delta, 5.0) * V_e * rho_e
            / (mu_e * std::pow(0.37, 5.0)
                * std::pow(ratio_turb, N_MU_LAYER + 1.0)),
            0.25);

        const int lam_pts = 2;
        const int turb_pts = N_POINTS - lam_pts;

        F.points[0] = ComputeLaminarPoint(F.x_cr / 2.0, V_e, rho_e, mu_e, F.ref_lam);
        F.points[1] = lam_cr;

        double l_turb = L - F.x_cr;
        for (int i = 0; i < turb_pts; ++i) {
            double x_i = F.x_cr + l_turb * (i + 1) / turb_pts;
            double x_fict = F.delta_x + (x_i - F.x_cr);
            F.points[2 + i] = ComputeTurbulentPoint(x_i, x_fict,
                V_e, rho_e, mu_e, F.ref_turb);
        }

        double Cf_lam_avg = 1.3 / std::sqrt(F.Re_cr)
            * std::pow(ratio_lam, (N_MU_LAYER - 1.0) / 2.0);

        double Re_l = V_e * rho_e * (L - F.x_cr + F.delta_x) / mu_e;
        double Cf_turb_avg_1 = 0.074 / std::pow(Re_l, 0.2)
            * std::pow(ratio_turb, (N_MU_LAYER - 4.0) / 5.0);

        double Re_dx = V_e * rho_e * F.delta_x / mu_e;
        double Cf_turb_avg_2 = 0.074 / std::pow(Re_dx, 0.2)
            * std::pow(ratio_turb, (N_MU_LAYER - 4.0) / 5.0);

        F.Cf_avg = Cf_lam_avg * F.x_cr / L
            + Cf_turb_avg_1 * (L - F.x_cr + F.delta_x) / L
            - Cf_turb_avg_2 * F.delta_x / L;
        return F;
    }

    FaceResult SolveTurbulentFace(double L,
        double V_e, double rho_e, double mu_e,
        double T_e, double M_e,
        double Cp_e, double lambda_e,
        double T_wall,
        double delta_start) {
        FaceResult F{};
        F.L = L;
        F.has_laminar = false;
        F.ref_turb = ComputeReference(T_wall, T_e, M_e, Cp_e, mu_e, lambda_e, true);

        double ratio_turb = F.ref_turb.T_star / T_e;
        double Re_factor = V_e * rho_e / mu_e;

        F.delta_x = std::pow(
            std::pow(delta_start, 5.0) * Re_factor
            / (std::pow(0.37, 5.0) * std::pow(ratio_turb, N_MU_LAYER + 1.0)),
            0.25);

        for (int i = 0; i < N_POINTS; ++i) {
            double x_i = L * (i + 1) / N_POINTS;
            double x_fict = F.delta_x + x_i;
            F.points[i] = ComputeTurbulentPoint(x_i, x_fict,
                V_e, rho_e, mu_e, F.ref_turb);
        }

        double Re_l = V_e * rho_e * (F.delta_x + L) / mu_e;
        double Cf_turb_avg_1 = 0.074 / std::pow(Re_l, 0.2)
            * std::pow(ratio_turb, (N_MU_LAYER - 4.0) / 5.0);

        double Re_dx = V_e * rho_e * F.delta_x / mu_e;
        double Cf_turb_avg_2 = 0.074 / std::pow(Re_dx, 0.2)
            * std::pow(ratio_turb, (N_MU_LAYER - 4.0) / 5.0);

        F.Cf_avg = Cf_turb_avg_1 * (F.delta_x + L) / L
            - Cf_turb_avg_2 * F.delta_x / L;

        F.x_cr = 0.0;
        F.Re_cr = 0.0;
        return F;
    }

    SurfaceBL BuildSurfaceBL(double L_face,
        const Gas_dynamics::GasState& G1,
        const Gas_dynamics::GasState& G3,
        double T_wall, double Re_cr_base) {
        SurfaceBL S{};

        S.ref_lam = ComputeReference(T_wall, G1.T, G1.M, G1.Cp, G1.mu, G1.lambda, false);
        S.ref_turb_1 = ComputeReference(T_wall, G1.T, G1.M, G1.Cp, G1.mu, G1.lambda, true);
        S.ref_turb_3 = ComputeReference(T_wall, G3.T, G3.M, G3.Cp, G3.mu, G3.lambda, true);

        S.Re_cr = CorrectedCriticalRe(Re_cr_base, T_wall, S.ref_lam.Tr, G1.M);
        S.x_cr = S.Re_cr * G1.mu / (G1.v * G1.rho);
        if (S.x_cr > L_face) S.x_cr = L_face;

        // Точка 1: x = 0 — ламинарная (нулевая)
        S.points[0] = BLPoint{};
        S.points[0].x = 0.0;
        S.points[0].type = 0;

        // Точка 2: x = x_кр — ламинарная (конец ЛПС)
        S.points[1] = ComputeLaminarPoint(S.x_cr, G1.v, G1.rho, G1.mu, S.ref_lam);

        // Толщина ЛПС при x_кр — для расчёта Δx₁
        BLPoint lam_cr = S.points[1];

        // Δx₁ для грани 1
        double ratio1 = S.ref_turb_1.T_star / G1.T;
        S.delta_x_1 = std::pow(
            std::pow(lam_cr.delta, 5.0) * G1.v * G1.rho
            / (G1.mu * std::pow(0.37, 5.0)
                * std::pow(ratio1, N_MU_LAYER + 1.0)),
            0.25);

        // Точка 3: x = x_кр — начало ТПС (x_ф = Δx₁)
        S.points[2] = ComputeTurbulentPoint(S.x_cr, S.delta_x_1,
            G1.v, G1.rho, G1.mu, S.ref_turb_1);

        // Точки 4, 5, 6: турбулентные между x_кр и L
        for (int i = 1; i <= 3; ++i) {
            double x_i = S.x_cr + (L_face - S.x_cr) * i / 3.0;
            double x_f = S.delta_x_1 + (x_i - S.x_cr);
            S.points[2 + i] = ComputeTurbulentPoint(x_i, x_f,
                G1.v, G1.rho, G1.mu, S.ref_turb_1);
        }

        // Δx₃ из δ в конце грани 1
        double delta_end_1 = S.points[5].delta;
        double ratio3 = S.ref_turb_3.T_star / G3.T;
        S.delta_x_3 = std::pow(
            std::pow(delta_end_1, 5.0) * G3.v * G3.rho
            / (G3.mu * std::pow(0.37, 5.0)
                * std::pow(ratio3, N_MU_LAYER + 1.0)),
            0.25);

        // Точки 7, 8, 9: турбулентные на грани 3
        double xf[3] = { S.delta_x_3, S.delta_x_3 + L_face / 2.0, S.delta_x_3 + L_face };
        double xl[3] = { 0.0,          L_face / 2.0,              L_face };
        for (int i = 0; i < 3; ++i) {
            S.points[6 + i] = ComputeTurbulentPoint(xl[i], xf[i],
                G3.v, G3.rho, G3.mu, S.ref_turb_3);
        }

        return S;
    }

    FrictionResult CombineFriction(const std::array<FaceResult, 4>& faces,
        double L_face, double beta_deg,
        double q_inf, double b) {
        FrictionResult R{};
        double beta_rad = beta_deg * PI / 180.0;

        R.X_friction = 0.0;
        for (int i = 0; i < 4; ++i) {
            R.Cf_avg[i] = faces[i].Cf_avg;
            R.F_fr[i] = faces[i].Cf_avg * q_inf * L_face;
            R.X_friction += R.F_fr[i] * std::cos(beta_rad);
        }
        R.Cf_sigma = R.X_friction / (q_inf * b);
        return R;
    }
}