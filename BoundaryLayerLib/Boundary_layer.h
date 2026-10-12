#pragma once
#include "Gas_dynamics.h"
#include <array>

namespace Boundary_layer {

    constexpr double K_BAR = 7.0;
    constexpr double N_MU_LAYER = 0.76;
    constexpr double RE_CR_BASE = 5.0e6;
    constexpr int    N_POINTS = 9;

    struct ReferenceState {
        double T_e;
        double T_wall;
        double Tr;
        double T_star;
        double Cp_star;
        double mu_star;
        double lambda_star;
        double Pr_star;
        double r;
    };

    struct BLPoint {
        double x;
        double x_fictive;
        double Re_x;
        int    type;
        // сжимаемый слой
        double delta;
        double delta_star;
        double delta_star_star;
        double tau;
        double Cfx;
        double Cf_avg;
        // несжимаемый слой (для сравнения)
        double Re_x_ns;
        double delta_ns;
        double delta_star_ns;
        double delta_star_star_ns;
        double tau_ns;
        double Cfx_ns;
        double Cf_ns;
    };

    struct FaceResult {
        double L;
        double Re_cr;
        double x_cr;
        double delta_x;
        double Cf_avg;
        double F_friction;
        ReferenceState ref_lam;
        ReferenceState ref_turb;
        std::array<BLPoint, N_POINTS> points;
        bool has_laminar;
    };

    struct SurfaceBL {
        ReferenceState ref_lam;
        ReferenceState ref_turb_1;
        ReferenceState ref_turb_3;
        double Re_cr;
        double x_cr;
        double delta_x_1;
        double delta_x_3;
        std::array<BLPoint, 9> points;
    };

    struct FrictionResult {
        std::array<double, 4> F_fr;
        std::array<double, 4> Cf_avg;
        double X_friction;
        double Cf_sigma;
    };

    ReferenceState ComputeReference(double T_wall, double T_e, double M_e,
        double Cp_e, double mu_e, double lambda_e,
        bool turbulent);

    double CorrectedCriticalRe(double Re_cr_base, double T_wall, double Tr,
        double M_e);

    BLPoint ComputeLaminarPoint(double x, double V_e, double rho_e, double mu_e,
        const ReferenceState& ref);

    BLPoint ComputeTurbulentPoint(double x, double x_fictive,
        double V_e, double rho_e, double mu_e,
        const ReferenceState& ref);

    FaceResult SolveMixedFace(double L,
        double V_e, double rho_e, double mu_e,
        double T_e, double M_e,
        double Cp_e, double lambda_e,
        double T_wall, double Re_cr_base);

    FaceResult SolveTurbulentFace(double L,
        double V_e, double rho_e, double mu_e,
        double T_e, double M_e,
        double Cp_e, double lambda_e,
        double T_wall,
        double delta_start);

    SurfaceBL BuildSurfaceBL(double L_face,
        const Gas_dynamics::GasState& G1,
        const Gas_dynamics::GasState& G3,
        double T_wall, double Re_cr_base);

    FrictionResult CombineFriction(const std::array<FaceResult, 4>& faces,
        const std::array<double, 4>& q_local,
        double beta_deg,
        double q_inf, double b);
}