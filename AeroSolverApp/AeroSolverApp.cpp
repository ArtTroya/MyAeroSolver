#include "Earth_atmosphere.h"
#include "Gas_dynamics.h"
#include "Flow_regions.h"
#include "Aerodynamics.h"
#include "Boundary_layer.h"
#include "Report.h"

#include <iostream>
#include <iomanip>
#include <array>
#include <cmath>
#include <string>
#include <windows.h>

using namespace Earth_atmosphere;
using namespace Gas_dynamics;
using namespace Flow_regions;
using namespace Aerodynamics;
using namespace Boundary_layer;
using namespace Report;

int main() {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);

    // ============================================================
    //  НАСТРОЙКИ ВАРИАНТА
    // ============================================================

    const int    N = 18;
    const int    GROUP = 72;

    const std::string FILE_REPORT = "report18.csv";
    const std::string FILE_GRAPHS = "graphs18.csv";

    // ============================================================
    //  Формулы варианта
    // ============================================================

    double c_mm, b_mm;
    if (N >= 1 && N <= 6) { c_mm = 80.0; b_mm = 1000.0; }
    else if (N >= 7 && N <= 12) { c_mm = 90.0; b_mm = 1100.0; }
    else if (N >= 13 && N <= 18) { c_mm = 100.0; b_mm = 1200.0; }
    else if (N >= 19 && N <= 24) { c_mm = 110.0; b_mm = 1300.0; }
    else { c_mm = 100.0; b_mm = 1200.0; }

    double H_km, M_inf, T_wall;
    if (GROUP == 71) {
        H_km = 10.0 + 0.6 * N;
        M_inf = 4.0 + 0.17 * N;
        T_wall = 400.0;
    }
    else {
        H_km = 10.0 + 0.5 * N;
        M_inf = 4.5 + 0.12 * N;
        T_wall = 390.0;
    }

    const double alpha_deg = 2.0 + 0.1 * N;
    const double h = H_km * 1000.0;
    const double c_m = c_mm * 1e-3;
    const double b_m = b_mm * 1e-3;
    const double beta_deg = std::atan(c_m / b_m) * 180.0 / PI;
    const double L_face = std::sqrt(0.25 * (c_m * c_m + b_m * b_m));

    // ============================================================
    //  Расчёт
    // ============================================================

    AtmosphereParams A = Calculate(h);
    GasState Ginf = FromMach(A.p, A.T, M_inf);

    FlowField F = ComputeAll(M_inf, A.p, A.T, alpha_deg, beta_deg);
    AeroCoefficients aero = Compute(F, A.p, Ginf.rho, Ginf.v, c_m, b_m, alpha_deg);

    SurfaceBL upper = BuildSurfaceBL(L_face,
        F.regions[0], F.regions[2], T_wall, RE_CR_BASE);
    SurfaceBL lower = BuildSurfaceBL(L_face,
        F.regions[1], F.regions[3], T_wall, RE_CR_BASE);

    FaceResult face1 = SolveMixedFace(L_face,
        F.regions[0].v, F.regions[0].rho, F.regions[0].mu,
        F.regions[0].T, F.regions[0].M, F.regions[0].Cp, F.regions[0].lambda,
        T_wall, RE_CR_BASE);

    FaceResult face2 = SolveMixedFace(L_face,
        F.regions[1].v, F.regions[1].rho, F.regions[1].mu,
        F.regions[1].T, F.regions[1].M, F.regions[1].Cp, F.regions[1].lambda,
        T_wall, RE_CR_BASE);

    double delta_start_3 = face1.points[N_POINTS - 1].delta;
    double delta_start_4 = face2.points[N_POINTS - 1].delta;

    FaceResult face3 = SolveTurbulentFace(L_face,
        F.regions[2].v, F.regions[2].rho, F.regions[2].mu,
        F.regions[2].T, F.regions[2].M, F.regions[2].Cp, F.regions[2].lambda,
        T_wall, delta_start_3);

    FaceResult face4 = SolveTurbulentFace(L_face,
        F.regions[3].v, F.regions[3].rho, F.regions[3].mu,
        F.regions[3].T, F.regions[3].M, F.regions[3].Cp, F.regions[3].lambda,
        T_wall, delta_start_4);

    std::array<FaceResult, 4> faces = { face1, face2, face3, face4 };

    std::array<double, 4> q_local = {
        0.5 * F.regions[0].rho * F.regions[0].v * F.regions[0].v,
        0.5 * F.regions[1].rho * F.regions[1].v * F.regions[1].v,
        0.5 * F.regions[2].rho * F.regions[2].v * F.regions[2].v,
        0.5 * F.regions[3].rho * F.regions[3].v * F.regions[3].v
    };

    FrictionResult FR = CombineFriction(faces, q_local, beta_deg, aero.q_inf, b_m);

    double alpha_rad = alpha_deg * PI / 180.0;
    double Cx_sum = aero.Cx + FR.Cf_sigma;
    double Cxa_sum = aero.Cy * std::sin(alpha_rad) + Cx_sum * std::cos(alpha_rad);
    double Cya_sum = aero.Cy * std::cos(alpha_rad) - Cx_sum * std::sin(alpha_rad);

    Summary S{};
    S.q_inf = aero.q_inf;
    S.Cx = aero.Cx; S.Cy = aero.Cy; S.mz = aero.mz;
    S.Cxa = aero.Cxa; S.Cya = aero.Cya; S.K = aero.K; S.Cd = aero.Cd;
    S.Cf_sigma = FR.Cf_sigma;
    S.X_friction = FR.X_friction;
    S.Cx_sum = Cx_sum;
    S.Cxa_sum = Cxa_sum;
    S.Cya_sum = Cya_sum;
    S.Cf_avg = FR.Cf_avg;
    S.F_fr = FR.F_fr;

    std::cout << std::fixed << std::setprecision(4);

    std::cout << "Курсовая работа: обтекание ромбовидного профиля\n";
    std::cout << "Вариант N = " << N << ", группа СМЗ-" << GROUP << "\n";
    std::cout << "h = " << h << " м, M = " << M_inf
        << ", alpha = " << alpha_deg << " град, beta = " << beta_deg
        << " град, T_wall = " << T_wall << " K\n\n";

    std::cout << "Атмосфера\n";
    std::cout << "H = " << A.H << " м, T = " << A.T << " K, p = " << A.p
        << " Па, rho = " << A.rho << " кг/м^3, a = " << A.a
        << " м/с, g = " << A.g << " м/с^2\n\n";

    std::cout << "Область  M      p,Па       T,K       rho       V,м/с     theta_c  beta_flow\n";
    for (int i = 0; i < 6; ++i) {
        const auto& G = F.regions[i];
        std::cout << "  " << (i + 1) << "    "
            << std::setw(7) << G.M << "  "
            << std::setw(9) << G.p << "  "
            << std::setw(7) << G.T << "  "
            << std::setw(8) << G.rho << "  "
            << std::setw(8) << G.v << "  "
            << std::setw(7) << F.theta_c_deg[i] << "  "
            << std::setw(7) << F.beta_flow_deg[i] << "\n";
    }
    std::cout << "\nУгол delta = " << F.delta_deg << " град\n\n";

    std::cout << "АДХ без трения\n";
    std::cout << "Cx = " << aero.Cx << ", Cy = " << aero.Cy
        << ", mz = " << aero.mz << ", Cxa = " << aero.Cxa
        << ", Cya = " << aero.Cya << ", K = " << aero.K
        << ", Cd = " << aero.Cd << "\n\n";

    std::cout << "Re_cr и x_cr\n";
    std::cout << "  Верхняя: Re_cr = " << upper.Re_cr
        << ", x_cr = " << upper.x_cr << " м\n";
    std::cout << "  Нижняя:  Re_cr = " << lower.Re_cr
        << ", x_cr = " << lower.x_cr << " м\n\n";

    std::cout << "Средние Cf и силы трения\n";
    std::cout << "Грань  Cf_avg    F_тр,Н\n";
    for (int i = 0; i < 4; ++i) {
        std::cout << "  " << (i + 1) << "    "
            << std::setw(8) << FR.Cf_avg[i] << "  "
            << std::setw(8) << FR.F_fr[i] << "\n";
    }
    std::cout << "X_трение = " << FR.X_friction << " Н\n";
    std::cout << "Cf_sigma = " << FR.Cf_sigma << "\n\n";

    std::cout << "Итоговые АДХ\n";
    std::cout << "              без трения   с трением    разница,%\n";
    std::cout << "Cx  =  " << std::setw(11) << aero.Cx
        << std::setw(12) << Cx_sum
        << std::setw(12) << (Cx_sum / aero.Cx - 1.0) * 100.0 << "\n";
    std::cout << "Cxa =  " << std::setw(11) << aero.Cxa
        << std::setw(12) << Cxa_sum
        << std::setw(12) << (Cxa_sum / aero.Cxa - 1.0) * 100.0 << "\n";
    std::cout << "Cya =  " << std::setw(11) << aero.Cya
        << std::setw(12) << Cya_sum
        << std::setw(12) << (Cya_sum / aero.Cya - 1.0) * 100.0 << "\n\n";

    CaseInput in{};
    in.N = N; in.h = h; in.M_inf = M_inf;
    in.alpha_deg = alpha_deg; in.beta_deg = beta_deg;
    in.T_wall = T_wall; in.c_m = c_m; in.b_m = b_m; in.L_face = L_face;

    WriteReportCsv(FILE_REPORT, in, A, F, aero, upper, lower, S);
    WriteGraphsCsv(FILE_GRAPHS, upper, lower, L_face);

    std::cout << "Отчёт сохранён: " << FILE_REPORT << "\n";
    std::cout << "Данные для графиков: " << FILE_GRAPHS << "\n\n";
    std::cout << "Нажмите Enter для выхода...";
    std::cin.get();
    return 0;
}