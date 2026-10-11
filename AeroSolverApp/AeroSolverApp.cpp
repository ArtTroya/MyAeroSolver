#include "Earth_atmosphere.h"
#include "Gas_dynamics.h"
#include "Flow_regions.h"
#include "Aerodynamics.h"
#include "Boundary_layer.h"
#include "Report.h"

#include <iostream>
#include <iomanip>
#include <array>
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

    const int    N = 18;
    const double h = (10.0 + 0.5 * N) * 1000.0;
    const double M_inf = 4.5 + 0.12 * N;
    const double alpha_deg = 2.0 + 0.1 * N;
    const double T_wall = 390.0;
    const double c_m = 0.100;
    const double b_m = 1.200;
    const double beta_deg = std::atan(c_m / b_m) * 180.0 / PI;
    const double L_face = std::sqrt(0.25 * (c_m * c_m + b_m * b_m));

    AtmosphereParams A = Calculate(h);
    GasState Ginf = FromMach(A.p, A.T, M_inf);

    FlowField F = ComputeAll(M_inf, A.p, A.T, alpha_deg, beta_deg);
    AeroCoefficients aero = Compute(F, A.p, Ginf.rho, Ginf.v, c_m, b_m, alpha_deg);

    // ---- Пограничный слой: 9 точек на каждую поверхность ----
    SurfaceBL upper = BuildSurfaceBL(L_face,
        F.regions[0], F.regions[2], T_wall, RE_CR_BASE);
    SurfaceBL lower = BuildSurfaceBL(L_face,
        F.regions[1], F.regions[3], T_wall, RE_CR_BASE);

    // ---- Силы трения и средние Cf по граням ----
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
    FrictionResult FR = CombineFriction(faces, L_face, beta_deg, aero.q_inf, b_m);

    // ---- Итоговые АДХ ----
    double alpha_rad = alpha_deg * PI / 180.0;
    double Cx_sum = aero.Cx + FR.Cf_sigma;
    double Cxa_sum = aero.Cy * std::sin(alpha_rad) + Cx_sum * std::cos(alpha_rad);
    double Cya_sum = aero.Cy * std::cos(alpha_rad) - Cx_sum * std::sin(alpha_rad);

    // ---- Сводка ----
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

    // ---- Консольный вывод ----
    std::cout << std::fixed << std::setprecision(4);

    std::cout << "Курсовая работа: обтекание ромбовидного профиля\n";
    std::cout << "Вариант N = " << N << ", группа СМЗ-72\n";
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

    std::cout << "Температуры восстановления и определяющие\n";
    std::cout << "Поверхность  Tr(lam)  T*(lam)  Tr(turb)  T*(turb)\n";
    std::cout << "  Верхняя     "
        << std::setw(7) << upper.ref_lam.Tr << "  "
        << std::setw(7) << upper.ref_lam.T_star << "  "
        << std::setw(7) << upper.ref_turb_3.Tr << "  "
        << std::setw(7) << upper.ref_turb_3.T_star << "\n";
    std::cout << "  Нижняя      "
        << std::setw(7) << lower.ref_lam.Tr << "  "
        << std::setw(7) << lower.ref_lam.T_star << "  "
        << std::setw(7) << lower.ref_turb_3.Tr << "  "
        << std::setw(7) << lower.ref_turb_3.T_star << "\n\n";

    std::cout << "Re_cr и x_cr\n";
    std::cout << "  Верхняя: Re_cr = " << upper.Re_cr
        << ", x_cr = " << upper.x_cr << " м\n";
    std::cout << "  Нижняя:  Re_cr = " << lower.Re_cr
        << ", x_cr = " << lower.x_cr << " м\n\n";

    auto printSurface = [&](const char* name, const SurfaceBL& s) {
        std::cout << "Поверхность " << name << "\n";
        std::cout << "  N   x,м     x_ф,м   тип   delta,м   tau,Па    Cfx\n";
        for (int i = 0; i < N_POINTS; ++i) {
            const auto& p = s.points[i];
            std::cout << "  " << (i + 1) << "  "
                << std::setw(6) << p.x << "  "
                << std::setw(6) << p.x_fictive << "  "
                << (p.type == 0 ? "лам" : "турб") << "  "
                << std::setw(8) << p.delta << "  "
                << std::setw(8) << p.tau << "  "
                << std::setw(7) << p.Cfx << "\n";
        }
        std::cout << "\n";
        };

    printSurface("верхняя", upper);
    printSurface("нижняя", lower);

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

    // ---- CSV ----
    CaseInput in{};
    in.N = N; in.h = h; in.M_inf = M_inf;
    in.alpha_deg = alpha_deg; in.beta_deg = beta_deg;
    in.T_wall = T_wall; in.c_m = c_m; in.b_m = b_m; in.L_face = L_face;

    WriteReportCsv("report.csv", in, A, F, aero, upper, lower, S);
    WriteGraphsCsv("graphs.csv", upper, lower, L_face);

    std::cout << "Отчёт сохранён: report.csv\n";
    std::cout << "Данные для графиков: graphs.csv\n\n";
    std::cout << "Нажмите Enter для выхода...";
    std::cin.get();
    return 0;
}