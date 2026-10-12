#include "Report.h"
#include "Gas_dynamics.h"
#include <fstream>
#include <iomanip>

namespace Report {

    static void WriteBom(std::ofstream& f) {
        const unsigned char bom[3] = { 0xEF, 0xBB, 0xBF };
        f.write(reinterpret_cast<const char*>(bom), 3);
    }

    static void Sec(std::ofstream& f, const char* n) {
        f << "\n" << n << "\n";
    }

    void WriteReportCsv(const std::string& path,
        const CaseInput& in,
        const Earth_atmosphere::AtmosphereParams& atm,
        const Flow_regions::FlowField& F,
        const Aerodynamics::AeroCoefficients& aero,
        const Boundary_layer::SurfaceBL& upper,
        const Boundary_layer::SurfaceBL& lower,
        const Summary& S)
    {
        std::ofstream f(path, std::ios::out | std::ios::binary);
        if (!f) return;
        WriteBom(f);
        f << std::fixed << std::setprecision(4);
        f.imbue(std::locale::classic());

        Sec(f, "1. Исходные данные");
        f << "Параметр;Значение;Единица\n";
        f << "Вариант N;" << in.N << ";\n";
        f << "h;" << in.h << ";м\n";
        f << "M_inf;" << in.M_inf << ";\n";
        f << "alpha;" << in.alpha_deg << ";град\n";
        f << "beta;" << in.beta_deg << ";град\n";
        f << "T_wall;" << in.T_wall << ";K\n";
        f << "c;" << in.c_m << ";м\n";
        f << "b;" << in.b_m << ";м\n";
        f << "L грани;" << in.L_face << ";м\n";

        Sec(f, "2. Стандартная атмосфера");
        f << "Параметр;Значение;Единица\n";
        f << "H геопотенциальная;" << atm.H << ";м\n";
        f << "T;" << atm.T << ";K\n";
        f << "p;" << atm.p << ";Па\n";
        f << "rho;" << atm.rho << ";кг/м^3\n";
        f << "a;" << atm.a << ";м/с\n";
        f << "g;" << atm.g << ";м/с^2\n";

        Sec(f, "3. Параметры потока в 6 областях");
        f << "Область;M;p, 10^4 Па;p0, 10^7 Па;rho, кг/м^3;rho0, кг/м^3;"
            "T, K;T0, K;a, м/с;V, м/с;"
            "Cp, Дж/(кг*K);mu, 10^-5 Па*с;lambda, 10^-2 Вт/(м*K);"
            "theta_c, град;beta_flow, град;mu_mach, град\n";
        for (int i = 0; i < 6; ++i) {
            const auto& G = F.regions[i];
            f << (i + 1) << ";"
                << G.M << ";"
                << G.p * 1e-4 << ";"
                << G.p0 * 1e-7 << ";"
                << G.rho << ";"
                << G.rho0 << ";"
                << G.T << ";"
                << G.T0 << ";"
                << G.a << ";"
                << G.v << ";"
                << G.Cp << ";"
                << G.mu * 1e5 << ";"
                << G.lambda * 100.0 << ";"
                << F.theta_c_deg[i] << ";"
                << F.beta_flow_deg[i] << ";"
                << G.mu_mach << "\n";
        }
        f << "delta, град;" << F.delta_deg << "\n";

        Sec(f, "4. АДХ без учёта трения");
        f << "Параметр;Значение;Единица\n";
        f << "Cx;" << aero.Cx * 1e3 << ";10^-3\n";
        f << "Cy;" << aero.Cy * 1e2 << ";10^-2\n";
        f << "mz;" << aero.mz * 1e2 << ";10^-2\n";
        f << "Cxa;" << aero.Cxa * 1e3 << ";10^-3\n";
        f << "Cya;" << aero.Cya * 1e2 << ";10^-2\n";
        f << "K;" << aero.K << ";\n";
        f << "Cd;" << aero.Cd << ";\n";
        f << "q_inf;" << aero.q_inf << ";Па\n";

        Sec(f, "5. Температуры восстановления и определяющие");
        f << "№;Tr(lam), K;T*(lam), K;Tr(turb), K;T*(turb), K\n";
        f << "1;"
            << upper.ref_lam.Tr << ";" << upper.ref_lam.T_star << ";"
            << upper.ref_turb_1.Tr << ";" << upper.ref_turb_1.T_star << "\n";
        f << "2;"
            << lower.ref_lam.Tr << ";" << lower.ref_lam.T_star << ";"
            << lower.ref_turb_1.Tr << ";" << lower.ref_turb_1.T_star << "\n";
        f << "3;-;-;"
            << upper.ref_turb_3.Tr << ";" << upper.ref_turb_3.T_star << "\n";
        f << "4;-;-;"
            << lower.ref_turb_3.Tr << ";" << lower.ref_turb_3.T_star << "\n";

        Sec(f, "6. Re критическое и координата перехода");
        f << "Поверхность;Re_cr, 10^6;x_cr, 10^-1 м;"
            "delta_x_1, 10^-2 м;delta_x_3, 10^-2 м\n";
        f << "Верхняя;"
            << upper.Re_cr * 1e-6 << ";" << upper.x_cr * 10.0 << ";"
            << upper.delta_x_1 * 100.0 << ";" << upper.delta_x_3 * 100.0 << "\n";
        f << "Нижняя;"
            << lower.Re_cr * 1e-6 << ";" << lower.x_cr * 10.0 << ";"
            << lower.delta_x_1 * 100.0 << ";" << lower.delta_x_3 * 100.0 << "\n";

        Sec(f, "7. Параметры ПС по расчётным точкам");

        auto writeSurface = [&](const char* name,
            const Boundary_layer::SurfaceBL& s) {
                f << "\nПоверхность " << name << "\n";
                f << "N;x, м;x_фикт, м;тип;Re_x, 10^6;"
                    "delta, 10^-3 м;delta*, 10^-3 м;delta**, 10^-3 м;"
                    "tau, Па;Cfx, 10^-3;Cf_avg, 10^-3\n";
                for (int i = 0; i < Boundary_layer::N_POINTS; ++i) {
                    const auto& p = s.points[i];
                    f << (i + 1) << ";"
                        << p.x << ";" << p.x_fictive << ";"
                        << (p.type == 0 ? "лам" : "турб") << ";"
                        << p.Re_x * 1e-6 << ";"
                        << p.delta * 1e3 << ";"
                        << p.delta_star * 1e3 << ";"
                        << p.delta_star_star * 1e3 << ";"
                        << p.tau << ";"
                        << p.Cfx * 1e3 << ";"
                        << p.Cf_avg * 1e3 << "\n";
                }
            };

        writeSurface("верхняя", upper);
        writeSurface("нижняя", lower);

        Sec(f, "8. Средние коэффициенты трения");
        f << "Грань;Cf_avg, 10^-3;F_тр, Н\n";
        for (int i = 0; i < 4; ++i) {
            f << (i + 1) << ";"
                << S.Cf_avg[i] * 1e3 << ";"
                << S.F_fr[i] << "\n";
        }
        f << "X_трение, Н;" << S.X_friction << "\n";
        f << "Cf_sigma, 10^-3;" << S.Cf_sigma * 1e3 << "\n";

        Sec(f, "9. Итоговые АДХ профиля");
        f << "Параметр;Без трения;С трением;Изменение, %\n";
        f << "Cx, 10^-3;"
            << aero.Cx * 1e3 << ";"
            << S.Cx_sum * 1e3 << ";"
            << (S.Cx_sum / aero.Cx - 1.0) * 100.0 << "\n";
        f << "Cxa, 10^-3;"
            << aero.Cxa * 1e3 << ";"
            << S.Cxa_sum * 1e3 << ";"
            << (S.Cxa_sum / aero.Cxa - 1.0) * 100.0 << "\n";
        f << "Cya, 10^-2;"
            << aero.Cya * 1e2 << ";"
            << S.Cya_sum * 1e2 << ";"
            << (S.Cya_sum / aero.Cya - 1.0) * 100.0 << "\n";
    }

    void WriteGraphsCsv(const std::string& path,
        const Boundary_layer::SurfaceBL& upper,
        const Boundary_layer::SurfaceBL& lower,
        double L_face) {
        std::ofstream f(path, std::ios::out | std::ios::binary);
        if (!f) return;
        WriteBom(f);
        f << std::fixed << std::setprecision(4);
        f.imbue(std::locale::classic());

        auto writeSurf = [&](const char* name,
            const Boundary_layer::SurfaceBL& s) {
                f << name << "\n";
                f << "x_global, м;x_local, м;delta, 10^-3 м;"
                    "delta*, 10^-3 м;delta**, 10^-3 м;tau, Па;"
                    "Cfx, 10^-3;тип\n";
                for (int i = 0; i < Boundary_layer::N_POINTS; ++i) {
                    const auto& p = s.points[i];
                    double xg = (i < 6) ? p.x : (L_face + p.x);
                    f << xg << ";"
                        << p.x << ";"
                        << p.delta * 1e3 << ";"
                        << p.delta_star * 1e3 << ";"
                        << p.delta_star_star * 1e3 << ";"
                        << p.tau << ";"
                        << p.Cfx * 1e3 << ";"
                        << (p.type == 0 ? "лам" : "турб") << "\n";
                }
                f << "\n";
            };

        writeSurf("Верхняя поверхность", upper);
        writeSurf("Нижняя поверхность", lower);
    }
}