#pragma once
#include "Flow_regions.h"

namespace Aerodynamics {

    struct AeroCoefficients {
        double Cx;          // коэффициент продольной силы (связанная СК)
        double Cy;          // коэффициент нормальной силы
        double mz;          // коэффициент момента тангажа
        double Cxa;         // коэффициент лобового сопротивления (поточная СК)
        double Cya;         // коэффициент подъёмной силы (поточная СК)
        double K;           // аэродинамическое качество Cya/Cxa
        double Cd;          // коэффициент центра давления
        double q_inf;       // скоростной напор, Па
        double X;           // суммарная продольная сила, Н
        double Y;           // суммарная нормальная сила, Н
        double Mz;          // суммарный момент тангажа, Н*м
    };

    // F — поле течения (6 областей), p_inf, rho_inf, V_inf — параметры набегающего потока
    // c, b — толщина и хорда профиля в метрах
    // alpha_deg — угол атаки в градусах
    AeroCoefficients Compute(const Flow_regions::FlowField& F,
        double p_inf, double rho_inf, double V_inf,
        double c, double b, double alpha_deg);
}