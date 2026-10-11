#pragma once
#include "Gas_dynamics.h"
#include <array>

namespace Flow_regions {

    struct FlowField {
        std::array<Gas_dynamics::GasState, 6> regions;  // 0..5 — области 1..6
        double beta_deg;                                // полураствор профиля
        double alpha_deg;                               // угол атаки
        double delta_deg;                               // угол линии тангенциального разрыва
        std::array<double, 6> beta_flow_deg;            // углы отклонения потока
        std::array<double, 6> theta_c_deg;              // углы наклона скачков (0 если нет)
    };

    struct DeltaSolution {
        double delta_deg;
        Gas_dynamics::GasState G5;
        Gas_dynamics::GasState G6;
    };

    // M_inf — число Маха набегающего потока
    // p_inf, T_inf — статические параметры набегающего потока
    // alpha_deg, beta_deg — углы атаки и полураствора, градусы
    FlowField ComputeAll(double M_inf, double p_inf, double T_inf,
        double alpha_deg, double beta_deg);

    // G3, G4 — состояния за веером разрежения, beta_deg — полураствор
    DeltaSolution SolveDelta(const Gas_dynamics::GasState& G3,
        const Gas_dynamics::GasState& G4,
        double beta_deg);
}