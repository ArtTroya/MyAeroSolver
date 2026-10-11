#pragma once
#include "Gas_dynamics.h"

namespace Prandtl_meyer {

    // M — число Маха, M >= 1
    // Возвращает фиктивный угол поворота потока omega, градусы
    double OmegaAngleDeg(double M);

    // omega_deg — угол в градусах
    // Возвращает число Маха, соответствующее этому углу
    double SolveMachFromOmega(double omega_deg);

    // G — состояние газа до расширения
    // delta_beta_deg — угол поворота потока (положительный)
    // Возвращает состояние газа за веером разрежения
    Gas_dynamics::GasState ExpandFlow(const Gas_dynamics::GasState& G,
        double delta_beta_deg);
}