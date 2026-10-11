#pragma once
#include "Gas_dynamics.h"

namespace Oblique_shock {

    enum class SolutionType {
        AttachedWeak,       // присоединённый слабый скачок (Q < 0)
        AttachedStrong,     // присоединённый сильный скачок (Q < 0)
        Detached            // присоединённого нет (Q > 0)
    };

    struct ShockSolution {
        double theta_c_deg;         // угол наклона скачка, градусы
        SolutionType type;
        bool ok;                    // true если решение найдено
    };

    // M — число Маха до скачка, beta_deg — угол отклонения потока, градусы
    ShockSolution SolveShockAngle(double M, double beta_deg);

    // M — число Маха; возвращает максимальный угол отклонения beta_max, градусы
    double MaxDeflectionAngleDeg(double M);

    // G — состояние газа до скачка, beta_deg — угол отклонения, градусы
    // Возвращает состояние газа за скачком (то же M получит новое значение)
    Gas_dynamics::GasState JumpRelations(const Gas_dynamics::GasState& G,
        double beta_deg,
        double* theta_c_out = nullptr);
}