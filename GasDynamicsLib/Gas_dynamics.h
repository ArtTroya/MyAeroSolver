#pragma once

namespace Gas_dynamics {

    constexpr double PI = 3.14159265358979323846;

    // Опорные значения свойств воздуха
    constexpr double Cp_REF = 1004.7;      // Дж/(кг*K) при T_NORM
    constexpr double MU_REF = 1.79e-5;     // Па*с при T_NORM
    constexpr double LAMBDA_REF = 0.0232;      // Вт/(м*K) при T_LAMBDA_REF
    constexpr double T_NORM = 288.15;      // K
    constexpr double T_LAMBDA_REF = 261.0;       // K

    // Показатели степенных зависимостей
    constexpr double PHI_CP = 0.1;
    constexpr double N_MU = 0.76;
    constexpr double XI_LAMBDA = 0.85;

    struct GasState {
        double p;          // статическое давление, Па
        double p0;         // полное давление, Па
        double rho;        // плотность, кг/м^3
        double rho0;       // плотность торможения, кг/м^3
        double T;          // статическая температура, K
        double T0;         // температура торможения, K
        double v;          // скорость, м/с
        double a;          // скорость звука, м/с
        double M;          // число Маха
        double Cp;         // теплоёмкость, Дж/(кг*K)
        double mu;         // динамическая вязкость, Па*с
        double lambda;     // теплопроводность, Вт/(м*K)
        double mu_mach;    // угол Маха, градусы
        double theta_shock;// угол скачка уплотнения, градусы (0 если нет)
        double theta_flow; // угол направления потока, градусы
    };

    // T — температура, K
    double Cp(double T);
    double Viscosity(double T);
    double ThermalConductivity(double T);
    double SpeedOfSound(double T);

    // M — число Маха; p, rho, T — статические параметры
    double StagnationPressure(double p, double M);
    double StagnationDensity(double rho, double M);
    double StagnationTemperature(double T, double M);

    // M — число Маха
    double MachAngleDeg(double M);

    // Заполнить Cp, mu, lambda, a, mu_mach по текущей T
    void FillTransport(GasState& s);

    // Заполнить p0, rho0, T0 по текущим p, rho, T, M
    void FillStagnation(GasState& s);

    // Собрать полное состояние газа по p, T, M и углам
    GasState FromMach(double p, double T, double M,
        double theta_shock_deg = 0.0,
        double theta_flow_deg = 0.0);
}