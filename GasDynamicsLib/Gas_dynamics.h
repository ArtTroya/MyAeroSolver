#pragma once

namespace Gas_dynamics {

    constexpr double PI = 3.14159265358979323846;

    // Опорные значения свойств воздуха
    constexpr double Cp_REF = 1004.7;       
    constexpr double MU_REF = 1.79e-5;      
    constexpr double LAMBDA_REF = 0.0232;   
    constexpr double T_NORM = 288.15;      
    constexpr double T_LAMBDA_REF = 261.0;       

    // Показатели степенных зависимостей
    constexpr double PHI_CP = 0.1;
    constexpr double N_MU = 0.76;
    constexpr double XI_LAMBDA = 0.85;

    struct GasState {
        double p;
        double p0;
        double rho;
        double rho0;
        double T;
        double T0;
        double v;
        double a;
        double M;
        double Cp;
        double mu;
        double lambda;
        double mu_mach;
        double theta_shock;
        double theta_flow;
    };

    double Cp(double T);
    double Viscosity(double T);
    double ThermalConductivity(double T);
    double SpeedOfSound(double T);

    double StagnationPressure(double p, double M);
    double StagnationDensity(double rho, double M);
    double StagnationTemperature(double T, double M);

    double MachAngleDeg(double M);

    void FillTransport(GasState& s);
    void FillStagnation(GasState& s);

    GasState FromMach(double p, double T, double M,
        double theta_shock_deg = 0.0,
        double theta_flow_deg = 0.0);
}