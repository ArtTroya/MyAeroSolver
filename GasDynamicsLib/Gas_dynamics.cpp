#include "Gas_dynamics.h"
#include "Earth_atmosphere.h"
#include <cmath>

namespace Gas_dynamics {

    double Cp(double T) {
        return Cp_REF * std::pow(T / T_NORM, PHI_CP);
    }

    double Viscosity(double T) {
        return MU_REF * std::pow(T / T_NORM, N_MU);
    }

    double ThermalConductivity(double T) {
        return LAMBDA_REF * std::pow(T / T_LAMBDA_REF, XI_LAMBDA);
    }

    double SpeedOfSound(double T) {
        return std::sqrt(Earth_atmosphere::GAMMA * Earth_atmosphere::R_AIR * T);
    }

    double StagnationPressure(double p, double M) {
        const double k = Earth_atmosphere::GAMMA;
        return p * std::pow(1.0 + 0.5 * (k - 1.0) * M * M, k / (k - 1.0));
    }

    double StagnationDensity(double rho, double M) {
        const double k = Earth_atmosphere::GAMMA;
        return rho * std::pow(1.0 + 0.5 * (k - 1.0) * M * M, 1.0 / (k - 1.0));
    }

    double StagnationTemperature(double T, double M) {
        const double k = Earth_atmosphere::GAMMA;
        return T * (1.0 + 0.5 * (k - 1.0) * M * M);
    }

    double MachAngleDeg(double M) {
        if (M <= 1.0) return 90.0;
        return std::asin(1.0 / M) * 180.0 / PI;
    }

    void FillTransport(GasState& s) {
        s.Cp      = Cp(s.T);
        s.mu      = Viscosity(s.T);
        s.lambda  = ThermalConductivity(s.T);
        s.a       = SpeedOfSound(s.T);
        s.mu_mach = (s.M > 1.0) ? std::asin(1.0 / s.M) * 180.0 / PI : 90.0;
    }

    void FillStagnation(GasState& s) {
        s.p0   = StagnationPressure(s.p, s.M);
        s.rho0 = StagnationDensity(s.rho, s.M);
        s.T0   = StagnationTemperature(s.T, s.M);
    }

    GasState FromMach(double p, double T, double M,
                      double theta_shock_deg,
                      double theta_flow_deg) {
        GasState s{};
        s.p = p;
        s.T = T;
        s.M = M;
        s.rho = p / (Earth_atmosphere::R_AIR * T);
        s.a   = SpeedOfSound(T);
        s.v   = M * s.a;
        s.theta_shock = theta_shock_deg;
        s.theta_flow  = theta_flow_deg;
        FillTransport(s);
        FillStagnation(s);
        return s;
    }
}