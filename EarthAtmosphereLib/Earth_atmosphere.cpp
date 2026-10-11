#include "Earth_atmosphere.h"
#include <cmath>

namespace Earth_atmosphere {

    double GeometricToGeopotential(double h) {
        return (R_EARTH * h) / (R_EARTH + h);
    }

    double CalculateGravity(double h) {
        double ratio = R_EARTH / (R_EARTH + h);
        return G_SEA_LEVEL * ratio * ratio;
    }

    double CalculateSpeedOfSound(double T) {
        return std::sqrt(GAMMA * R_AIR * T);
    }

    double CalculateTemperature(double H) {
        if (H <= 11000.0) {
            return T_SEA_LEVEL - 0.0065 * H;
        }
        else if (H <= 20000.0) {
            return 216.65;
        }
        else if (H <= 32000.0) {
            return 216.65 + 0.0010 * (H - 20000.0);
        }
        else if (H <= 47000.0) {
            return 228.65 + 0.0028 * (H - 32000.0);
        }
        else if (H <= 51000.0) {
            return 270.65;
        }
        else if (H <= 71000.0) {
            return 270.65 - 0.0028 * (H - 51000.0);
        }
        else if (H <= 85000.0) {
            return 214.65 - 0.0020 * (H - 71000.0);
        }
        else {
            return 186.65;
        }
    }

    double CalculatePressure(double H, double T) {
        if (H <= 11000.0) {
            double beta = -0.0065;
            return P_SEA_LEVEL * std::pow(T_SEA_LEVEL / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else if (H <= 20000.0) {
            double p11 = CalculatePressure(11000.0, 216.65);
            return p11 * std::exp(-G_SEA_LEVEL * (H - 11000.0) / (R_AIR * 216.65));
        }
        else if (H <= 32000.0) {
            double p20 = CalculatePressure(20000.0, 216.65);
            double beta = 0.0010;
            return p20 * std::pow(216.65 / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else if (H <= 47000.0) {
            double p32 = CalculatePressure(32000.0, 228.65);
            double beta = 0.0028;
            return p32 * std::pow(228.65 / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else if (H <= 51000.0) {
            double p47 = CalculatePressure(47000.0, 270.65);
            return p47 * std::exp(-G_SEA_LEVEL * (H - 47000.0) / (R_AIR * 270.65));
        }
        else if (H <= 71000.0) {
            double p51 = CalculatePressure(51000.0, 270.65);
            double beta = -0.0028;
            return p51 * std::pow(270.65 / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else if (H <= 85000.0) {
            double p71 = CalculatePressure(71000.0, 214.65);
            double beta = -0.0020;
            return p71 * std::pow(214.65 / T, G_SEA_LEVEL / (R_AIR * beta));
        }
        else {
            double p85 = CalculatePressure(85000.0, 186.65);
            return p85 * std::exp(-G_SEA_LEVEL * (H - 85000.0) / (R_AIR * 186.65));
        }
    }

    double CalculateDensity(double p, double T) {
        return p / (R_AIR * T);
    }

    AtmosphereParams Calculate(double h) {
        AtmosphereParams A;
        A.h = h;
        A.H = GeometricToGeopotential(h);
        A.T = CalculateTemperature(A.H);
        A.p = CalculatePressure(A.H, A.T);
        A.rho = CalculateDensity(A.p, A.T);
        A.a = CalculateSpeedOfSound(A.T);
        A.g = CalculateGravity(h);
        return A;
    }

    std::vector<AtmosphereParams> CalculateRange(double h_min, double h_max, double step) {
        std::vector<AtmosphereParams> result;
        for (double h = h_min; h <= h_max; h += step)
            result.push_back(Calculate(h));
        return result;
    }
}