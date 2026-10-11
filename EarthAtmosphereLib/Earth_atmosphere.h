#pragma once
#include <vector>

namespace Earth_atmosphere {

    constexpr double R_EARTH = 6356767.0;
    constexpr double G_SEA_LEVEL = 9.80665;
    constexpr double P_SEA_LEVEL = 101325.0;
    constexpr double T_SEA_LEVEL = 288.15;
    constexpr double R_AIR = 287.05287;
    constexpr double GAMMA = 1.4;

    struct AtmosphereParams {
        double h;       // геометрическая высота, м
        double H;       // геопотенциальная высота, м
        double T;       // температура, K
        double p;       // давление, Па
        double rho;     // плотность, кг/м^3
        double a;       // скорость звука, м/с
        double g;       // ускорение свободного падения, м/с^2
    };

    // h  — геометрическая высота, м
    AtmosphereParams Calculate(double h);

    // h  — геометрическая высота, м
    double GeometricToGeopotential(double h);

    // h  — геометрическая высота, м
    double CalculateGravity(double h);

    // H  — геопотенциальная высота, м
    double CalculateTemperature(double H);

    // H — геопотенциальная высота, м; T — температура, K
    double CalculatePressure(double H, double T);

    // p — давление, Па; T — температура, K
    double CalculateDensity(double p, double T);

    // T — температура, K
    double CalculateSpeedOfSound(double T);

    // h_min, h_max — границы, м; step — шаг, м
    std::vector<AtmosphereParams> CalculateRange(double h_min, double h_max, double step);
}