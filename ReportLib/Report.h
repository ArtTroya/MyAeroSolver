#pragma once
#include "Earth_atmosphere.h"
#include "Flow_regions.h"
#include "Aerodynamics.h"
#include "Boundary_layer.h"
#include <array>
#include <string>

namespace Report {

    struct CaseInput {
        int    N;
        double h;
        double M_inf;
        double alpha_deg;
        double beta_deg;
        double T_wall;
        double c_m;
        double b_m;
        double L_face;
    };

    struct Summary {
        double q_inf;
        double Cx, Cy, mz;
        double Cxa, Cya, K, Cd;
        double Cf_sigma;
        double X_friction;
        double Cx_sum, Cxa_sum, Cya_sum;
        std::array<double, 4> Cf_avg;
        std::array<double, 4> F_fr;
    };

    void WriteReportCsv(const std::string& path,
        const CaseInput& in,
        const Earth_atmosphere::AtmosphereParams& atm,
        const Flow_regions::FlowField& F,
        const Aerodynamics::AeroCoefficients& aero,
        const Boundary_layer::SurfaceBL& upper,
        const Boundary_layer::SurfaceBL& lower,
        const Summary& S);

    void WriteGraphsCsv(const std::string& path,
        const Boundary_layer::SurfaceBL& upper,
        const Boundary_layer::SurfaceBL& lower,
        double L_face);
}