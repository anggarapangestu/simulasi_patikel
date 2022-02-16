#ifndef LSMPSA_H
#define LSMPSA_H
#include "linearAlgebra.hpp"
#include <vector>
#include <cmath>
#include <iostream>

#ifndef PARTICLE
#define PARTICLE
#include "particle.hpp"
#endif

class LSMPSA{
    private:
        double rs = 1;
        std::vector<double> Df_x{std::vector<double>(5,0)};
        std::vector<std::vector<double>> Hrs{std::vector<std::vector<double>>(5, std::vector<double>(5,0))};
        std::vector<std::vector<double>> M{std::vector<std::vector<double>>(5, std::vector<double>(5,0))};
        std::vector<double> b{std::vector<double>(5, 0)};

        std::vector<std::vector<double>> Px(Particle p_j, Particle p_i);
    public:
        // Calculation of LSMPS
        LSMPSA(std::vector<Particle>& p_n, std::vector<double>& w, Particle p, double Re);
        LSMPSA(){};
        ~LSMPSA(){};

        // Get the LSMPS result
        double getDx(){return Df_x[0];};
        double getDy(){return Df_x[1];};
        double getDx2(){return Df_x[2];};
        double getDxy(){return Df_x[3];};
        double getDy2(){return Df_x[4];};

        //void resetLSMPSA(){Df_x.clear();weight.clear();};  // Reset the neighbor and weight
        
};

#endif