#include "LSMPSA.hpp"

std::vector<std::vector<double>> LSMPSA::Px(Particle p_j, Particle p_i){
    std::vector<std::vector<double>> P(5, std::vector<double>(1,0));
    int ctr=0;
    for (int alp=1; alp<3; alp++){
        for (int xi=alp; xi>=0; xi--){
            int yi = alp - xi;
            P[ctr][0] = pow(p_j.getX()-p_i.getX(),xi) * pow(p_j.getY()-p_i.getY(),yi)/pow(rs,alp);
            ctr++;
        }
    }
    return P;
}

LSMPSA::LSMPSA(std::vector<Particle>& p_n, std::vector<double>& w, Particle p, double Re){
    rs = 0.8*Re;
    // std::cout << "#Log: Performing LSMPS::Number of neighbor: " << p_n.size() << std::endl << std::endl;
    // Calculating Hrs
    int ctr=0;
    for (int alp=1; alp<3; alp++){
        for (int xi=alp; xi>-1; xi--){
            int yi = alp - xi;
            // Forcing factorial
            int X = (xi==0) ? 1 : xi;
            int Y = (yi==0) ? 1 : yi;
            
            Hrs[ctr][ctr] = Y*X/pow(rs,alp);
            ctr++;
        }
    }
    
    // Summation calculation of each neighboring particle
    for (int i=0; i<p_n.size(); i++){
        std::vector<std::vector<double>> PX = Px(p_n[i], p);
        std::vector<std::vector<double>> PXt = linalg::transpose(PX);
        std::vector<std::vector<double>> Pm = linalg::multiMat(PX, PXt);   // Must be checked
        
        // Calculation for M
        for (int r=0; r<5; r++){
            for (int c=0; c<5; c++){
                M[r][c] += w[i] * Pm[r][c];
            }
        };
        
        // Calculation for b
        for (int r=0; r<5; r++){
            b[r] += w[i] * PX[r][0] * (p_n[i].getF() - p.getF());
        };
    }

    std::vector<double> temporary = linalg::GJE(M, b);
    this->Df_x = linalg::multiMat(Hrs, temporary);
}