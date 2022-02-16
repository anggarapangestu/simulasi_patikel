#include <iostream>
#include <vector>
#include <fstream>
#include "LSMPSA.hpp"
#include "spatialHash.hpp"

// Test Function f(x,y) = x^3 - y^3 + x^2 y^2
double fun_f(double x, double y){
    return (pow(x,3) - pow(y,3) + pow(x,2)*pow(y,2));
}
// Analytical solution
double fun_dfx(double x, double y){
    return (3*pow(x,2) + 2*x*pow(y,2));
}
double fun_dfy(double x, double y){
    return ( -3*pow(y,2) + 2*y*pow(x,2));
}
double fun_dfx2(double x, double y){
    return (6*x + 2*pow(y,2));
}
double fun_dfxy(double x, double y){
    return (4*x*y);
}
double fun_dfy2(double x, double y){
    return (-6*y + 2*pow(x,2));
}

// Write file:: (1) Derivative to one substance, (2) all derivatives
void write_File(std::vector<Particle> P, double (&DFA)[][5], double (&DFL)[][5], int index){
    std::cout << "#Log: Writting file !!" << std::endl;
    std::string name1;
    std::ofstream ofs;
    
    std::string caseType;
    switch (index){
        case 0:caseType = "dx";break;
        case 1:caseType = "dy";break;
        case 2:caseType = "dx2";break;
        case 3:caseType = "dxy";break;
        case 4:caseType = "dy2";break;
    }

    name1.append("Hasil_Perhitungan_Df");
    name1.append(caseType);
    name1.append(".csv");
	
    ofs.open(name1.c_str());
	ofs << "x" << "," << "y" << "," << "f(x)" 
        << "," <<"df(x)/dx"<< "," << "df(x)/dx LSMPS" 
        << "," << "difference\n"
    ;
	
    for (int i = 0; i < P.size(); i++){
        ofs << P[i].getX() << "," << P[i].getY() << "," << P[i].getF()
            << "," << DFA[i][index] << "," << DFL[i][index]
            << "," << pow(DFA[i][index] - DFL[i][index], 2)
            << "\n"
        ;
	}
	
    ofs.close();
}

void write_File_All(std::vector<Particle> P, double (&DFA)[][5], double (&DFL)[][5]){
    std::cout << "#Log: Writting file !!" << std::endl;
    
    std::string name1;
    std::ofstream ofs;
    name1.append("Hasil_Perhitungan_LSMPSA");
    name1.append(".csv");
	
    ofs.open(name1.c_str());
	ofs << "x" << "," << "y" << "," << "f(x)" << "," 
        <<"df/dx_A"<< "," << "df/dx_LSMPS" << "," << "Ddx" << ","
        <<"df/dy_A"<< "," << "df/dy_LSMPS" << "," << "Ddy" << ","
        <<"df/dx2_A"<< "," << "df/dx2_LSMPS" << "," << "Ddx2" << ","
        <<"df/dxy_A"<< "," << "df/dxy_LSMPS" << "," << "Ddxy" << ","
        <<"df/dy2_A"<< "," << "df/dy2_LSMPS" << "," << "Ddy2" << "\n"
	;

    for (int i = 0; i < P.size(); i++){
        ofs << P[i].getX() << "," << P[i].getY() << "," << P[i].getF()
            << "," << DFA[i][0] << "," << DFL[i][0] << "," << pow(DFA[i][0] - DFL[i][0], 2)
            << "," << DFA[i][1] << "," << DFL[i][1] << "," << pow(DFA[i][1] - DFL[i][1], 2)
            << "," << DFA[i][2] << "," << DFL[i][2] << "," << pow(DFA[i][2] - DFL[i][2], 2)
            << "," << DFA[i][3] << "," << DFL[i][3] << "," << pow(DFA[i][3] - DFL[i][3], 2)
            << "," << DFA[i][4] << "," << DFL[i][4] << "," << pow(DFA[i][4] - DFL[i][4], 2)
            << "\n";
	}
	
    ofs.close();
}


int main(){
    int Nx;
    int Ny;
    float divider;
    std::cout << "Input the number of particle in x direction: ";
    std::cin >> Nx;
    std::cout << "Input the number of particle in y direction: ";
    std::cin >> Ny;

    std::cout << "Input the divider: ";
    std::cin >> divider;

    int N = Nx*Ny;
    double Re = 0.8;
    double Dfl[N][5] = {};
    double Dfa[N][5] = {};
    double err[N][5] = {};

    std::vector<Particle> par(N);
    
    // Particle generation
    int pos=0;
    double Xx=0;
    
    for (int i=0; i<Nx; i++){
        double Yy=0;
        for (int j=0; j<Ny; j++){
            par[pos].setParticle(Xx, Yy, 0, fun_f(Xx,Yy));
            Dfa[pos][0] = fun_dfx(Xx, Yy);
            Dfa[pos][1] = fun_dfy(Xx, Yy);
            Dfa[pos][2] = fun_dfx2(Xx, Yy);
            Dfa[pos][3] = fun_dfxy(Xx, Yy);
            Dfa[pos][4] = fun_dfy2(Xx, Yy);
            Yy += Re/divider;
            pos++;
        }
        Xx += Re/divider;
    }

    // Spatial Hashing Table
    SpatialHash sh(&par[0], N, Re);
    
    // Calculate each particle    
    for (int i=0; i<N; i++){
        if (i%100 == 0){std::cout << "LSMPS A: Evaluating on particle " << i << std::endl;}
        sh.neighborSearch(&par[0], i);
        std::vector<Particle> p_n = sh.getNeigborPar();
        std::vector<double> W = sh.getNeigborWeight();
        
        // Catch a problem in LSMPSA Code
        LSMPSA lsmpsa(p_n, W, par[i], Re);

        Dfl[i][0] = lsmpsa.getDx();
        Dfl[i][1] = lsmpsa.getDy();
        Dfl[i][2] = lsmpsa.getDx2();
        Dfl[i][3] = lsmpsa.getDxy();
        Dfl[i][4] = lsmpsa.getDy2();

        sh.resetNeighbor();
    }
    
    // write_File(par, Dfa, Dfl, 1);
    write_File_All(par, Dfa, Dfl);
    
    /*/ === Neighbor Search Debugger === // Remove * to debug
    // Check neighborhood
    int index = 8;
    sh.neighborSearch(&par[0],index);
    std::cout << "For particle " << index << ", x= " << par[index].getX() << ", y= " << par[index].getY() << "The neighborhood are below:" << std::endl;
    
    std::vector<Particle> NP = sh.getNeigborPar();
    std::vector<double> weight = sh.getNeigborWeight();

    for (int i=0; i<weight.size(); i++){
        std::cout << " > Particle " << i+1 << " x= " << NP[i].getX() << ", y= " << NP[i].getY() << ", w= " << weight[i] << std::endl;
    }
    std::cout << std::endl;
    // === Neighbor Search Debugger === /*/

    /*/ === Spatial Hash Debugger === // Remove * to debug
    for (int i=0; i<N; i++){
        std::cout << "Particle " << i << std::endl;
        std::cout << " > x=" << par[i].x << std::endl;
        std::cout << " > y=" << par[i].y << std::endl;
        std::cout << " > Xgrid=" << par[i].grid_x << std::endl;
        std::cout << " > Ygrid=" << par[i].grid_y << std::endl << std::endl;
    }
    
    for (int i=0; i<sh.getn_x(); i++){
        for (int j=0; j<sh.getn_y(); j++){
            std::cout << "All Particle in grid " << i << ", " << j << std::endl;
            int counter=1;
            for (auto& val:sh.getParticleIndexes(i,j)){
                std::cout << " > index " << counter << " : " << val << std::endl;
                counter++;
            }
        }
    }
    // === Spatial Hash Debugger === /*/
}