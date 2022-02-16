#ifndef SH_H
#define SH_H
#include <vector>
// #include <iostream>
#include <cmath>

#ifndef PARTICLE
#define PARTICLE
#include "particle.hpp"
#endif

// ================== Special Notes: ==================
// Spatial Hash improve the running time by not evaluating all particle in the domain for neighbor search
// This class store the index data of particle in the domain at one time
// This class is used only when finding neighboring particle

// This class have a function to give the grid position of the particle
// This class have a function to give particle index inside a grid

class SpatialHash{
    private:
        // The spatial hash algortihm
        int n_x;    // number of x-directed grid
        int n_y;    // number of y-directed grid
        double re;  // radius of influence
        double s;   // size of the spatial tabel box (square shape)
        std::vector<std::vector<std::vector<int>>> hashTable;   // Store all particle index at corresponding grid

        // Neighbor search algortihm + Weight calculation
        std::vector<Particle> neighborPar;   // store the neighboring particle data for LSMPS A Calculation
        std::vector<double> weight;     // store the weight of neighboring particle
        double weightCalc (Particle p1, Particle p2);   // Calculating weight
    public:
        // The spatial hash algortihm
        SpatialHash(Particle* P, int N, double Re);    // resize the object Index
        void setValue(Particle* P, int N, double Re);  // resize the object Index
        SpatialHash(){};    // Default constructor
        ~SpatialHash(){};   // Default deconstructor
        std::vector<int> getParticleIndexes(int posX, int posY);  // Give all index of the particles in grid x, y
        int getn_x(){return n_x;};
        int getn_y(){return n_y;};

        // Neighbor search algortihm + Weight calculation
        void neighborSearch(Particle* P, int ind);    // Perform a neighbor finding + Weight calculation
        std::vector<Particle> getNeigborPar();      // Return the vector of neighbor particle index
        std::vector<double> getNeigborWeight();     // Return the vector of weight
        void resetNeighbor(){neighborPar.clear();weight.clear();};  // Reset the neighbor and weight
};

#endif
// The neighbor search: storing the neighboring particle and calculating the weight of neighboring particle