#include "spatialHash.hpp"

// ==== Spatial Hash Function ==== //
// resize the object Index
SpatialHash::SpatialHash(Particle* p, int n, double r){
    re = r;
    s = re*1.05;
    int ymax=0, ymin=0, xmax=0, xmin=0;

    // Find the index of max and min of x and y --> Algorithm, may be optimized
    for (int i=0; i<n; i++){
        xmin = ((*(p+xmin)).getX() > (*(p+i)).getX()) ? i : xmin;
        xmax = ((*(p+xmax)).getX() < (*(p+i)).getX()) ? i : xmax;
        ymin = ((*(p+ymin)).getY() > (*(p+i)).getY()) ? i : ymin;
        ymax = ((*(p+ymax)).getY() < (*(p+i)).getY()) ? i : ymax;
    }

    // Find the grid length in x and y direction
    // The expand factor is set 1.05, this value may be changed
    double lenGridX = ((*(p+xmax)).getX() - (*(p+xmin)).getX())*1.05;
    double lenGridY = ((*(p+ymax)).getY() - (*(p+ymin)).getY())*1.05;

    // Find the hash Table origin
    double X0 = (*(p+xmin)).getX() - lenGridX*(0.025/1.05);
    double Y0 = (*(p+ymin)).getY() - lenGridY*(0.025/1.05);

    // Calculate the number of grid
    // The + 1 term is to ceil the value
    n_x = lenGridX/s + 1;
    n_y = lenGridY/s + 1;
    
    // Resize the grid
    hashTable.resize(n_x);
    for (int i=0; i<n_x; i++){hashTable[i].resize(n_y);};

    // Assign each particle grid index
    for (int i=0; i<n; i++){
        // i is the particle index
        int posX = ((*(p+i)).getX() - X0)/s;
        int posY = ((*(p+i)).getY() - Y0)/s;
        
        // Add the index of particle i in the hashTable grid
        hashTable[posX][posY].push_back(i);
        
        // objectIndex of particle is assigned
        (*(p+i)).setGrid(posX, posY, 0);
    }
}

// resize the object Index
void SpatialHash::setValue(Particle* p, int n, double r){
    re = r;
    s = re*1.05;
    int ymax=0, ymin=0, xmax=0, xmin=0;

    // Find the index of max and min of x and y --> Algorithm, may be optimized
    for (int i=0; i<n; i++){
        xmin = ((*(p+xmin)).getX() > (*(p+i)).getX()) ? i : xmin;
        xmax = ((*(p+xmax)).getX() < (*(p+i)).getX()) ? i : xmax;
        ymin = ((*(p+ymin)).getY() > (*(p+i)).getY()) ? i : ymin;
        ymax = ((*(p+ymax)).getY() < (*(p+i)).getY()) ? i : ymax;
    }

    // Find the grid length in x and y direction
    // The expand factor is set 1.05, this value may be changed
    double lenGridX = ((*(p+xmax)).getX() - (*(p+xmin)).getX())*1.05;
    double lenGridY = ((*(p+ymax)).getY() - (*(p+ymin)).getY())*1.05;

    // Find the hash Table origin
    double X0 = (*(p+xmin)).getX() - lenGridX*(0.025/1.05);
    double Y0 = (*(p+ymin)).getY() - lenGridY*(0.025/1.05);

    // Calculate the number of grid
    // The + 1 term is to ceil the value
    n_x = lenGridX/s + 1;
    n_y = lenGridY/s + 1;
    
    // Resize the grid
    hashTable.resize(n_x);
    for (int i=0; i<n_x; i++){hashTable[i].resize(n_y);};

    // Assign each particle grid index
    for (int i=0; i<n; i++){
        // i is the particle index
        int posX = ((*(p+i)).getX() - X0)/s;
        int posY = ((*(p+i)).getY() - Y0)/s;
        
        // Add the index of particle i in the hashTable grid
        hashTable[posX][posY].push_back(i);
        
        // objectIndex of particle is assigned
        (*(p+i)).setGrid(posX, posY, 0);
    }
}

// Give all index of the particles in grid x, y
// Return a vector by value
std::vector<int> SpatialHash::getParticleIndexes(int xpos, int ypos){
    return hashTable[xpos][ypos];
}

// ==== Neighbor Search Function ==== //
// Private weight function 3D applied
double SpatialHash::weightCalc (Particle p_j, Particle p_i){
    double W = 0;   // Initial weight
    double R = sqrt(pow(p_j.getX() - p_i.getX(),2) 
             + pow(p_j.getY() - p_i.getY(),2)
             + pow(p_j.getZ() - p_i.getZ(),2));   // Distance between particle
    // Calculating squared weight
    W = (R < re) ? pow(1 - R/re, 2) : 0;
    return W;
}

// Calculating weight function
void SpatialHash::neighborSearch(Particle* p, int id){
    // Perform a neighbor finding + Weight calculation
    
    // Get the corresponding particle Hash grid
    int gridX = (*(p+id)).getGridX();
    int gridY = (*(p+id)).getGridY();
    
    std::vector<int> neighborIndexes;
    // Store all particle candidate index
    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            int gridPosX = gridX + i-1;
            int gridPosY = gridY + j-1;
            
            // Confined between the grid domain
            if (gridPosX<0 || gridPosX >= n_x){continue;}
            if (gridPosY<0 || gridPosY >= n_y){continue;}

            // Assign the neighboring particle index
            for (auto& val:this->getParticleIndexes(gridPosX, gridPosY)){
                if (val==id){continue;}
                neighborIndexes.push_back(val);
            }            
        }
    }
    
    // Perform the weight function
    // ================= CODE ================= //
    // Check the weight
    // Store the particle and the correspond weight if fulfill the neighbor condition
    for (auto& val:neighborIndexes){
        double w = this->weightCalc(*(p+val), *(p+id));
        if (w != 0){
            neighborPar.push_back(*(p+val));
            weight.push_back(w);
        }
    }
}

// Give the neighboring particle
std::vector<Particle> SpatialHash::getNeigborPar(){
    return neighborPar;
}

// Give the weight of neighboring particle
std::vector<double> SpatialHash::getNeigborWeight(){
    return weight;
}