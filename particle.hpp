class Particle{
    private:
        double x, y, z, f;
        double dfx, dfy, dfx2, dfxy, dfy2;  // Analytical
        double Dfx, Dfy, Dfx2, Dfxy, Dfy2;  // LSMPS calculation
        int grid_x, grid_y, grid_z;
    public:
        // Particle position, value, and grid position definition
        void setParticle(double X, double Y, double Z, double F){
            x = X;
            y = Y;
            z = Z;
            f = F;
        }
        void setGrid(int gridX, int gridY, int gridZ){
            grid_x = gridX;
            grid_y = gridY;
            grid_z = gridZ;
        }
        
        // Particle position and value getter
        double getX(){
            return x;
        }
        double getY(){
            return y;
        }
        double getZ(){
            return z;
        }
        double getF(){
            return f;
        }

        // Particle grid position getter
        int getGridX(){
            return grid_x;
        }
        int getGridY(){
            return grid_y;
        }
        int getGridZ(){
            return grid_z;
        }
        
};