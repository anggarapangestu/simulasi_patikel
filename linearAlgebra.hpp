#ifndef LINALG_H
#define LINALG_H
#include <vector>
#include <cmath>
#include <iostream>

namespace linalg{
    // Transpose :: Tested
    template <typename T>   
    std::vector<std::vector<T>> transpose(std::vector<std::vector<T>>& A){
        int row = A.size();
        int col = A[0].size();

        std::vector<std::vector<T>> A_t(col, std::vector<T>(row,0));

        // Calculating the value
        for (int i = 0; i < row; i++){
            for (int j = 0; j < col; j++){
                A_t[j][i] = A[i][j];
            }
        }
        return A_t;
    }

    // Print matrix 1D :: Tested
    template <typename T>
    void printMatrix(std::vector<T>& A){
        for (int i=0; i<A.size(); i++){
            std::cout << "[" << A[i] << "]" << std::endl;
        }
        std::cout << std::endl;
    }

    // Print matrix 2D :: Tested
    template <typename T>
    void printMatrix(std::vector<std::vector<T>>& A){
        for (int i=0; i<A.size(); i++){
            std::cout << "[";
            for (int j=0; j<A[i].size(); j++){
                if (j!=0){std::cout << ", ";}
                std::cout << A[i][j];
            }
            std::cout << "]" << std::endl;
        }
        std::cout << std::endl;
    }

    // Multiplication of A x B :: Tested
    template <typename T>
    std::vector<std::vector<T>> multiMat(std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& B){
        std::vector<std::vector<T>> multi;
        int r1 = A.size();
        int c1 = A[0].size();
        int r2 = B.size();
        int c2 = B[0].size();
        
        if (c1 != r2){
            std::cout << "Math Error: The matrix size violate !!" << std::endl;
            return multi;
        }
        multi.resize(r1);
        for (int i=0; i<r1; i++){multi[i].resize(c2,0);}

        // Calculating the value
        for (int i = 0; i < r1; i++){
            for (int j = 0; j < c2; j++){
                for (int k = 0; k < c1; k++){
                    multi[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return multi;
    }

    template <typename T>
    std::vector<T> multiMat(std::vector<std::vector<T>>& A, std::vector<T>& B){
        std::vector<T> multi;
        int r1 = A.size();
        int c1 = A[0].size();
        int r2 = B.size();
        
        if (c1 != r2){
            std::cout << "Math Error: The matrix size violate !!" << std::endl;
            return multi;
        }
        
        multi.resize(r1);
        // Calculating the value
        for (int i = 0; i < r1; i++){
            for (int k = 0; k < c1; k++){
                multi[i] += A[i][k] * B[k];
            }
        }
        return multi;
    }

    // Determinant by Row Expansion :: Tested
    template <typename T>
    T det(std::vector<std::vector<T>>& A){
        if (A.size() != A[0].size()){
            std::cout << "Math Error: Matrix determinant is violated, the matrix is not square!" << std::endl;
            return 0.0;
        }

        int s = A.size();
        T det = 0;
        T val;

        // Set the first rows as the multiplication row
        for (int i=0; i<s; i++){
            val = pow(-1,i) * A[0][i];
            
            // Get the determinant
            if (s-1 == 1){
                val *= A[1][1-i];   // value of the cross element
            }else{
                // Define the new component determinant matrix
                std::vector<std::vector<T>> detMat(s-1, std::vector<T>(s-1,0));
                
                // Assign each component
                for (int j = 0; j < s-1; j++){
                    for (int k = 0; k < s-1; k++){
                        if (k >= i){
                            detMat[j][k] = A[j+1][k+1];
                        }else{
                            detMat[j][k] = A[j+1][k];
                        }
                    }
                }
                val *= linalg::det(detMat);
            }
            // Add the determinant
            det += val;
        }
        return det;
    }

    // Lin Alg solution by Gauss-Jordan elemination method :: Tested
    template <typename T>
    std::vector<T> GJE(std::vector<std::vector<T>> A, std::vector<T> b){
        if (A.size() != A[0].size()){
            std::cout << "Math Error: The matrix is not square!" << std::endl;
            return b;
        }

        int s = A.size();

        // Check the singularity of matrix
        double det = linalg::det(A);
        if (std::abs(det) < pow(10,-10)){
            std::cout << "Math Error: The matrix is singular" << std::endl;
            return b;
        }
        
        double tempDat;     // Storing temporary data for change row
        double mul;         // Storing the multiplication factor for determinant
        int counter;        // Signify the position of row interchange
        double red;         // Storing the reducting factor for determinant

        for(int i=0; i<s; i++){
            // Loop check the pivot, no zero allowed
            counter = i+1;
            while (A[i][i] == 0.0){
                // Interchange all row element
                for (int k=0; k<s; k++){        // Change from 1 to s
                    tempDat = A[i][k];          // Store the pivot row element value
                    A[i][k] = A[counter][k];    // Store the pivot row element value
                    A[counter][k] = tempDat;
                }
                
                tempDat = b[i];         // Store the pivot row element value
                b[i] = b[counter];      // Store the pivot row element value
                b[counter] = tempDat;

                // Update the counter
                counter++;
            }
            
            // Multiply all coresponding row by the invers determinant factor
            mul = A[i][i];
            for (int k=i; k<s; k++){
                A[i][k] /= mul;    // Normalize all pivot row element
            }
            b[i] /= mul;    // Normalize all pivot row element

            // Remember the pivot index is 'i'
            // Reduce the other rows
            for(int j=0; j<s; j++){         // Start from the one next pivot row
                if (j == i){
                    continue;
                }
                if (A[j][i] != 0){     // Check whether the value is already zero or not
                    red = A[j][i];
                    for (int k=0; k<s; k++){
                        // Reduce all element
                        A[j][k] -= A[i][k]*red;
                    }
                    b[j] -= b[i]*red;
                }
            }
        }
        return b;
    }

    // Inverse by Gauss-Jordan elemination method :: Tested
    template <typename T>
    std::vector<std::vector<T>> inverse(std::vector<std::vector<T>> A){
        if (A.size() != A[0].size()){
            std::cout << "Math Error: The matrix is not square!" << std::endl;
            return A;
        }

        int s = A.size();

        // Initiation of the inverse matrix //
        std::vector<std::vector<T>> inv(s, std::vector<T>(s,0));
        for (int i=0; i<s; i++){inv[i][i]=1;}

        // Check the singularity of matrix
        double det = linalg::det(A);
        if (std::abs(det) < pow(10,-10)){
            std::cout << "Math Error: The matrix is singular" << std::endl;
            return A;
        }
        
        double tempDat;     // Storing temporary data for change row
        double mul;         // Storing the multiplication factor for determinant
        int counter;        // Signify the position of row interchange
        double red;         // Storing the reducting factor for determinant

        for(int i=0; i<s; i++){
            // Loop check the pivot, no zero allowed
            counter = i+1;
            while (A[i][i] == 0.0){
                // Interchange all row element
                for (int k=0; k<s; k++){        // Change from 1 to s
                    tempDat = A[i][k];          // Store the pivot row element value
                    A[i][k] = A[counter][k];    // Store the pivot row element value
                    A[counter][k] = tempDat;

                    tempDat = inv[i][k];          // Store the pivot row element value
                    inv[i][k] = inv[counter][k];    // Store the pivot row element value
                    inv[counter][k] = tempDat;
                }

                // Update the counter
                counter++;
            }
            
            // Multiply all coresponding row by the invers determinant factor
            mul = A[i][i];
            for (int k=0; k<s; k++){
                A[i][k] /= mul;    // Normalize all pivot row element
                inv[i][k] /= mul;    // Normalize all pivot row element
            }

            // Remember the pivot index is 'i'
            // Reduce the other rows
            for(int j=0; j<s; j++){         // Start from the one next pivot row
                if (j == i){
                    continue;
                }
                if (A[j][i] != 0){     // Check whether the value is already zero or not
                    red = A[j][i];
                    for (int k=0; k<s; k++){
                        // Reduce all element
                        A[j][k] -= A[i][k]*red;
                        inv[j][k] -= inv[i][k]*red;
                    }
                }
            }
        }
        return inv;
    }

}
#endif