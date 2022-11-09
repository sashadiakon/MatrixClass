#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <vector>
#include <cmath>
#include<cstdlib>

using namespace std;

class Matrix {
    public:
    
        Matrix(int, int);
        Matrix(int);
        Matrix(const Matrix&);
        Matrix(vector <float>);

        vector<float>& operator[](int index);

        Matrix& operator+=(const Matrix&);
        Matrix& operator+=(double);
        Matrix& operator-=(const Matrix&);
        Matrix& operator-=(double);
        Matrix& operator*=(const Matrix&);
        Matrix& operator*=(double);
        Matrix& operator/=(double);
        Matrix operator~();
        void operator=(const Matrix&);
        bool operator==(const Matrix&);
        
        friend ostream& operator<<(ostream&, const Matrix&);
        friend istream& operator>>(istream&, Matrix&);

        void print();
        int getRows();
        bool is_square();
        bool is_diag();
        void Identity();
        Matrix Inverse();
        Matrix getCofactor(int n, int p, int q);
        float Minor(int n, int i, int j);
        float Det(int n);
        vector<float> GaussianElimination();
        pair<int, int> nondiag();
        Matrix Triangular();
        vector<float> JacobiAlg(float eps);
        void RandomSymetric(float min, float max);

    private:
        int rows, cols;
        vector<vector<float>> values;
        

};

Matrix operator+(const Matrix&, const Matrix&);
Matrix operator+(const Matrix&, double);
Matrix operator+(double, const Matrix&);
Matrix operator-(const Matrix&, const Matrix&);
Matrix operator-(const Matrix&, double);
Matrix operator-(double, const Matrix&);
Matrix operator*(const Matrix&, const Matrix&);
Matrix operator*(const Matrix&, double);
Matrix operator*(double, const Matrix&);
Matrix operator/(const Matrix&, double);



#endif
