#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <vector>

using namespace std;

class Matrix {
    public:
        Matrix(int, int);
        Matrix(int);
        Matrix(const Matrix&);

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


    private:
        vector<vector<float>> values;
        int rows, cols;

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
