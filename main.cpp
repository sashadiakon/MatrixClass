#include "main.h"


Matrix::Matrix(int r, int c){
    rows = r;
    cols = c;
    for (int i = 0; i< rows; i++){
        values.push_back(vector<float> {});
        for (int j = 0; j< cols; j++){
            values[i].push_back(0);
        }
    }
}


Matrix::Matrix(const Matrix& m){
    values = m.values;
    rows = m.rows;
    cols = m.cols;
}


Matrix Matrix::operator~(){
    // Transpose function
    Matrix ret(cols, rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            ret.values[j][i] = values[i][j];
        }
    }
    return ret;
}


Matrix& Matrix::operator+=(double real){
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            values[i][j] += real;
        }
    }
    return *this;
}


Matrix& Matrix::operator+=(const Matrix& m){
    if (cols == m.cols && rows == m.rows){
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                values[i][j] += m.values[i][j];
            }
        }
    }
    else cout<<"\nMatrices have different shapes\n";
    return *this;
}



Matrix& Matrix::operator-=(double real){
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            values[i][j] -= real;
        }
    }
    return *this;
}


Matrix& Matrix::operator-=(const Matrix& m){
    if (cols == m.cols && rows == m.rows){
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                values[i][j] -= m.values[i][j];
            }
        }
    }
    else cout<<"\nMatrices have different shapes\n";
    return *this;
}


Matrix& Matrix::operator*=(double real){
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            values[i][j] *= real;
        }
    }
    return *this;
}


Matrix& Matrix::operator*=(const Matrix& m){
    Matrix temp(rows, m.cols);
    if (cols == m.rows){
        for (int i = 0; i < temp.rows; i++) {
            for (int j = 0; j < temp.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    temp.values[i][j] += (values[i][k] * m.values[k][j]);
                }
            }
        }
    }
    else cout<<"Not appropriate sizes for matrix multiplication\n";
    *this = temp;
    return *this;
}


Matrix& Matrix::operator/=(double real){
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            values[i][j] /= real;
            
        }
    }
    return *this;
}


vector<float> &Matrix::operator[](int index){
    return values[index];
}


void Matrix::operator=(const Matrix& m){
    values = m.values;
    rows = m.rows;
    cols = m.cols;
}


bool Matrix::operator==(const Matrix& m){
    if (values == m.values){
        return true;
    }
    return false;
}


Matrix operator+(const Matrix& m1, const Matrix& m2)
{
    Matrix temp(m1);
    return (temp += m2);
}

Matrix operator+(const Matrix& m1, double real)
{
    Matrix temp(m1);
    return (temp += real);
}

Matrix operator+(double real, const Matrix& m1)
{
    return (m1 + real);
}


Matrix operator-(const Matrix& m1, const Matrix& m2)
{
    Matrix temp(m1);
    return (temp -= m2);
}

Matrix operator-(const Matrix& m1, double real)
{
    Matrix temp(m1);
    return (temp -= real);
}


Matrix operator*(const Matrix& m1, const Matrix& m2)
{
    Matrix temp(m1);
    return (temp *= m2);
}

Matrix operator*(const Matrix& m, double num)
{
    Matrix temp(m);
    return (temp *= num);
}

Matrix operator*(double num, const Matrix& m)
{
    return (m * num);
}

Matrix operator/(const Matrix& m, double num)
{
    Matrix temp(m);
    return (temp /= num);
}

ostream& operator<<(ostream& os, const Matrix& m)
{
    for (int i = 0; i < m.rows; i++) {
        os << m.values[i][0];
        for (int j = 1; j < m.cols; ++j) {
            os << " " << m.values[i][j];
        }
        os << endl;
    }
    return os;
}


istream& operator>>(istream& is, Matrix& m)
{
    for (int i = 0; i < m.rows; i++) {
        for (int j = 0; j < m.cols; ++j) {
            is >> m.values[i][j];
        }
    }
    return is;
}


void Matrix::print(){
    for (auto row : values){
        cout<<endl;
        for (float el: row){
            cout<<el<<" ";
        }
    }
    cout<<endl;
}

bool Matrix::is_square(){
    return rows == cols;
}

Matrix Matrix::Triangular(){
    Matrix tr_m = *this;
    for (int b = 0; b<rows-1; b++){
        for (int i = b+1; i<rows; i++){
            if (tr_m[b][b] == 0){
                cout<<"Mathematical Error";
                exit(0);
            }
            float f = tr_m[i][b]/tr_m[b][b];
            for (int j = 0; j<cols; j++){
                // cout<<value[i][j]<<"  "<<values[]<<"\n";
                tr_m[i][j]-=tr_m[b][j]*f;
            }
        }
    }
    return tr_m;
}

float* Matrix::GaussianElimination(){
    // Check size
    
    if (rows+1 != cols){
        cout<<"Not apropriate size, matrix should have size - (n, n+1)\n";
        exit(0);
    }
    float x[rows];
    Matrix tr_m = this->Triangular();
    // Find solutions:
    if (tr_m[rows-1][cols-2] == 0){
        cout<<"Mathematical Error";
        exit(0);
    }

    x[rows-1] = tr_m[rows-1][cols-1]/tr_m[rows-1][cols-2];

    for(int i=rows-2;i>=0;i--)
    {
        x[i] = tr_m[i][cols-1];
        for(int j=i+1;j<cols;j++)
        {
            x[i] = x[i] - tr_m[i][j]*x[j];
        }

        if (tr_m[i][i] == 0){
            cout<<"Mathematical Error";
            exit(0);
        }    
        x[i] = x[i]/tr_m[i][i];
    }
    return x;
}


int main(){
    Matrix matrix(3, 4);
    matrix+=3;
    matrix[0][1] = 5;
    matrix[1][2] = 7;
    // matrix[][]
    Matrix m2 = ~matrix;
    matrix+=4*matrix;
    matrix /= 4;
    cout<<"\nOrig: ";
    matrix.print();
    cout<<"\nGaussian: ";
    matrix.GaussianElimination();
    matrix *= m2;
    cout<<"\nMultiplied:";
    matrix.print();
    // cout<<matrix;
    // Matrix a(2,2);
    // cin>>a;
    // cout<<a;
    int a;
    cin>>a;
    cout<<float(00)/a;
    return 0;
}
