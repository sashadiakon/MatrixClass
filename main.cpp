#include "main.h"



float get_rand(float min, float max){
    return min + rand()%int(max-min);
}

int sign(double a) {
    if (a > 0) return 1;
    if (a < 0) return -1;
    return 0;
}

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

int Matrix::getRows(){
    return rows;
}

void Matrix::Identity(){
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (i==j){
                values[i][j] = 1;
            }
            else {
                values[i][j] = 0;
            }
        }
    }
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

bool Matrix::is_diag(){
    bool res = true;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; ++j) {
            if (i!=j){
                res &= (values[i][j] == 0);
            }
        }
    }
    return res;
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

vector<float> Matrix::GaussianElimination(){
    // Check size
    
    if (rows+1 != cols){
        cout<<"Not apropriate size, matrix should have size - (n, n+1)\n";
        exit(0);
    }
    vector<float> x(rows,0);
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

// function to find max non diagonal element
pair<int, int> Matrix::nondiag(){
    float max = 0;
    pair<int, int> ij;
    for (int i =0; i < rows; i++){
        for(int j = 1; j < cols; j++){
            if ((i != j) && (fabs(values[i][j])>max)){
                max = fabs(values[i][j]);
                ij.first = i;
                ij.second = j;
            }
        }
    }
    return ij;
}

void Matrix::Rotate(int k, int l, int i, int j) {   
    cout<<1;
}


vector<float> Matrix::JacobiAlg(float eps =0.01) {
    int i, j, p, q, flag;
    
    float theta, t;
    t = INFINITY;
    Matrix A = *this;
    Matrix U(rows, rows);
    
    while (t>eps){
        pair<int, int> coor = A.nondiag();
        i = coor.first;
        j = coor.second;
        if (A[i][i] == A[j][j]){
            if (A[i][j]> 0){
                theta = M_PI/4;
            }
            else{
                theta = -M_PI/4;
            }
        }
        else{
            theta = 0.5*atan(2*A[i][j]/(A[i][i]-A[j][j]));
        }
        U.Identity();
        U[i][i] = cos(theta);
        U[j][j] = cos(theta);
        U[i][j] = sin(theta);
        U[j][i] = -sin(theta);
        A = U*A*~U;
        t = 0;
        for (int i =0; i < rows; i++){
            for(int j = 0; j < cols; j++){
                if (i!=j){
                    t+= A[i][j]*A[i][j];
                }
            }
        }
        // cout<<"U: \n"<<U<<endl;
        // cout<<"A:\n"<<A<<endl;
    }
    vector<float> res;
    for (i = 0; i<rows; i++){
        res.push_back(A[i][i]); 
    }
    return res;
}

void Matrix::RandomSymetric(float min = 0, float max = 100){
    for (int i = 0; i < rows; i++){
        for (int j = 0; j<i; j++){
            values[i][j] = values[j][i];
        }
        for (int j = i; j < cols; j++){
            values[i][j] = get_rand(min, max);
        }
    }
}

int main(){
    // Set seed for random gen
    srand((unsigned) time(NULL));

    Matrix matrix(4, 4);
    matrix[0] = vector<float> {2,-1,0,0};
    matrix[1] = vector<float> {-1, 2, -1, 0};
    matrix[2] = vector<float> {0, -1, 2, -1};
    matrix[3] = vector<float> {0,0,-1, 2};
    vector<float> eigenvalues = matrix.JacobiAlg(0.1);
    for (float x: eigenvalues){
        cout<<x<<" ";
    }

    // Test for random 5x5 matrix
    Matrix m(5,5);
    m.RandomSymetric();
    cout<<m;
    vector<float> eigenvalues2 = m.JacobiAlg(0.1);
    for (float x: eigenvalues2){
        cout<<x<<" ";
    }

    // Test of GaussianElimination
    Matrix m2(3,4);
    m2+=3;
    m2[0][1] +=4; m2[1][3] += 2; m2[2][0] = 7;
    cout<<"\nGaussian: ";
    vector<float> xs = m2.GaussianElimination();
    for (int i =0; i < m2.getRows(); i++){
        cout<<xs[i]<<" ";
    }
    return 0;
}
