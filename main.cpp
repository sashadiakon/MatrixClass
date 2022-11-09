#include "main.h"



int get_rand(int min = 0, int max = 10){
    return min + rand()%int(max-min);
}

float rand_float(int min = 0, int max = 7001, int m = 7001){
    return float(get_rand(0, 7001))/7001;
}

float rand_diviation(int min = 0, int max = 7001, int m = 7001){
    return rand_float(min, max, m)-0.5;
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

Matrix::Matrix(vector <float> m){
    values = vector<vector<float>>(1, m);
    rows = 1;
    cols = m.size();
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

// Calculate triangular form for matrix
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
                tr_m[i][j]-=tr_m[b][j]*f;
            }
        }
    }
    return tr_m;
}

// return Cofactor of the matrix 
Matrix Matrix::getCofactor(int n, int p, int q)
{
    int i = 0, j = 0;
    Matrix temp(rows, rows);
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those
            //  element which are not in given row and
            //  column
            if (row != p && col != q)
            {
                temp[i][j++] = values[row][col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
    return temp;
}

// Calculate Minor
float Matrix::Minor(int n, int i, int j){
    return this->getCofactor(n, i, j).Det(n-1); 
}

// Calculate Determinant
float Matrix::Det(int n){
    if (not this->is_square()){
        cout<<"We cannot calculate it for non square matrix";
        exit(0);
    }
    float D = 0;
    if (n == 1){
        return values[0][0];
    }
    for (int i =0; i< n; i++){
        D += pow(-1, i)*values[0][i]*Minor(n, 0, i);
    }
    return D;
}

// Calc inverse matrix
Matrix Matrix::Inverse(){
    if (not this->is_square()){
        cout<<"We cannot calculate it for non square matrix";
        exit(0);
    }
    Matrix inv(rows, cols);
    float det = this->Det(rows);
    for (int i =0; i<rows; i++){
        for(int j=0; j<cols; j++){
            inv[i][j] = pow(-1, i+j)*Minor(rows, i, j);
        }
    }
    inv /= det;
    return inv;
}

// Gaussian elimination
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

// find eigenvalues with Jacobi method
vector<float> Matrix::JacobiAlg(float eps =0.01) {
    if (not this->is_square()){
        cout<<"We cannot calculate it for non square matrix";
        exit(0);
    }
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
    }
    vector<float> res;
    for (i = 0; i<rows; i++){
        res.push_back(A[i][i]); 
    }
    return res;
}

// transform matrix to symetric matrix with random coefs
void Matrix::RandomSymetric(float min = 0, float max = 100){
    if (not this->is_square()){
        cout<<"We cannot calculate it for non square matrix";
        exit(0);
    }
    for (int i = 0; i < rows; i++){
        for (int j = 0; j<i; j++){
            values[i][j] = values[j][i];
        }
        for (int j = i; j < cols; j++){
            values[i][j] = get_rand(min, max);
        }
    }
}

// generate random coefficients for linear function
vector <float> GenerateLinearFunc(int n = 1, int min = -100, int max = 100){
    vector <float> coefs;
    for (int i = 0; i <= n; i++){
        coefs.push_back(get_rand(min, max));
    }
    return coefs;
    
}

// func to generate data for linear function with small deviation from real values
auto TrainDataGenerator(vector <float> func, int n_samples = 30){
    Matrix coef(func);
    coef = ~coef;
    Matrix X(n_samples, func.size());
    Matrix y(n_samples, 1);
    for (int i =0; i<n_samples; i++){
        for(int x = 0; x<func.size();x++){
            if (x == 0){
                X[i][x] = 1;
            }
            else{
                X[i][x] = get_rand(-100, 100);
            }
        }
        float div = rand_diviation();
        y[i][0] = (Matrix(X[i])*coef)[0][0]+div;
    }
    auto data = make_pair(X, y);
    data.first = X;
    data.second = y;
    return data;
}

// Use LeastSquares method to approximate the linear func coefs
Matrix LeastSquares(Matrix X_train, Matrix Y_train){
    Matrix weights = (~X_train*X_train).Inverse()*~X_train*Y_train;
    return weights;
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
    cout<<"Jakobi algorithm tests: \n";
    for (float x: eigenvalues){
        cout<<x<<" ";
    }

    // Test for random 5x5 matrix
    Matrix m(3,3);
    m.RandomSymetric();
    cout<<endl<<m;
    vector<float> eigenvalues2 = m.JacobiAlg(0.1);
    for (float x: eigenvalues2){
        cout<<x<<" ";
    }

    // Test of GaussianElimination
    Matrix m2(3,4);
    m2+=3;
    m2[0][1] +=4; m2[1][3] += 2; m2[2][0] = 7;
    cout<<"\nGaussian: \n";
    vector<float> xs = m2.GaussianElimination();
    for (int i =0; i < m2.getRows(); i++){
        cout<<xs[i]<<" ";
    }

    // Test of linear regresion 
    cout<<"\nLinearReg:\n";
    vector <float> func = GenerateLinearFunc(3);
    auto data = TrainDataGenerator(func);
    Matrix x_train = data.first;
    Matrix y_train = data.second;
    cout<<"Real func: f= "<<func[0];
    for (int i = 1; i < func.size(); i++){
        if (func[i] >= 0){
            cout<<" + ";
        }
        else cout<<" - ";
        cout<<fabs(func[i])<<"x"<<i;
    }
    // predict function parameters
    Matrix weights(LeastSquares(x_train, y_train));
    cout<<"\nPredicted: \n f = "<<weights[0][0];
    for (int i =1; i < weights.getRows(); i++){
        if (weights[i][0] >= 0){
            cout<<" + ";
        }
        else cout<<" - ";
        cout<<fabs(weights[i][0])<<"x"<<i;
    }
    return 0;
}
