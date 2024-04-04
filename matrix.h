#pragma once
#include <assert.h>
#include <initializer_list>
#include <cmath>
#include <string>
#include <fstream>

using ldouble = long double;
const int MAXLEN=1000000;



class Matrix{
    protected:
        int nrows;
        int ncols;
        ldouble* matrix;
    public:
        Matrix();
        Matrix(int, int);
        Matrix(int, int, ldouble);
        Matrix(std::initializer_list <std::initializer_list<ldouble> > m);
        Matrix(const Matrix&); // copy constructor
        Matrix(Matrix&&); // move constructor
        int rows() const;
        int cols() const;
        const ldouble& operator() (int, int) const; // access element at position (no mod)
        ldouble& operator() (int, int); // access element at position (modify)
        friend Matrix operator+(Matrix&, Matrix&); // addition
        friend Matrix operator-(Matrix&, Matrix&); // subtraction
        friend Matrix operator*(Matrix&, Matrix&); // multiplication
        friend Matrix operator*(ldouble, Matrix&); // scaling
        friend Matrix operator/(Matrix&, ldouble); // division by scalar
        friend Matrix operator*=(Matrix&, ldouble);
        friend Matrix operator/=(Matrix&, ldouble);
        friend void multiply(Matrix&, Matrix&, Matrix&); // multiply in-place
        Matrix& T(); // transpose
        ldouble sum() const; // return sum of the matrix
        bool is_sym() const; // returns true if matrix is symmetric
        bool is_antisym() const; // check antisym
        void resize(int, int); // change dims
        void operator=(std::initializer_list <std::initializer_list<ldouble> > m); // replace values with a list
        void operator=(Matrix&);
        void print() const; // print matrix
        void output(std::string) const; // output to a file
        ~Matrix();
};

//derived class of diagonal matrices
class Diagonal : public Matrix{
    public:
        Diagonal();
        Diagonal(int);
        Diagonal(int, ldouble);
};

// derived class of 1d matrices (i.e. vectors)
class Vector : public Matrix{
    public:
        Vector();
        Vector(int);
        Vector(int, ldouble);
        const ldouble& operator() (int) const; 
        ldouble& operator() (int); 
        void operator=(Vector&);
        void resize(int);
};

Diagonal::Diagonal() : Matrix(){};
Diagonal::Diagonal(int x): Matrix(x, x){};
Diagonal::Diagonal(int x, ldouble val): Matrix(x,x,.0){
    assert(-1 < x && x < MAXLEN);
    for (int i=0; i < this->rows(); i++)
        (*this)(i,i) = val;
}
    
Vector::Vector() : Matrix(){};
Vector::Vector(int x): Matrix(1, x){};
Vector::Vector(int x, ldouble val): Matrix(1,x,val){};

inline const ldouble& Vector::operator()(int i) const{
    assert(-1<i && i < this->cols());
    return matrix[i];
}
ldouble& Vector::operator()(int i){
    assert(-1<i && i < this->cols());
    return matrix[i];
}

void Vector::operator=(Vector& V){
    assert(ncols == V.ncols);
    for (int i=0; i<V.ncols; i++)
        matrix[i] = V(i);
}

void Vector::resize(int x){
    assert(x > 0 && x < MAXLEN);
    if (matrix == NULL){
        nrows = 1;
        ncols = x;
        matrix = new ldouble[x];
    }
    else{
        ldouble* newmatrix = new ldouble[x];
        for (int i=0; i<x; i++)
            newmatrix[i] = .0;
        for (int i=0; i<ncols; i++)
            newmatrix[i] = matrix[i];
        delete[] matrix;
        ncols = x;
        matrix = newmatrix;
    }
}

static Matrix Ainv = Matrix({{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             {2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                             { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
                             { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                             { 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0},
                             { 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0},
                             {-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
                             { 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0},
                             { 9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1},
                             {-6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1},
                             { 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
                             { 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0},
                             {-6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1},
                             { 4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1}});


Matrix::Matrix(){
    ldouble* matrix = NULL;
    nrows = 0;
    ncols = 0;
}

Matrix::~Matrix(){
    delete[] matrix;
}

// construct zero matrix of given size
Matrix::Matrix(int nr, int nc){
    assert(0<nr && 0 < nc && nr < MAXLEN && nc < MAXLEN);
    matrix = new ldouble[nr*nc];
    for (int i=0; i<nr*nc; i++)
        matrix[i] = .0;
    nrows = nr;
    ncols = nc;
}
    
// construct const matrix of given size
Matrix::Matrix(int nr, int nc, ldouble val){
    assert(0<nr && 0 < nc && nr < MAXLEN && nc < MAXLEN);
    matrix = new ldouble[nr*nc];
    for (int i=0; i<nr*nc; i++)
        matrix[i] = val;
    nrows = nr;
    ncols = nc;
}

// construct matrix with initialization values
Matrix::Matrix(std::initializer_list <std::initializer_list<ldouble> > m){
    nrows = m.size();
    ncols = m.begin()->size();
    assert(nrows < MAXLEN && ncols < MAXLEN);
    matrix = new ldouble[nrows*ncols];
    int i=0, j=0;
    for (auto v : m){
        assert(v.size() == ncols);
        for (auto e: v){
            matrix[i*ncols+j] = e;
            j++;
        }
        i++;
        j=0;
    }
 
}

// copy constructor
Matrix::Matrix(const Matrix& M){
    nrows = M.rows();
    ncols = M.cols();
    matrix = new ldouble[nrows*ncols];
    for(int i=0; i<nrows*ncols;i++)
        matrix[i] = M.matrix[i];
}

// move constructor
Matrix::Matrix(Matrix&& M){
    nrows = M.rows();
    ncols = M.cols();
    matrix = M.matrix;
    M.matrix = nullptr;
}

int Matrix::rows() const{
    return nrows;
}

int Matrix::cols() const{
    return ncols;
}

// get the value of Matrix at pos(i,j), const
inline const ldouble& Matrix::operator()(int i, int j) const{
    assert(-1<i && -1<j && i < nrows && j < ncols);
    return matrix[i*ncols + j];
}

// get the value of Matrix at pos(i,j)
ldouble& Matrix::operator()(int i, int j){
    assert(-1<i && -1<j && i < nrows && j < ncols);
    return matrix[i*ncols + j];
}

Matrix operator+(Matrix& M1, Matrix& M2){
    assert(M1.nrows==M2.nrows && M1.ncols==M2.ncols);
    Matrix result(M1.nrows, M2.nrows);
    for (int i=0; i<M1.nrows; i++)
        for (int j=0; j<M1.ncols; j++)
            result(i,j) = M1(i,j)+M2(i,j);
    return result;
}

Matrix operator-(Matrix& M1, Matrix& M2){
    assert(M1.nrows==M2.nrows && M1.ncols==M2.ncols);
    Matrix result(M1.nrows, M2.ncols);
    for (int i=0; i<M1.nrows; i++)
        for (int j=0; j<M1.ncols; j++)
            result(i,j) = M1(i,j)-M2(i,j);
    return result;
}

Matrix operator*(Matrix& M1, Matrix& M2){
    int nr1=M1.nrows;
    int nc1=M1.ncols;
    int nr2=M2.nrows;
    int nc2=M2.ncols;
    assert(nc1==nr2);
    Matrix result(nr1, nc2);
    for (int i=0; i<nr1; i++){
        for (int j=0; j<nc2; j++){
            for (int k=0; k<nc1; k++){
                result(i,j) += M1(i,k)*M2(k,j);
            }
        }
    }
    return result;
}

Matrix operator*(ldouble lambda, Matrix& M){
   Matrix result(M.rows(), M.cols());
   for (int i=0; i<M.rows(); i++)
       for (int j=0; j<M.cols(); j++)
        result(i,j) = lambda*M(i,j);
   return result;
}

Matrix operator/(Matrix& M, ldouble lambda){
    Matrix result(M.rows(), M.cols());
    for (int i=0; i<M.rows(); i++)
        for (int j=0; j<M.cols(); j++)
            result(i,j) = M(i,j)/lambda;
    return result;
}


Matrix operator*=(Matrix& M, ldouble lambda){
    for (int i=0; i<M.rows(); i++)
       for (int j=0; j<M.cols(); j++)
           M(i,j) *= lambda;
   return M;
}

Matrix operator/=(Matrix& M, ldouble lambda){
    for (int i=0; i<M.rows(); i++)
       for (int j=0; j<M.cols(); j++)
           M(i,j) /= lambda;
   return M;
}

// multiply in place (no new res matrix)
void multiply(Matrix& M1, Matrix& M2, Matrix& result){
    int nr1=M1.nrows;
    int nc1=M1.ncols;
    int nr2=M2.nrows;
    int nc2=M2.ncols;
    assert(nc1==nr2); 
    assert(result.rows() == nr1 && result.cols() == nc2);
    for (int i=0; i<nr1*nc2; i++)
        result.matrix[i] = 0;
    for (int i=0; i<nr1; i++){
        for (int j=0; j<nc2; j++){
            for (int k=0; k<nc1; k++){
                result(i,j) += M1(i,k)*M2(k,j);
            }
        }
    }
}

void Matrix::resize(int r, int c){
    assert(r*c == nrows*ncols && r > 0 && c > 0);
    nrows = r;
    ncols = c;
}

Matrix& Matrix::T(){
    ldouble* newptr = new ldouble[ncols*nrows];
    for(int i=0; i<nrows; i++)
        for (int j=0; j<ncols; j++)
            newptr[j*nrows+i] = matrix[i*ncols+j];
    delete[] matrix;
    matrix = newptr;
    return *this;
}

ldouble Matrix::sum() const{
    ldouble sum = .0;
    for (int i=0; i<nrows*ncols; i++)
        sum += matrix[i];
    return sum;
}

bool Matrix::is_sym() const{
    assert(ncols == nrows);
    ldouble usum = .0, lsum = .0, tr = .0;
    for (int i=0; i<nrows; i++){
        for (int j=0; j<ncols; j++){
            if (i > j)
                usum += matrix[i*ncols+j];
            else if (i < j)
                lsum += matrix[i*ncols+j];
            else
                tr += matrix[i*ncols+j];
        }
    }
    std::cout << "lsum: " << lsum << " usum: " << usum << " diff: " << lsum - usum << std::endl;
    if (lsum == usum){
        std::cout << "Matrix is symmetric." << std::endl;
        return true;
    }
    else{
        std::cout << "Matrix is NOT symmetric." << std::endl;
        return false;
    }
}

bool Matrix::is_antisym() const{
    assert(ncols == nrows);
    ldouble rsum = .0, lsum = .0;
    for (int i=0; i<nrows; i++){
        for (int j=0; j<ncols; j++){
            if ( (i+j) < nrows-1)
                lsum += matrix[i*ncols+j];
            else if ( (i+j) >= nrows)
                rsum += matrix[i*ncols+j];
        }
    }
    std::cout << "lsum: " << lsum << " rsum: " << rsum << " diff: " << lsum - rsum << std::endl;
    if (lsum == rsum){
        std::cout << "Matrix is anti-symmetric." << std::endl;
        return true;
    }
    else{
        std::cout << "Matrix is NOT anti-symmetric." << std::endl;
        return false;
    }
}



void Matrix::operator=(std::initializer_list <std::initializer_list<ldouble> > m){
    int nr = m.size();
    int nc = m.begin()->size();
    assert(nr < MAXLEN && nc < MAXLEN);
    assert(nr <= this->nrows && nc <= this->ncols);
    int i=0, j=0;
    for (auto v : m){
        assert(v.size() == nc);
        for (auto e: v){
            matrix[i*nc+j] = e;
            j++;
        }
        i++;
        j=0;
    }
}

void Matrix::operator=(Matrix& M){
    assert(nrows >= M.nrows && ncols >= M.ncols);
    for (int i=0; i<M.nrows; i++){
        for (int j=0; j<M.ncols; j++){
            matrix[i*ncols+j] = M(i,j);
        }
    }
}

void Matrix::print() const{
    for (int i=0; i<nrows; i++){
        for (int j=0; j<ncols; j++)
            std::cout << matrix[i*ncols+j] << " ";
        std::cout << std::endl;
    }
}

void Matrix::output(std::string str) const{
    std::ofstream out;
    out.open(str, std::ios::out);
    for(int i=0; i<nrows; i++){
        for (int j=0; j<ncols; j++)
            out << matrix[i*ncols+j] << " ";
        out << std::endl;
    }

}

ldouble bic_int(ldouble x, ldouble y, Matrix a){
    Matrix vx = {{1, x, pow(x,2), pow(x,3)}};
    Matrix vy = {{1},{y},{pow(y,2)}, {pow(y,3)}};
    Matrix tmp = vx*a;
    return (tmp*vy)(0,0);
}


Matrix bicubic(Matrix& M, int xscale, int yscale){
    int nr = M.rows();
    int nc = M.cols();
    assert(xscale >=1 && yscale >=1);
    int res_nr = nr+(nr-1)*(xscale-1);
    int res_nc = nc+(nc-1)*(yscale-1);
    Matrix result(res_nr, res_nc, -1);
    std::cout << "result matrix size: " << res_nr << "," << res_nc << std::endl;
    Matrix a(16,1);
    Matrix x(16,1, -1.0); // f
    for (int i=0; i<nr-1; i++){
        for (int j=0; j<nc-1; j++){
            x = { {M(i, j)}, //f00
                  {M(i+1, j)},
                  {M(i, j+1)},
                  {M(i+1, j+1)},
                  
                  
                  {(i > 0) ? (M(i+1,j)-M(i-1,j))/2 : (M(i+1,j)-M(i,j))}, //fx00
                  {(i < nr-2) ? (M(i+2,j)-M(i,j))/2 : (M(i+1,j)-M(i,j))}, //fx10
                  {(i > 0) ? (M(i+1,j+1)-M(i-1,j+1))/2 : (M(i+1,j+1)-M(i,j+1))}, //fx01
                  {(i < nr-2) ? (M(i+2,j+1)-M(i,j+1))/2 : (M(i+1,j+1)-M(i,j+1))}, //fx11
                  
                  /*
                  {(M(i+1,j)-M(i,j))}, //fx00
                  {(M(i+1,j)-M(i,j))}, //fx10
                  {(M(i+1,j+1)-M(i,j+1))}, //fx01
                  {(M(i+1,j+1)-M(i,j+1))}, //fx11
                  */

                  
                  {(j > 0) ? (M(i,j+1)-M(i,j-1))/2 : (M(i,j+1)-M(i,j))}, //fy00
                  {(j > 0) ? (M(i+1,j+1)-M(i+1,j-1))/2 : (M(i+1,j+1)-M(i+1,j))}, //fy10
                  {(j < nc-2) ? (M(i,j+2)-M(i,j))/2 : (M(i,j+1)-M(i,j))}, //fy01
                  {(j < nc-2) ? (M(i+1,j+2)-M(i+1,j))/2 : (M(i+1,j+1)-M(i+1,j))}, //fy11
                  

                  /*
                  {(M(i,j+1)-M(i,j))}, //fy00
                  {(M(i+1,j+1)-M(i+1,j))}, //fy10
                  {(M(i,j+1)-M(i,j))}, //fy01
                  { (M(i+1,j+1)-M(i+1,j))}, //fy11
                  */
                  
                  
                  {((i > 0) && (j > 0)) ? (M(i+1,j+1)+M(i-1,j-1)-M(i-1,j+1)-M(i+1,j-1))/4 : (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) }, //fxy00
                  //{((i < nr-2) && (j > 0)) ? (M(i+2,j+1)+M(i,j-1)-M(i+2,j-1)-M(i,j+1))/4 : ( (j==0) ? (M(i+1,j+1)+M(i,j)-M(i,j+1)-M(i+1,j)) : (M(i+1,j)+M(i,j-1)-M(i,j)-M(i+1,j-1)) )}, //fxy10
                  {((i < nr-2) && (j > 0)) ? (M(i+2,j+1)+M(i,j-1)-M(i+2,j-1)-M(i,j+1))/4 : ( (j==0) ? (M(i+1,j+1)+M(i,j)-M(i,j+1)-M(i+1,j)) : (M(i,j-1)+M(i+1,j+1)-M(i,j+1)-M(i+1,j-1)) )}, //fxy10
                  //{((i > 0) && (j < nc-2)) ? (M(i+1,j+2)+M(i-1,j)-M(i-1,j+2)-M(i+1,j))/4 : ( (i==0) ? (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) : (M(i,j+1)+M(i-1,j)-M(i,j)-M(i-1,j+1)) )}, //fxy01
                  {((i > 0) && (j < nc-2)) ? (M(i+1,j+2)+M(i-1,j)-M(i-1,j+2)-M(i+1,j))/4 : ( (i==0) ? (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) : (M(i-1,j)+M(i+1,j+1)-M(i+1,j)-M(i-1,j+1)) )}, //fxy01
                  {((i < nr-2) && (j < nc-2)) ? (M(i+2,j+2)+M(i,j)-M(i+2,j)-M(i,j+2))/4 : (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) } //fxy11
                  

                  /*
                  {(M(i+1,j+1)+M(i,j)-M(i,j+1)-M(i+1,j)) }, //fxy00
                  {(M(i+1,j+1)+M(i,j)-M(i,j+1)-M(i+1,j)) }, //fxy10
                  {(M(i+1,j+1)+M(i,j)-M(i,j+1)-M(i+1,j)) }, //fxy01
                  {(M(i+1,j+1)+M(i,j)-M(i,j+1)-M(i+1,j)) }, //fxy11
                  */
                  

                  /*
                  {((i > 0) && (j > 0)) ? (M(i+1,j+1)+M(i-1,j-1)-M(i-1,j+1)-M(i+1,j-1))/4 : (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) }, //fxy00
                  {((i < nr-2) && (j > 0)) ? (M(i+2,j+1)+M(i,j-1)-M(i+2,j-1)-M(i,j+1))/4 : (M(i+1,j+1)+M(i,j)-M(i,j+1)-M(i+1,j)) }, //fxy10
                  {((i > 0) && (j < nc-2)) ? (M(i+1,j+2)+M(i-1,j)-M(i-1,j+2)-M(i+1,j))/4 : (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) }, //fxy01
                  {((i < nr-2) && (j < nc-2)) ? (M(i+2,j+2)+M(i,j)-M(i+2,j)-M(i,j+2))/4 : (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) } //fxy11
                  */
            };
            //x.print();
            multiply(Ainv, x, a);
            a.resize(4,4);
            a.T();

            // now interpolate any points inside the grid
            ldouble r_range=(ldouble)xscale, c_range=(ldouble)yscale;
            if (i==nr-2)
                r_range += 1.0;            
            if (j==nc-2)
                c_range += 1.0;

            for (ldouble r=0; r<r_range; r++){
                for (ldouble c=0; c<c_range; c++){
                    result(i*xscale+r, j*yscale+c) = bic_int(r/(ldouble)xscale, c/(ldouble)yscale, a);
                }
            }
            a.resize(16,1);
        }
    }

    return result;
}

// mem efficient version, simply calculate the sum instead of the matrix
ldouble bicubic_mem_eff(Matrix& M, int xscale, int yscale){
    int nr = M.rows();
    int nc = M.cols();
    assert(xscale >=1 && yscale >=1);
    ldouble res_sum = .0; // calc sum of result
    Matrix a(16,1);
    Matrix x(16,1, -1.0); // f
    for (int i=0; i<nr-1; i++){
        for (int j=0; j<nc-1; j++){
            x = { {M(i, j)}, //f00
                  {M(i+1, j)},
                  {M(i, j+1)},
                  {M(i+1, j+1)},

                  {(i > 0) ? (M(i+1,j)-M(i-1,j))/2 : (M(i+1,j)-M(i,j))}, //fx00
                  {(i < nr-2) ? (M(i+2,j)-M(i,j))/2 : (M(i+1,j)-M(i,j))}, //fx10
                  {(i > 0) ? (M(i+1,j+1)-M(i-1,j+1))/2 : (M(i+1,j+1)-M(i,j+1))}, //fx01
                  {(i < nr-2) ? (M(i+2,j+1)-M(i,j+1))/2 : (M(i+1,j+1)-M(i,j+1))}, //fx11

                  {(j > 0) ? (M(i,j+1)-M(i,j-1))/2 : (M(i,j+1)-M(i,j))}, //fy00
                  {(j > 0) ? (M(i+1,j+1)-M(i+1,j-1))/2 : (M(i+1,j+1)-M(i+1,j))}, //fy10
                  {(j < nc-2) ? (M(i,j+2)-M(i,j))/2 : (M(i,j+1)-M(i,j))}, //fy01
                  {(j < nc-2) ? (M(i+1,j+2)-M(i+1,j))/2 : (M(i+1,j+1)-M(i+1,j))}, //fy11
                  
                  {((i > 0) && (j > 0)) ? (M(i+1,j+1)+M(i-1,j-1)-M(i-1,j+1)-M(i+1,j-1))/4 : (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) }, //fxy00
                  {((i < nr-2) && (j > 0)) ? (M(i+2,j+1)+M(i,j-1)-M(i+2,j-1)-M(i,j+1))/4 : ( (j==0) ? (M(i+1,j+1)+M(i,j)-M(i,j+1)-M(i+1,j)) : (M(i,j-1)+M(i+1,j+1)-M(i,j+1)-M(i+1,j-1)) )}, //fxy10
                  {((i > 0) && (j < nc-2)) ? (M(i+1,j+2)+M(i-1,j)-M(i-1,j+2)-M(i+1,j))/4 : ( (i==0) ? (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) : (M(i-1,j)+M(i+1,j+1)-M(i+1,j)-M(i-1,j+1)) )}, //fxy01
                  {((i < nr-2) && (j < nc-2)) ? (M(i+2,j+2)+M(i,j)-M(i+2,j)-M(i,j+2))/4 : (M(i+1,j+1)+M(i,j)-M(i+1,j)-M(i,j+1)) } //fxy11
            };
            multiply(Ainv, x, a);
            a.resize(4,4);
            a.T();

            // now interpolate any points inside the grid
            ldouble r_range=(ldouble)xscale, c_range=(ldouble)yscale;
            if (i==nr-2)
                r_range += 1.0;
            if (j==nc-2)
                c_range += 1.0;

            for (ldouble r=0; r<r_range; r++){
                for (ldouble c=0; c<c_range; c++){
                    res_sum += bic_int(r/(ldouble)xscale, c/(ldouble)yscale, a)/xscale/yscale;
                }
            }
            a.resize(16,1);
        }
    }
    return res_sum;
}

/*
Vector cubic(Vector& V, int scale){
}

// cubic interpolation with their X-coordinates (i.e. non-uniform intervals)
Vector cubic(Vector& V, Vector&X, int scale){

}
*/

