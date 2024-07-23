#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <chrono>
#include <functional>
#include <cmath>
#include "./eigen-3.4.0/Eigen/Dense"
using namespace std;

typedef long double ld;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> ld_mat;
typedef Eigen::Matrix<complex<ld>, Eigen::Dynamic, Eigen::Dynamic> complex_ld_mat;

unsigned int seed=std::random_device{}();
std::mt19937 rdm(seed);

ld_mat matrixPower(ld_mat A, int n){
    ld_mat ans(A.cols(),A.rows()), tmp=A;
    ans.setIdentity();
    for(;n;n>>=1, tmp=tmp*tmp) if(n&1) ans=ans*tmp;
    return ans;
}

complex<ld> fastPower(complex<ld> a, int b){
    complex<ld> ans=a, tmp=a;
    b--;
    for(;b;b>>=1, tmp=tmp*tmp) if(b&1) ans=ans*tmp;
    return ans;
}

class UndiagonalizedMatrix{
    public:
        UndiagonalizedMatrix(int n){
            std::mt19937 rdm(seed);
            mtx.resize(n,n);
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++) mtx(i,j)=rdm()%13;
            }
        }

        void calculate(int n){
            //cout<<"Undiagonalized原矩陣：\n"<<mtx<<endl<<endl;

            auto start = std::chrono::high_resolution_clock::now();
            ans=matrixPower(mtx,n);
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

            cout << "冪運算時長：\n" << duration.count() << " ms" << endl;
            //cout<<"Undiagonalized計算結果：\n\n"<<ans<<endl<<endl;
            cout<<"校驗：\n"<<ans(0,0)<<endl<<endl;
        }

    private:
        ld_mat mtx;
        ld_mat ans;
};


class DiagonalizedMatrix{
    public:
        DiagonalizedMatrix(int n){
            std::mt19937 rdm(seed);
            mtx.resize(n,n);
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++) mtx(i,j)=rdm()%13;
            }
        }

        void diagonalize(){
            auto start = std::chrono::high_resolution_clock::now();

            Eigen::EigenSolver<ld_mat> eigensolver(mtx);
            if(eigensolver.info() != Eigen::Success) abort();

            eigen_val.setZero(mtx.rows(),mtx.rows());
            for(int i=0;i<mtx.rows();i++) eigen_val(i,i)=eigensolver.eigenvalues()[i];
            eigen_vec = eigensolver.eigenvectors();

            //cout<<"Diagonalized原矩陣：\n"<<mtx<<endl<<endl;
            //cout<<"eigenvalues對角矩陣：\n"<<eigen_val<<endl<<endl;
            //cout<<"eigenvector：\n"<<eigen_vec<<endl<<endl;

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            cout << "對角化時長：\n" << duration.count() << " ms" << endl;
        }

        void caluclate(int n, std::function<complex<ld> (complex<ld>, int)> powFunc){
            auto start = std::chrono::high_resolution_clock::now();

            for(int i=0;i<mtx.rows();i++){
                eigen_val(i,i)=powFunc(eigen_val(i,i), n);
            }
            ans=eigen_vec*eigen_val*eigen_vec.inverse();

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

            cout << "冪運算時長：\n" << duration.count() << " ms" << endl;
            //cout<<"Diagonalized計算結果：\n"<<ans<<endl;
            cout<<"校驗：\n"<<ans(0,0)<<endl;
        }
    
    private:
        ld_mat mtx;
        complex_ld_mat eigen_val;
        complex_ld_mat eigen_vec;
        complex_ld_mat ans;
};


signed main(){

    UndiagonalizedMatrix undiagonalized_matrix(2);
    undiagonalized_matrix.calculate(4500);

    DiagonalizedMatrix diagonalized_matrix(2);
    diagonalized_matrix.diagonalize();
    diagonalized_matrix.caluclate(4500,fastPower);

    return 0;
}