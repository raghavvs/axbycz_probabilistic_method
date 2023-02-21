#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <vector>

using namespace Eigen;
using namespace std;

template <typename DerivedA, typename DerivedB, typename DerivedC>
void batchSolveXY(const MatrixBase<DerivedA>& A, const MatrixBase<DerivedB>& B, const string& opt, const double nstd1, const double nstd2,
                   MatrixBase<DerivedC>& X, MatrixBase<DerivedC>& MeanA, MatrixBase<DerivedC>& MeanB, MatrixBase<DerivedC>& MeanC) {

    const int Num = A.derived().cols();

    // Calculate B_inv and A_inv
    vector<Matrix4d> A_inv(Num);
    vector<Matrix4d> B_inv(Num);
    for (int i = 0; i < Num; ++i) {
        A_inv[i] = A.col(i).inverse();
        B_inv[i] = B.col(i).inverse();
    }

    X.resize(4, 4, 4);
    MeanA.resize(4, 4);
    MeanB.resize(4, 4);
    MeanC.resize(4, 4);

    if (opt == "med") {
        // Calculate median values of A, B_inv
        for (int i = 0; i < 4; ++i) {
            vector<double> av, bv;
            for (int j = 0; j < Num; ++j) {
                av.push_back(A(i, j));
                bv.push_back(B_inv[j](i, i));
            }
            nth_element(av.begin(), av.begin() + av.size() / 2, av.end());
            MeanA(i, i) = av[av.size() / 2];
            nth_element(bv.begin(), bv.begin() + bv.size() / 2, bv.end());
            MeanB(i, i) = 1.0 / bv[bv.size() / 2];
        }

        // Calculate mean value of C
        MeanC = Matrix4d::Zero();
        for (int j = 0; j < Num; ++j) {
            MeanC += A.col(j) * MeanB * B.col(j).transpose();
        }
        MeanC /= Num;
        MeanC(3, 3) = 1.0;

        // Compute X
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < Num; ++k) {
                    X(i, j, k) = A(i, k) * MeanB(j, j);
                }
            }
        }

    } else if (opt == "mean") {
        // Calculate mean values of A, B_inv
        for (int i = 0; i < 4; ++i) {
            double sumA = 0.0, sumB = 0.0;
            for (int j = 0; j < Num; ++j) {
                sumA += A(i, j);
                sumB += B_inv[j](i, i);
            }
            MeanA(i, i) = sumA / Num;
            MeanB(i, i) = 1.0 / (sumB / Num);
        }

        // Calculate mean value of C
        MeanC = Matrix4d::Zero();
        for (int j = 0; j < Num; ++j) {
            MeanC += A.col(j) * MeanB * B.col(j).transpose
