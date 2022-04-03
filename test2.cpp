#include "Math.h"
#include <iostream>

int main()
{
    using P = ScalarCfv::point;
    ScalarCfv::tensor2D<ScalarCfv::real, 7, 6> data;
    ScalarCfv::tensor1D<ScalarCfv::real, 7> moment;
    ScalarCfv::cellFieldData c;
    c.cellType_ = ScalarCfv::Quadrilateral;
    c.cellNode.resize(5);
    c.cellNode[1].second = P(0, 0);
    c.cellNode[2].second = P(0, 1);
    c.cellNode[3].second = P(1, 1);
    c.cellNode[4].second = P(1, 0);
    int nx = 20;
    int ny = 20;
    std::cout << '[';
    for (int ix = 0; ix <= nx; ix++)
        for (int iy = 0; iy <= ny; iy++)
        {
            CfvMath::getDiffBaseValueRBFB1_7_6(P(ix * 1.0 / nx, iy * 1.0 / nx),
                                               P(0.5, 0.5),
                                               P(1, 1),
                                               moment, data, c);
            for (int a = 0; a < 7; a++)
                for (int b = 0; b < 6; b++)
                    std::cout << data[a][b] << '\t';
            std::cout << ";\n";
        }
    std::cout << ']';

    Eigen::MatrixXd A(5, 5);
    A << 0.543404941790965, 0.121569120783114, 0.891321954312264, 0.978623784707370, 0.431704183663122,
        0.278369385093796, 0.670749084726779, 0.209202122117190, 0.811683149089323, 0.940029819622375,
        0.424517590749133, 0.825852755105048, 0.185328219550075, 0.171941012732594, 0.817649378776727,
        0.844776132319904, 0.136706589684953, 0.108376890464255, 0.816224748725840, 0.336111950120899,
        0.00471885619097257, 0.575093329427250, 0.219697492624992, 0.274073747041699, 0.175410453742337;
    A(Eigen::all, 0).setZero();
    A(0, Eigen::all).setZero();
    Eigen::MatrixXd B(5, 5);
    CfvMath::EigenLeastSquareInverse(A, B);
    std::cout << B << std::endl;

    A.resize(4, 6);
    A << 0.598843376928493, 0.0364760565925689, 0.890545944728504, 0.581842192398778, 0.769115171105652, 0.975006493606588,
        0.603804539042854, 0.890411563442076, 0.576901499400033, 0.0204391320269232, 0.250695229138396, 0.884853293491106,
        0.105147685412056, 0.980920857012312, 0.742479689097977, 0.210026577672861, 0.285895690406865, 0.359507843936902,
        0.381943444943110, 0.0599419888180373, 0.630183936475376, 0.544684878178648, 0.852395087841306, 0.598858945875747;
    B.resize(6, 4);
    CfvMath::EigenLeastSquareInverse(A, B);
    std::cout << "Ainv : \n"
              << B << std::endl;

    CfvMath::EigenLeastSquareInverse(A.transpose() * A, B);
    std::cout << "ATAinv : \n"
              << B << std::endl;

    return 0;
}