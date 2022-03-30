#include "Math.h"


namespace CfvMath
{
    void EigenLeastSquareInverse(const Eigen::MatrixXd &A, Eigen::MatrixXd &AI)
    {
        auto SVDResult = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
        AI = SVDResult.solve(Eigen::MatrixXd::Identity(AI.cols(), AI.rows()));
    }
}
