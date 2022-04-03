#pragma once

#include "Math.h"

#include "RBF.h"

#include "Point.h"

namespace CfvMath
{
    ScalarCfv::real rbfScale = 1;
    /*

    3-3

    */

    bool getMomentRBFB1_POLY(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 3> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {
#ifdef RBFB1_GlobalPoly
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        A[1] = delta.x / scale.x;
        A[2] = delta.y / scale.y;
#else
        A[1] = p.x - 0.5;
        A[2] = p.y - 0.5;
#ifdef TRIAL
        // std::cout << "Moment in " << p.x << p.y << std::endl;
#endif

#endif
        return true;
    }

    bool getBaseValueRBFB1_POLY(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 3> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 3> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
#ifdef RBFB1_GlobalPoly
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        A[1] = delta.x / scale.x - moment[1];
        A[2] = delta.y / scale.y - moment[2];
#else
        A[1] = p.x - 0.5 - moment[1];
        A[2] = p.y - 0.5 - moment[2];

#endif
        return true;
    }

    bool getDiffBaseValueRBFB1_POLY(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 3> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 3, 3> &A,
        ScalarCfv::cellFieldData &cell) // adding 0
    {

        //(1-x)(1-y) y(1-x) x*y x*(1-y)
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};

#ifdef RBFB1_GlobalPoly
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        A[1][0] = delta.x / scale.x - moment[1];
        A[2][0] = delta.y / scale.y - moment[2];

        A[1][1] = 1 / scale.x;
        A[1][2] = 0 / scale.x;
        A[2][1] = 0 / scale.y;
        A[2][2] = 1 / scale.y;
#else

        A[1][0] = p.x - 0.5 - moment[1];
        A[2][0] = p.y - 0.5 - moment[2];

        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj

        // diff 1:
        Eigen::Vector2d dphidetaj, dphidxj;
        dphidetaj << 1, 0;
        dphidxj = iJacobi * dphidetaj;
        A[1][1] = dphidxj(0);
        A[1][2] = dphidxj(1);

        dphidetaj << 0, 1;
        dphidxj = iJacobi * dphidetaj;
        A[2][1] = dphidxj(0);
        A[2][2] = dphidxj(1);

        auto a = iJacobi.determinant();
        auto ia = Jacobi.determinant();
#endif
        // if (abs(baryCenter.x - 0.552) < 1e-2 && abs(baryCenter.y - 3.056) < 1e-1)
        //     return false;
        // if (abs(baryCenter.x - 0.0833) < 1e-2 && abs(baryCenter.y - 3.038) < 1e-1)
        //     return false;

        return true;
    }

    /*

    4-6

    */

    bool getMomentRBFB1_4_6(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {
        // Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        // Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
        //                                      {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        // Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
        //                                  {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        // Eigen::Vector2d pp = XiNj * Nj;
        // ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        // A[1] = delta.x / scale.x;
        // A[2] = delta.y / scale.y;
        auto pc = p - ScalarCfv::point(0.5, 0.5);
        A[1] = pc.x;
        A[2] = pc.y;
        A[3] = RBF::RBF0(pc, crbf, RBF::F0);
        return true;
    }

    bool getBaseValueRBFB1_4_6(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {

        auto pc = p - ScalarCfv::point(0.5, 0.5);
        A[1] = pc.x - moment[1];
        A[2] = pc.y - moment[2];
        A[3] = RBF::RBF0(pc, crbf, RBF::F0) - moment[3];

        return true;
    }

    bool getDiffBaseValueRBFB1_4_6(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 4, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {
        auto pc = p - ScalarCfv::point(0.5, 0.5);
        A[1][0] = pc.x - moment[1];
        A[2][0] = pc.y - moment[2];
        A[3][0] = RBF::RBF0(pc, crbf, RBF::F0) - moment[3];

        assert(cell.cellType_ == ScalarCfv::Quadrilateral);

        //(1-x)(1-y) y(1-x) x*y x*(1-y)
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};

        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};

        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj

        // Eigen::Vector2d pp = XiNj * Nj;
        // ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        // A[1][0] = delta.x / scale.x;
        // A[1][1] = 1 / scale.x;
        // A[1][2] = 0;
        // A[2][0] = delta.y / scale.y;
        // A[2][1] = 0;
        // A[2][2] = 1 / scale.y;

        // diff 1:
        Eigen::Vector2d dphidetaj, dphidxj;
        dphidetaj << 1, 0;
        dphidxj = iJacobi * dphidetaj;
        A[1][1] = dphidxj[0];
        A[1][2] = dphidxj[1];

        dphidetaj << 0, 1;
        dphidxj = iJacobi * dphidetaj;
        A[2][1] = dphidxj[0];
        A[2][2] = dphidxj[1];

        dphidetaj << RBF::RBF1x(pc, crbf, RBF::F0, RBF::F1), RBF::RBF1y(pc, crbf, RBF::F0, RBF::F1);
        dphidxj = iJacobi * dphidetaj;
        A[3][1] = dphidxj[0];
        A[3][2] = dphidxj[1];

        // diff 2:
        A[1][3] = A[1][4] = A[1][5] = 0;
        A[2][3] = A[2][4] = A[2][5] = 0;
        Eigen::Matrix2d dphidetaidetaj, dphidxidxj;
        ScalarCfv::real rbfxx, rbfxy, rbfyy;
        rbfxx = RBF::RBF2xx(pc, crbf, RBF::F0, RBF::F1, RBF::F2);
        rbfxy = RBF::RBF2xy(pc, crbf, RBF::F0, RBF::F1, RBF::F2);
        rbfyy = RBF::RBF2yy(pc, crbf, RBF::F0, RBF::F1, RBF::F2);
        dphidetaidetaj << rbfxx, rbfxy, rbfxy, rbfyy;
        dphidxidxj = iJacobi * dphidetaidetaj * iJacobi.transpose();
        A[3][3] = dphidxidxj(0, 0);
        A[3][4] = dphidxidxj(0, 1);
        A[3][5] = dphidxidxj(1, 1);

#ifdef DEBUG_RBFB
        std::cout << "XiNj\n"
                  << XiNj << std::endl;
        std::cout << "Jacobi\n"
                  << Jacobi << std::endl;
        std::cout << Jacobi.determinant() << std::endl;
#endif

        return true;
    }

    bool getMomentRBFB1_4_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {
        // Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        // Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
        //                                      {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        // Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
        //                                  {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        // Eigen::Vector2d pp = XiNj * Nj;
        // ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        // A[1] = delta.x / scale.x;
        // A[2] = delta.y / scale.y;

#ifdef RBFB1_RBF_USE_MEANIJ
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d iJacobi = cell.MeanIJ;
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
        Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
        ScalarCfv::point pc(ppc(0), ppc(1));
#else
        auto pc = p - ScalarCfv::point(0.5, 0.5);
#endif

        A[1] = pc.x;
        A[2] = pc.y;
        A[3] = RBF::RBF0(pc, crbf, RBF::F0);
#ifdef RBFB1_RBF_USE_MEANIJ
        auto meanJ = iJacobi.inverse();
        auto L0 = meanJ(0, Eigen::all).norm();
        auto L1 = meanJ(1, Eigen::all).norm();
        if (L0 > RBFB1_AR0 * L1 || L1 > RBFB1_AR0 * L0)
            A[3] *= 0;
        else
            std::cout << "LS" << L0 << '\t' << L1 << std::endl;
#endif
        return true;
    }

    bool getBaseValueRBFB1_4_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 4> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {

#ifdef RBFB1_RBF_USE_MEANIJ
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d iJacobi = cell.MeanIJ;
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
        Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
        ScalarCfv::point pc(ppc(0), ppc(1));
#else
        auto pc = p - ScalarCfv::point(0.5, 0.5);
#endif
        A[1] = pc.x - moment[1];
        A[2] = pc.y - moment[2];
        A[3] = RBF::RBF0(pc, crbf, RBF::F0) - moment[3];
#ifdef RBFB1_RBF_USE_MEANIJ
        auto meanJ = iJacobi.inverse();
        auto L0 = meanJ(0, Eigen::all).norm();
        auto L1 = meanJ(1, Eigen::all).norm();
        if (L0 > RBFB1_AR0 * L1 || L1 > RBFB1_AR0 * L0)
            A[3] *= 0;
#endif

        return true;
    }

    bool getDiffBaseValueRBFB1_4_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 4> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 4, 3> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);

        //(1-x)(1-y) y(1-x) x*y x*(1-y)
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};

        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};

        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj

#ifdef RBFB1_RBF_USE_MEANIJ
        iJacobi = cell.MeanIJ;
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
        Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
        ScalarCfv::point pc(ppc(0), ppc(1));
#else
        auto pc = p - ScalarCfv::point(0.5, 0.5);
#endif
        A[1][0] = pc.x - moment[1];
        A[2][0] = pc.y - moment[2];
        A[3][0] = RBF::RBF0(pc, crbf, RBF::F0) - moment[3];

        // Eigen::Vector2d pp = XiNj * Nj;
        // ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        // A[1][0] = delta.x / scale.x;
        // A[1][1] = 1 / scale.x;
        // A[1][2] = 0;
        // A[2][0] = delta.y / scale.y;
        // A[2][1] = 0;
        // A[2][2] = 1 / scale.y;

        // diff 1:
        Eigen::Vector2d dphidetaj, dphidxj;
        dphidetaj << 1, 0;
        dphidxj = iJacobi * dphidetaj;
        A[1][1] = dphidxj[0];
        A[1][2] = dphidxj[1];

        dphidetaj << 0, 1;
        dphidxj = iJacobi * dphidetaj;
        A[2][1] = dphidxj[0];
        A[2][2] = dphidxj[1];

        dphidetaj << RBF::RBF1x(pc, crbf, RBF::F0, RBF::F1), RBF::RBF1y(pc, crbf, RBF::F0, RBF::F1);
        dphidxj = iJacobi * dphidetaj;
        A[3][1] = dphidxj[0];
        A[3][2] = dphidxj[1];
#ifdef RBFB1_RBF_USE_MEANIJ
        auto meanJ = iJacobi.inverse();
        auto L0 = meanJ(0, Eigen::all).norm();
        auto L1 = meanJ(1, Eigen::all).norm();
        if (L0 > RBFB1_AR0 * L1 || L1 > RBFB1_AR0 * L0)
            A[3][0] *= 0, A[3][1] *= 0, A[3][2] *= 0;
#endif

        return true;
    }

    bool getMomentRBFB1_Interp_15_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 15> &A,
        ScalarCfv::cellFieldData &cell,
        std::vector<ScalarCfv::faceFieldData> &faceVec) // adding 1
    {
        // Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        // Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
        //                                      {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        // Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
        //                                  {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        // Eigen::Vector2d pp = XiNj * Nj;
        // ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        // A[1] = delta.x / scale.x;
        // A[2] = delta.y / scale.y;
        assert(fO == 5);
        if (cell.cellType_ == ScalarCfv::Quadrilateral)
        {

#ifdef RBFB1_RBF_USE_MEANIJ

            Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
            Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                             {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
            Eigen::Matrix2d iJacobi = cell.MeanIJ;
            Eigen::Vector2d pp = XiNj * Nj;
            ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
            Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
            Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
            ScalarCfv::point pc(ppc(0), ppc(1));
#else
            auto pc = p - ScalarCfv::point(0.5, 0.5);
#endif

            A[1] = pc.x;
            A[2] = pc.y;
            // if (cell.boundaryCellType_ == ScalarCfv::boundaryCellType::WallCell)
            //     return false;
            int ic = 3;
            for (int cf = 1; cf <= cell.cellFaceNumber; cf++)
            {
                auto &face = faceVec[cell.cellFaceIndex[cf] - 1];
                for (int gg = 0; gg < 3; gg++)
                {
                    auto &pg = face.gaussPairVector_[gg].p;
                    ScalarCfv::point delta = ScalarCfv::point(pg.x - baryCenter.x, pg.y - baryCenter.y);
                    Eigen::Vector2d deltap = {delta.x, delta.y};
                    Eigen::Vector2d ppcg = iJacobi.transpose() * deltap;
                    ScalarCfv::point pcg(ppcg(0), ppcg(1));
                    A[ic] = RBF::RBF0(pc - pcg, crbf, RBF::F0);
                    ic++;
                }
            }
        }
        else if (cell.cellType_ == ScalarCfv::Triangle)
        {
            assert(false);
        }
        return true;
    }

    bool getBaseValueRBFB1_Interp_15_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 15> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 15> &A,
        ScalarCfv::cellFieldData &cell,
        std::vector<ScalarCfv::faceFieldData> &faceVec) // adding 1
    {
        assert(fO == 5);
        if (cell.cellType_ == ScalarCfv::Quadrilateral)
        {

#ifdef RBFB1_RBF_USE_MEANIJ
            Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
            Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                             {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
            Eigen::Matrix2d iJacobi = cell.MeanIJ;
            Eigen::Vector2d pp = XiNj * Nj;
            ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
            Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
            Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
            ScalarCfv::point pc(ppc(0), ppc(1));
#else
            auto pc = p - ScalarCfv::point(0.5, 0.5);
#endif
            A[1] = pc.x - moment[1];
            A[2] = pc.y - moment[2];
            int ic = 3;
            // if (cell.boundaryCellType_ == ScalarCfv::boundaryCellType::WallCell)
            //     return false;
            for (int cf = 1; cf <= cell.cellFaceNumber; cf++)
            {
                auto &face = faceVec[cell.cellFaceIndex[cf] - 1];
                for (int gg = 0; gg < 3; gg++)
                {
                    auto &pg = face.gaussPairVector_[gg].p;
                    ScalarCfv::point delta = ScalarCfv::point(pg.x - baryCenter.x, pg.y - baryCenter.y);
                    Eigen::Vector2d deltap = {delta.x, delta.y};
                    Eigen::Vector2d ppcg = iJacobi.transpose() * deltap;
                    ScalarCfv::point pcg(ppcg(0), ppcg(1));
                    A[ic] = RBF::RBF0(pc - pcg, crbf, RBF::F0) - moment[ic];
                    ic++;
                }
            }
        }
        else if (cell.cellType_ == ScalarCfv::Triangle)
        {
            assert(false);
        }

        return true;
    }

    bool getDiffBaseValueRBFB1_Interp_15_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 15> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 15, 3> &A,
        ScalarCfv::cellFieldData &cell,
        std::vector<ScalarCfv::faceFieldData> &faceVec) // adding 1
    {
        assert(fO == 5);
        if (cell.cellType_ == ScalarCfv::Quadrilateral)
        {
            //(1-x)(1-y) y(1-x) x*y x*(1-y)
            Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};

            Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                                 {-(1 - p.x), (1 - p.x), p.x, -p.x}};
            Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                             {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};

            Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
            Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj

#ifdef RBFB1_RBF_USE_MEANIJ
            iJacobi = cell.MeanIJ;
            Eigen::Vector2d pp = XiNj * Nj;
            ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
            Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
            Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
            ScalarCfv::point pc(ppc(0), ppc(1));
#else
            auto pc = p - ScalarCfv::point(0.5, 0.5);
#endif
            A[1][0] = pc.x - moment[1];
            A[2][0] = pc.y - moment[2];

            // diff 1:
            Eigen::Vector2d dphidetaj, dphidxj;
            dphidetaj << 1, 0;
            dphidxj = iJacobi * dphidetaj;
            A[1][1] = dphidxj[0];
            A[1][2] = dphidxj[1];

            dphidetaj << 0, 1;
            dphidxj = iJacobi * dphidetaj;
            A[2][1] = dphidxj[0];
            A[2][2] = dphidxj[1];
            // if (cell.boundaryCellType_ == ScalarCfv::boundaryCellType::WallCell)
            //     return false;
            int ic = 3;
            for (int cf = 1; cf <= cell.cellFaceNumber; cf++)
            {
                auto &face = faceVec[cell.cellFaceIndex[cf] - 1];
                for (int gg = 0; gg < 3; gg++)
                {
                    auto &pg = face.gaussPairVector_[gg].p;
                    ScalarCfv::point delta = ScalarCfv::point(pg.x - baryCenter.x, pg.y - baryCenter.y);
                    Eigen::Vector2d deltap = {delta.x, delta.y};
                    Eigen::Vector2d ppcg = iJacobi.transpose() * deltap;
                    ScalarCfv::point pcg(ppcg(0), ppcg(1));
                    dphidetaj << RBF::RBF1x(pc - pcg, crbf, RBF::F0, RBF::F1), RBF::RBF1y(pc - pcg, crbf, RBF::F0, RBF::F1);
                    dphidxj = iJacobi * dphidetaj;
                    A[ic][1] = dphidxj[0];
                    A[ic][2] = dphidxj[1];
                    A[ic][0] = RBF::RBF0(pc - pcg, crbf, RBF::F0) - moment[ic];
                    ic++;
                }
            }
        }
        else if (cell.cellType_ == ScalarCfv::Triangle)
        {
            assert(false);
        }

        return true;
    }

    /*
     */

    /*

    6-6

    */

    bool getMomentRBFB1_POLY(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 0 for rO == 2
    {
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d iJacobi = cell.MeanIJ;
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
        Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
#ifndef RBFB1_GlobalPoly
#ifdef RBFB1_POLY_USE_MEANIJ
        ScalarCfv::point pc(ppc(0), ppc(1));
#else
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
#endif
        A[1] = pc.x;
        A[2] = pc.y;
        A[3] = pc.x * pc.x;
        A[4] = pc.x * pc.y;
        A[5] = pc.y * pc.y;

#else
        A[1] = delta.x / scale.x;
        A[2] = delta.y / scale.y;
        A[3] = delta.x * delta.x / scale.x / scale.x;
        A[4] = delta.x * delta.y / scale.x / scale.y;
        A[5] = delta.y * delta.y / scale.y / scale.y;
#endif

        return true;
    }

    bool getBaseValueRBFB1_POLY(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 6> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 0 for rO = 2
    {
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d iJacobi = cell.MeanIJ;
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
        Eigen::Vector2d ppc = iJacobi.transpose() * deltap;

#ifndef RBFB1_GlobalPoly
#ifdef RBFB1_POLY_USE_MEANIJ
        ScalarCfv::point pc(ppc(0), ppc(1));
#else
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
#endif
        A[1] = pc.x - moment[1];
        A[2] = pc.y - moment[2];
        A[3] = pc.x * pc.x - moment[3];
        A[4] = pc.x * pc.y - moment[4];
        A[5] = pc.y * pc.y - moment[5];
#else

        A[1] = delta.x / scale.x - moment[1];
        A[2] = delta.y / scale.y - moment[2];
        A[3] = delta.x * delta.x / scale.x / scale.x - moment[3];
        A[4] = delta.x * delta.y / scale.x / scale.y - moment[4];
        A[5] = delta.y * delta.y / scale.y / scale.y - moment[5];
#endif

        return true;
    }

    bool getDiffBaseValueRBFB1_POLY(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 6> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 0 for rO = 2
    {

        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
#ifdef RBFB1_POLY_USE_MEANIJ
        iJacobi = cell.MeanIJ;
        Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
        Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
        ScalarCfv::point pc(ppc(0), ppc(1));
#else
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
#endif

#ifndef RBFB1_GlobalPoly
        Eigen::Vector2d dphidetaj, dphidxj;
        Eigen::Matrix2d ddphidetaidetaj, ddphidxidxj;
        Eigen::ETensorR3<ScalarCfv::real, 2, 2, 2> dddphidetaidetajdetak, dddphidxidxj;
        ScalarCfv::real dxxx, dxxy, dxyy, dyyy;

        A[1][0] = pc.x - moment[1];
        A[2][0] = pc.y - moment[2];
        A[3][0] = pc.x * pc.x - moment[3];
        A[4][0] = pc.x * pc.y - moment[4];
        A[5][0] = pc.y * pc.y - moment[5];
        dphidetaj << 1, 0;
        dphidxj = iJacobi * dphidetaj;
        A[1][1] = dphidxj(0);
        A[1][2] = dphidxj(1);
        A[1][3] = A[1][4] = A[1][5] = 0;

        dphidetaj << 0, 1;
        dphidxj = iJacobi * dphidetaj;
        A[2][1] = dphidxj(0);
        A[2][2] = dphidxj(1);
        A[2][3] = A[2][4] = A[2][5] = 0;

        dphidetaj << 2 * pc.x, 0;
        dphidxj = iJacobi * dphidetaj;
        A[3][1] = dphidxj(0);
        A[3][2] = dphidxj(1);
        ddphidetaidetaj << 2, 0, 0, 0;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[3][3] = ddphidxidxj(0, 0);
        A[3][4] = ddphidxidxj(0, 1);
        A[3][5] = ddphidxidxj(1, 1);

        dphidetaj << pc.y, pc.x;
        dphidxj = iJacobi * dphidetaj;
        A[4][1] = dphidxj(0);
        A[4][2] = dphidxj(1);
        ddphidetaidetaj << 0, 1, 1, 0;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[4][3] = ddphidxidxj(0, 0);
        A[4][4] = ddphidxidxj(0, 1);
        A[4][5] = ddphidxidxj(1, 1);

        dphidetaj << 0, 2 * pc.y;
        dphidxj = iJacobi * dphidetaj;
        A[5][1] = dphidxj(0);
        A[5][2] = dphidxj(1);
        ddphidetaidetaj << 0, 0, 0, 2;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[5][3] = ddphidxidxj(0, 0);
        A[5][4] = ddphidxidxj(0, 1);
        A[5][5] = ddphidxidxj(1, 1);
#else
        A[1][0] = delta.x / scale.x - moment[1];
        A[2][0] = delta.y / scale.y - moment[2];
        A[3][0] = delta.x * delta.x / scale.x / scale.x - moment[3];
        A[4][0] = delta.x * delta.y / scale.x / scale.y - moment[4];
        A[5][0] = delta.y * delta.y / scale.y / scale.y - moment[5];

        A[1][1] = 1 / scale.x;
        A[1][2] = 0;
        A[1][3] = A[1][4] = A[1][5] = 0;

        A[2][1] = 0;
        A[2][2] = 1 / scale.y;
        A[2][3] = A[2][4] = A[2][5] = 0;

        A[3][1] = 2 * delta.x / scale.x / scale.x;
        A[4][1] = delta.y / scale.x / scale.y;
        A[5][1] = 0;

        A[3][2] = 0;
        A[4][2] = delta.x / scale.x / scale.y;
        A[5][2] = 2 * delta.y / scale.y / scale.y;

        A[3][3] = 2 / scale.x / scale.x;
        A[4][3] = 0;
        A[5][3] = 0;

        A[3][4] = 0;
        A[4][4] = 1 / scale.x / scale.y;
        A[5][4] = 0;

        A[3][5] = 0;
        A[4][5] = 0;
        A[5][5] = 2 / scale.y / scale.y;

#endif

        return true;
    }

    bool getMomentRBFB1_POLY(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 10> &A,
        ScalarCfv::cellFieldData &cell) // adding 0 for rO == 2
    {
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d iJacobi = cell.MeanIJ;
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
        Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
        ScalarCfv::point pc(ppc(0), ppc(1));
#ifndef RBFB1_GlobalPoly
        // auto pc = p;
        // pc.x -= 0.5;
        // pc.y -= 0.5;
        A[1] = pc.x;
        A[2] = pc.y;
        A[3] = pc.x * pc.x;
        A[4] = pc.x * pc.y;
        A[5] = pc.y * pc.y;
        A[6] = pc.x * pc.x * pc.x;
        A[7] = pc.x * pc.x * pc.y;
        A[8] = pc.x * pc.y * pc.y;
        A[9] = pc.y * pc.y * pc.y;

#else

        A[1] = delta.x / scale.x;
        A[2] = delta.y / scale.y;
        A[3] = delta.x * delta.x / scale.x / scale.x;
        A[4] = delta.x * delta.y / scale.x / scale.y;
        A[5] = delta.y * delta.y / scale.y / scale.y;
        A[6] = delta.x * delta.x * delta.x / scale.x / scale.x / scale.x;
        A[7] = delta.x * delta.x * delta.y / scale.x / scale.x / scale.y;
        A[8] = delta.x * delta.y * delta.y / scale.x / scale.y / scale.y;
        A[9] = delta.y * delta.y * delta.y / scale.y / scale.y / scale.y;
#endif

        return true;
    }

    bool getBaseValueRBFB1_POLY(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 10> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 10> &A,
        ScalarCfv::cellFieldData &cell) // adding 0 for rO = 2
    {
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d iJacobi = cell.MeanIJ;
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
        Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
        ScalarCfv::point pc(ppc(0), ppc(1));
#ifndef RBFB1_GlobalPoly
        // auto pc = p;
        // pc.x -= 0.5;
        // pc.y -= 0.5;
        A[1] = pc.x - moment[1];
        A[2] = pc.y - moment[2];
        A[3] = pc.x * pc.x - moment[3];
        A[4] = pc.x * pc.y - moment[4];
        A[5] = pc.y * pc.y - moment[5];
        A[6] = pc.x * pc.x * pc.x - moment[6];
        A[7] = pc.x * pc.x * pc.y - moment[7];
        A[8] = pc.x * pc.y * pc.y - moment[8];
        A[9] = pc.y * pc.y * pc.y - moment[9];
#else

        A[1] = delta.x / scale.x - moment[1];
        A[2] = delta.y / scale.y - moment[2];
        A[3] = delta.x * delta.x / scale.x / scale.x - moment[3];
        A[4] = delta.x * delta.y / scale.x / scale.y - moment[4];
        A[5] = delta.y * delta.y / scale.y / scale.y - moment[5];
        A[6] = delta.x * delta.x * delta.x / scale.x / scale.x / scale.x - moment[6];
        A[7] = delta.x * delta.x * delta.y / scale.x / scale.x / scale.y - moment[7];
        A[8] = delta.x * delta.y * delta.y / scale.x / scale.y / scale.y - moment[8];
        A[9] = delta.y * delta.y * delta.y / scale.y / scale.y / scale.y - moment[9];
#endif

        return true;
    }

    bool getDiffBaseValueRBFB1_POLY(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 10> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 10, 10> &A,
        ScalarCfv::cellFieldData &cell) // adding 0 for rO = 2
    {
        // auto pc = p;
        // pc.x -= 0.5;
        // pc.y -= 0.5;

        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 3, 4> dNjdetai2{
            {0, 0, 0, 0},
            {1, -1, 1, -1},
            {0, 0, 0, 0}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj
        iJacobi = cell.MeanIJ;
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        Eigen::Vector2d deltap = {pp(0) - baryCenter.x, pp(1) - baryCenter.y};
        Eigen::Vector2d ppc = iJacobi.transpose() * deltap;
        ScalarCfv::point pc(ppc(0), ppc(1));

#ifndef RBFB1_GlobalPoly
        Eigen::Vector2d dphidetaj,
            dphidxj;
        Eigen::Matrix2d ddphidetaidetaj, ddphidxidxj;
        Eigen::ETensorR3<ScalarCfv::real, 2, 2, 2> dddphidetaidetajdetak, dddphidxidxj;
        ScalarCfv::real dxxx, dxxy, dxyy, dyyy;

        A[1][0] = pc.x - moment[1];
        A[2][0] = pc.y - moment[2];
        A[3][0] = pc.x * pc.x - moment[3];
        A[4][0] = pc.x * pc.y - moment[4];
        A[5][0] = pc.y * pc.y - moment[5];
        A[6][0] = pc.x * pc.x * pc.x - moment[6];
        A[7][0] = pc.x * pc.x * pc.y - moment[7];
        A[8][0] = pc.x * pc.y * pc.y - moment[8];
        A[9][0] = pc.y * pc.y * pc.y - moment[9];
        dphidetaj << 1, 0;
        dphidxj = iJacobi * dphidetaj;
        A[1][1] = dphidxj(0);
        A[1][2] = dphidxj(1);
        A[1][3] = A[1][4] = A[1][5] = A[1][6] = A[1][7] = A[1][8] = A[1][9] = 0;

        dphidetaj << 0, 1;
        dphidxj = iJacobi * dphidetaj;
        A[2][1] = dphidxj(0);
        A[2][2] = dphidxj(1);
        A[2][3] = A[2][4] = A[2][5] = A[2][6] = A[2][7] = A[2][8] = A[2][9] = 0;

        dphidetaj << 2 * pc.x, 0;
        dphidxj = iJacobi * dphidetaj;
        A[3][1] = dphidxj(0);
        A[3][2] = dphidxj(1);
        ddphidetaidetaj << 2, 0, 0, 0;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[3][3] = ddphidxidxj(0, 0);
        A[3][4] = ddphidxidxj(0, 1);
        A[3][5] = ddphidxidxj(1, 1);
        A[3][6] = A[3][7] = A[3][8] = A[3][9] = 0;

        dphidetaj << pc.y, pc.x;
        dphidxj = iJacobi * dphidetaj;
        A[4][1] = dphidxj(0);
        A[4][2] = dphidxj(1);
        ddphidetaidetaj << 0, 1, 1, 0;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[4][3] = ddphidxidxj(0, 0);
        A[4][4] = ddphidxidxj(0, 1);
        A[4][5] = ddphidxidxj(1, 1);
        A[4][6] = A[4][7] = A[4][8] = A[4][9] = 0;

        dphidetaj << 0, 2 * pc.y;
        dphidxj = iJacobi * dphidetaj;
        A[5][1] = dphidxj(0);
        A[5][2] = dphidxj(1);
        ddphidetaidetaj << 0, 0, 0, 2;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[5][3] = ddphidxidxj(0, 0);
        A[5][4] = ddphidxidxj(0, 1);
        A[5][5] = ddphidxidxj(1, 1);
        A[5][6] = A[5][7] = A[5][8] = A[5][9] = 0;

        dphidetaj << 3 * pc.x * pc.x, 0;
        dphidxj = iJacobi * dphidetaj;
        A[6][1] = dphidxj(0);
        A[6][2] = dphidxj(1);
        ddphidetaidetaj << 6 * pc.x, 0, 0, 0;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[6][3] = ddphidxidxj(0, 0);
        A[6][4] = ddphidxidxj(0, 1);
        A[6][5] = ddphidxidxj(1, 1);
        dxxx = 6, dxxy = 0, dxyy = 0, dyyy = 0;
        dddphidetaidetajdetak(0, 0, 0) = dxxx;
        dddphidetaidetajdetak(0, 0, 1) = dddphidetaidetajdetak(0, 1, 0) = dddphidetaidetajdetak(1, 0, 0) = dxxy;
        dddphidetaidetajdetak(0, 1, 1) = dddphidetaidetajdetak(1, 1, 0) = dddphidetaidetajdetak(1, 0, 1) = dxyy;
        dddphidetaidetajdetak(1, 1, 1) = dyyy;
        dddphidetaidetajdetak.MatTransform0(iJacobi.transpose());
        dddphidetaidetajdetak.MatTransform1(iJacobi.transpose());
        dddphidetaidetajdetak.MatTransform2(iJacobi.transpose());
        A[6][6] = dddphidetaidetajdetak(0, 0, 0);
        A[6][7] = dddphidetaidetajdetak(0, 0, 1);
        A[6][8] = dddphidetaidetajdetak(0, 1, 1);
        A[6][9] = dddphidetaidetajdetak(1, 1, 1);

        dphidetaj << 2 * pc.x * pc.y, pc.x * pc.x;
        dphidxj = iJacobi * dphidetaj;
        A[7][1] = dphidxj(0);
        A[7][2] = dphidxj(1);
        ddphidetaidetaj << 2 * pc.y, 2 * pc.x, 2 * pc.x, 0;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[7][3] = ddphidxidxj(0, 0);
        A[7][4] = ddphidxidxj(0, 1);
        A[7][5] = ddphidxidxj(1, 1);
        dxxx = 0, dxxy = 2, dxyy = 0, dyyy = 0;
        dddphidetaidetajdetak(0, 0, 0) = dxxx;
        dddphidetaidetajdetak(0, 0, 1) = dddphidetaidetajdetak(0, 1, 0) = dddphidetaidetajdetak(1, 0, 0) = dxxy;
        dddphidetaidetajdetak(0, 1, 1) = dddphidetaidetajdetak(1, 1, 0) = dddphidetaidetajdetak(1, 0, 1) = dxyy;
        dddphidetaidetajdetak(1, 1, 1) = dyyy;
        dddphidetaidetajdetak.MatTransform0(iJacobi.transpose());
        dddphidetaidetajdetak.MatTransform1(iJacobi.transpose());
        dddphidetaidetajdetak.MatTransform2(iJacobi.transpose());
        A[7][6] = dddphidetaidetajdetak(0, 0, 0);
        A[7][7] = dddphidetaidetajdetak(0, 0, 1);
        A[7][8] = dddphidetaidetajdetak(0, 1, 1);
        A[7][9] = dddphidetaidetajdetak(1, 1, 1);

        dphidetaj << pc.y * pc.y, 2 * pc.x * pc.y;
        dphidxj = iJacobi * dphidetaj;
        A[8][1] = dphidxj(0);
        A[8][2] = dphidxj(1);
        ddphidetaidetaj << 0, 2 * pc.y, 2 * pc.y, 2 * pc.x;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[8][3] = ddphidxidxj(0, 0);
        A[8][4] = ddphidxidxj(0, 1);
        A[8][5] = ddphidxidxj(1, 1);
        dxxx = 0, dxxy = 0, dxyy = 2, dyyy = 0;
        dddphidetaidetajdetak(0, 0, 0) = dxxx;
        dddphidetaidetajdetak(0, 0, 1) = dddphidetaidetajdetak(0, 1, 0) = dddphidetaidetajdetak(1, 0, 0) = dxxy;
        dddphidetaidetajdetak(0, 1, 1) = dddphidetaidetajdetak(1, 1, 0) = dddphidetaidetajdetak(1, 0, 1) = dxyy;
        dddphidetaidetajdetak(1, 1, 1) = dyyy;
        dddphidetaidetajdetak.MatTransform0(iJacobi.transpose());
        dddphidetaidetajdetak.MatTransform1(iJacobi.transpose());
        dddphidetaidetajdetak.MatTransform2(iJacobi.transpose());
        A[8][6] = dddphidetaidetajdetak(0, 0, 0);
        A[8][7] = dddphidetaidetajdetak(0, 0, 1);
        A[8][8] = dddphidetaidetajdetak(0, 1, 1);
        A[8][9] = dddphidetaidetajdetak(1, 1, 1);

        dphidetaj << 0, 3 * pc.y * pc.y;
        dphidxj = iJacobi * dphidetaj;
        A[9][1] = dphidxj(0);
        A[9][2] = dphidxj(1);
        ddphidetaidetaj << 0, 0, 0, 6 * pc.y;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[9][3] = ddphidxidxj(0, 0);
        A[9][4] = ddphidxidxj(0, 1);
        A[9][5] = ddphidxidxj(1, 1);
        dxxx = 0, dxxy = 0, dxyy = 0, dyyy = 6;
        dddphidetaidetajdetak(0, 0, 0) = dxxx;
        dddphidetaidetajdetak(0, 0, 1) = dddphidetaidetajdetak(0, 1, 0) = dddphidetaidetajdetak(1, 0, 0) = dxxy;
        dddphidetaidetajdetak(0, 1, 1) = dddphidetaidetajdetak(1, 1, 0) = dddphidetaidetajdetak(1, 0, 1) = dxyy;
        dddphidetaidetajdetak(1, 1, 1) = dyyy;
        dddphidetaidetajdetak.MatTransform0(iJacobi.transpose());
        dddphidetaidetajdetak.MatTransform1(iJacobi.transpose());
        dddphidetaidetajdetak.MatTransform2(iJacobi.transpose());
        A[9][6] = dddphidetaidetajdetak(0, 0, 0);
        A[9][7] = dddphidetaidetajdetak(0, 0, 1);
        A[9][8] = dddphidetaidetajdetak(0, 1, 1);
        A[9][9] = dddphidetaidetajdetak(1, 1, 1);
#else

        A[1][0] = delta.x / scale.x - moment[1];
        A[2][0] = delta.y / scale.y - moment[2];
        A[3][0] = delta.x * delta.x / scale.x / scale.x - moment[3];
        A[4][0] = delta.x * delta.y / scale.x / scale.y - moment[4];
        A[5][0] = delta.y * delta.y / scale.y / scale.y - moment[5];
        A[6][0] = delta.x * delta.x * delta.x / scale.x / scale.x / scale.x - moment[6];
        A[7][0] = delta.x * delta.x * delta.y / scale.x / scale.x / scale.y - moment[7];
        A[8][0] = delta.x * delta.y * delta.y / scale.x / scale.y / scale.y - moment[8];
        A[9][0] = delta.y * delta.y * delta.y / scale.y / scale.y / scale.y - moment[9];

        A[1][1] = 1 / scale.x;
        A[1][2] = 0;
        A[1][3] = A[1][4] = A[1][5] = A[1][6] = A[1][7] = A[1][8] = A[1][9] = 0;

        A[2][1] = 0;
        A[2][2] = 1 / scale.y;
        A[2][3] = A[2][4] = A[2][5] = A[2][6] = A[2][7] = A[2][8] = A[2][9] = 0;

        A[3][1] = 2 * delta.x / scale.x / scale.x;
        A[4][1] = delta.y / scale.x / scale.y;
        A[5][1] = 0;
        A[6][1] = 3 * delta.x * delta.x / scale.x / scale.x / scale.x;
        A[7][1] = 2 * delta.x * delta.y / scale.x / scale.x / scale.y;
        A[8][1] = delta.y * delta.y / scale.x / scale.y / scale.y;
        A[9][1] = 0;

        A[3][2] = 0;
        A[4][2] = delta.x / scale.x / scale.y;
        A[5][2] = 2 * delta.y / scale.y / scale.y;
        A[6][2] = 0;
        A[7][2] = delta.x * delta.x / scale.x / scale.x / scale.y;
        A[8][2] = 2 * delta.x * delta.y / scale.x / scale.y / scale.y;
        A[9][2] = 3 * delta.y * delta.y / scale.y / scale.y / scale.y;

        A[3][3] = 2 / scale.x / scale.x;
        A[4][3] = 0;
        A[5][3] = 0;
        A[6][3] = 6 * delta.x / scale.x / scale.x / scale.x;
        A[7][3] = 2 * delta.y / scale.x / scale.x / scale.y;
        A[8][3] = 0;
        A[9][3] = 0;

        A[3][4] = 0;
        A[4][4] = 1 / scale.x / scale.y;
        A[5][4] = 0;
        A[6][4] = 0;
        A[7][4] = 2 * delta.x / scale.x / scale.x / scale.y;
        A[8][4] = 2 * delta.y / scale.x / scale.y / scale.y;
        A[9][4] = 0;

        A[3][5] = 0;
        A[4][5] = 0;
        A[5][5] = 2 / scale.y / scale.y;
        A[6][5] = 0;
        A[7][5] = 0;
        A[8][5] = 2 * delta.x / scale.x / scale.y / scale.y;
        A[9][5] = 6 * delta.y / scale.y / scale.y / scale.y;

        A[3][6] = A[3][7] = A[3][8] = A[3][9] = 0;
        A[4][6] = A[4][7] = A[4][8] = A[4][9] = 0;
        A[5][6] = A[5][7] = A[5][8] = A[5][9] = 0;

        A[6][6] = 6 / scale.x / scale.x / scale.x;
        A[7][6] = 0;
        A[8][6] = 0;
        A[9][6] = 0;

        A[6][7] = 0;
        A[7][7] = 2 / scale.x / scale.x / scale.y;
        A[8][7] = 0;
        A[9][7] = 0;

        A[6][8] = 0;
        A[7][8] = 0;
        A[8][8] = 2 / scale.x / scale.y / scale.y;
        A[9][8] = 0;

        A[6][9] = 0;
        A[7][9] = 0;
        A[8][9] = 0;
        A[9][9] = 6 / scale.y / scale.y / scale.y;

#endif

        return true;
    }

    ScalarCfv::real B1_QuadP3[3][2] = {{0.2, 0.2}, {0.8, 0.2}, {0.5, 0.8}};

    bool getMomentRBFB1(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 3 for rO == 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
        A[1] = pc.x;
        A[2] = pc.y;
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 3; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP3[i][0], B1_QuadP3[i][1]);
            A[3 + i] = RBF::RBF0(pcc, crbf, RBF::F0);
        }

        return true;
    }

    bool getBaseValueRBFB1(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 6> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 0 for rO = 2
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
        A[1] = pc.x - moment[1];
        A[2] = pc.y - moment[2];
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 3; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP3[i][0], B1_QuadP3[i][1]);
            A[3 + i] = RBF::RBF0(pcc, crbf, RBF::F0) - moment[3 + i];
        }

        return true;
    }

    bool getDiffBaseValueRBFB1(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 6> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 6, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 0 for rO = 2
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
        A[1][0] = pc.x - moment[1];
        A[2][0] = pc.y - moment[2];
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 3; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP3[i][0], B1_QuadP3[i][1]);
            A[3 + i][0] = RBF::RBF0(pcc, crbf, RBF::F0) - moment[3 + i];
        }

        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj

        Eigen::Vector2d dphidetaj, dphidxj;
        Eigen::Matrix2d ddphidetaidetaj, ddphidxidxj;
        dphidetaj << 1, 0;
        dphidxj = iJacobi * dphidetaj;
        A[1][1] = dphidxj(0);
        A[1][2] = dphidxj(1);
        A[1][3] = 0;
        A[1][4] = 0;
        A[1][5] = 0;

        dphidetaj << 0, 1;
        dphidxj = iJacobi * dphidetaj;
        A[2][1] = dphidxj(0);
        A[2][2] = dphidxj(1);
        A[2][3] = 0;
        A[2][4] = 0;
        A[2][5] = 0;

        for (int i = 0; i < 3; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP3[i][0], B1_QuadP3[i][1]);
            A[3 + i][0] = 0 * RBF::RBF0(pcc, crbf, RBF::F0) - moment[3 + i];

            Eigen::Vector2d dphidetaj, dphidxj;
            dphidetaj << RBF::RBF1x(pcc, crbf, RBF::F0, RBF::F1), RBF::RBF1y(pcc, crbf, RBF::F0, RBF::F1);
            dphidxj = iJacobi * dphidetaj;
            A[3 + i][1] = dphidxj[0];
            A[3 + i][2] = dphidxj[1];

            Eigen::Matrix2d dphidetaidetaj, dphidxidxj;
            ScalarCfv::real rbfxx, rbfxy, rbfyy;
            rbfxx = RBF::RBF2xx(pcc, crbf, RBF::F0, RBF::F1, RBF::F2);
            rbfxy = RBF::RBF2xy(pcc, crbf, RBF::F0, RBF::F1, RBF::F2);
            rbfyy = RBF::RBF2yy(pcc, crbf, RBF::F0, RBF::F1, RBF::F2);
            dphidetaidetaj << rbfxx, rbfxy, rbfxy, rbfyy;
            dphidxidxj = iJacobi * dphidetaidetaj * iJacobi.transpose();
            A[3 + i][3] = dphidxidxj(0, 0);
            A[3 + i][4] = dphidxidxj(0, 1);
            A[3 + i][5] = dphidxidxj(1, 1);
        }

        return true;
    }

    // ScalarCfv::real B1_QuadP4[4][2] = {{0.45, 0.45}, {0.55, 0.45}, {0.45, 0.55}, {0.55, 0.55}};
    // ScalarCfv::real B1_QuadP4[4][2] = {{0.2, 0.2}, {0.2, 0.8}, {0.8, 0.2}, {0.8, 0.8}};
    ScalarCfv::real B1_QuadP4[4][2] = {{0.4, 0.4}, {0.4, 0.6}, {0.6, 0.4}, {0.6, 0.6}};
    /*

    7-6

    */

    bool getMomentRBFB1_7_6(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO == 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
#ifndef RBFB1_GlobalPoly
        A[1] = pc.x;
        A[2] = pc.y;
#else
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        A[1] = delta.x / scale.x;
        A[2] = delta.y / scale.y;

#endif

        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[3 + i] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0);
        }

        return true;
    }

    bool getBaseValueRBFB1_7_6(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO = 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
#ifndef RBFB1_GlobalPoly
        A[1] = pc.x - moment[1];
        A[2] = pc.y - moment[2];
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
#else
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        A[1] = delta.x / scale.x - moment[1];
        A[2] = delta.y / scale.y - moment[2];

#endif

        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[3 + i] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0) - moment[3 + i];
        }

        return true;
    }

    bool getDiffBaseValueRBFB1_7_6(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 7, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO = 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;

        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[3 + i][0] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0) - moment[3 + i];
        }

        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj

        Eigen::Vector2d dphidetaj, dphidxj;
        Eigen::Matrix2d ddphidetaidetaj, ddphidxidxj;

#ifndef RBFB1_GlobalPoly
        A[1][0] = pc.x - moment[1];
        A[2][0] = pc.y - moment[2];
        dphidetaj << 1, 0;
        dphidxj = iJacobi * dphidetaj;
        A[1][1] = dphidxj(0);
        A[1][2] = dphidxj(1);
        A[1][3] = 0;
        A[1][4] = 0;
        A[1][5] = 0;

        dphidetaj << 0, 1;
        dphidxj = iJacobi * dphidetaj;
        A[2][1] = dphidxj(0);
        A[2][2] = dphidxj(1);
        A[2][3] = 0;
        A[2][4] = 0;
        A[2][5] = 0;
#else
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        A[1][0] = delta.x / scale.x - moment[1];
        A[1][1] = 1 / scale.x;
        A[1][2] = 0;
        A[2][0] = delta.y / scale.y - moment[2];
        A[2][1] = 0;
        A[2][2] = 1 / scale.y;

        A[1][3] = 0;
        A[1][4] = 0;
        A[1][5] = 0;
        A[2][3] = 0;
        A[2][4] = 0;
        A[2][5] = 0;

#endif

        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[3 + i][0] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0) - moment[3 + i];

            Eigen::Vector2d dphidetaj, dphidxj;
            dphidetaj << RBF::RBF1x(pcc, crbf, RBF::F0, RBF::F1), RBF::RBF1y(pcc, crbf, RBF::F0, RBF::F1);
            dphidxj = iJacobi * dphidetaj;
            A[3 + i][1] = rbfScale * dphidxj[0];
            A[3 + i][2] = rbfScale * dphidxj[1];

            Eigen::Matrix2d dphidetaidetaj, dphidxidxj;
            ScalarCfv::real rbfxx, rbfxy, rbfyy;
            rbfxx = RBF::RBF2xx(pcc, crbf, RBF::F0, RBF::F1, RBF::F2);
            rbfxy = RBF::RBF2xy(pcc, crbf, RBF::F0, RBF::F1, RBF::F2);
            rbfyy = RBF::RBF2yy(pcc, crbf, RBF::F0, RBF::F1, RBF::F2);
            dphidetaidetaj << rbfxx, rbfxy, rbfxy, rbfyy;
            dphidxidxj = iJacobi * dphidetaidetaj * iJacobi.transpose();
            A[3 + i][3] = rbfScale * dphidxidxj(0, 0);
            A[3 + i][4] = rbfScale * dphidxidxj(0, 1);
            A[3 + i][5] = rbfScale * dphidxidxj(1, 1);
        }

        return true;
    }

    /*

    7-3

    */
    bool getMomentRBFB1_7_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO == 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
        A[1] = pc.x;
        A[2] = pc.y;
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[3 + i] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0);
        }

        return true;
    }

    bool getBaseValueRBFB1_7_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO = 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
        A[1] = pc.x - moment[1];
        A[2] = pc.y - moment[2];
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[3 + i] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0) - moment[3 + i];
        }

        return true;
    }

    bool getDiffBaseValueRBFB1_7_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 7, 3> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO = 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
        A[1][0] = pc.x - moment[1];
        A[2][0] = pc.y - moment[2];
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[3 + i][0] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0) - moment[3 + i];
        }

        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj

        Eigen::Vector2d dphidetaj, dphidxj;
        Eigen::Matrix2d ddphidetaidetaj, ddphidxidxj;
        dphidetaj << 1, 0;
        dphidxj = iJacobi * dphidetaj;
        A[1][1] = dphidxj(0);
        A[1][2] = dphidxj(1);

        dphidetaj << 0, 1;
        dphidxj = iJacobi * dphidetaj;
        A[2][1] = dphidxj(0);
        A[2][2] = dphidxj(1);

        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[3 + i][0] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0) - moment[3 + i];

            Eigen::Vector2d dphidetaj, dphidxj;
            dphidetaj << RBF::RBF1x(pcc, crbf, RBF::F0, RBF::F1), RBF::RBF1y(pcc, crbf, RBF::F0, RBF::F1);
            dphidxj = iJacobi * dphidetaj;
            A[3 + i][1] = rbfScale * dphidxj[0];
            A[3 + i][2] = rbfScale * dphidxj[1];
        }

        return true;
    }

    bool getMomentRBFB1_7_6_P2(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO == 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
#ifndef RBFB1_GlobalPoly
        A[1] = pc.x;
        A[2] = pc.y;
        A[3] = pc.x * pc.x;
        A[4] = pc.x * pc.y;
        A[5] = pc.y * pc.y;
#else
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        A[1] = delta.x / scale.x;
        A[2] = delta.y / scale.y;
        A[3] = delta.x * delta.x / scale.x / scale.x;
        A[4] = delta.x * delta.y / scale.x / scale.y;
        A[5] = delta.y * delta.y / scale.y / scale.y;

#endif

        A[6] = rbfScale * RBF::RBF0(pc, crbf, RBF::F0);
        return true;
    }

    bool getBaseValueRBFB1_7_6_P2(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO = 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
#ifndef RBFB1_GlobalPoly
        A[1] = pc.x - moment[1];
        A[2] = pc.y - moment[2];
        A[3] = pc.x * pc.x - moment[3];
        A[4] = pc.x * pc.y - moment[4];
        A[5] = pc.y * pc.y - moment[5];
#else
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        A[1] = delta.x / scale.x - moment[1];
        A[2] = delta.y / scale.y - moment[2];
        A[3] = delta.x * delta.x / scale.x / scale.x - moment[3];
        A[4] = delta.x * delta.y / scale.x / scale.y - moment[4];
        A[5] = delta.y * delta.y / scale.y / scale.y - moment[5];

#endif

        A[6] = rbfScale * RBF::RBF0(pc, crbf, RBF::F0);
        return true;
    }

    bool getDiffBaseValueRBFB1_7_6_P2(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 7> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 7, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO = 1
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;

        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj

        Eigen::Vector2d dphidetaj, dphidxj;
        Eigen::Matrix2d ddphidetaidetaj, ddphidxidxj;

#ifndef RBFB1_GlobalPoly
        A[1][0] = pc.x - moment[1];
        A[2][0] = pc.y - moment[2];
        A[3][0] = pc.x * pc.x - moment[3];
        A[4][0] = pc.x * pc.y - moment[4];
        A[5][0] = pc.y * pc.y - moment[5];

        dphidetaj << 1, 0;
        dphidxj = iJacobi * dphidetaj;
        A[1][1] = dphidxj(0);
        A[1][2] = dphidxj(1);
        A[1][3] = 0;
        A[1][4] = 0;
        A[1][5] = 0;

        dphidetaj << 0, 1;
        dphidxj = iJacobi * dphidetaj;
        A[2][1] = dphidxj(0);
        A[2][2] = dphidxj(1);
        A[2][3] = 0;
        A[2][4] = 0;
        A[2][5] = 0;

        dphidetaj << 2 * pc.x, 0;
        dphidxj = iJacobi * dphidetaj;
        A[3][1] = dphidxj(0);
        A[3][2] = dphidxj(1);
        ddphidetaidetaj << 2, 0, 0, 0;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[3][3] = ddphidxidxj(0, 0);
        A[3][4] = ddphidxidxj(0, 1);
        A[3][5] = ddphidxidxj(1, 1);

        dphidetaj << pc.y, pc.x;
        dphidxj = iJacobi * dphidetaj;
        A[4][1] = dphidxj(0);
        A[4][2] = dphidxj(1);
        ddphidetaidetaj << 0, 1, 1, 0;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[4][3] = ddphidxidxj(0, 0);
        A[4][4] = ddphidxidxj(0, 1);
        A[4][5] = ddphidxidxj(1, 1);

        dphidetaj << 0, 2 * pc.y;
        dphidxj = iJacobi * dphidetaj;
        A[5][1] = dphidxj(0);
        A[5][2] = dphidxj(1);
        ddphidetaidetaj << 0, 0, 0, 2;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[5][3] = ddphidxidxj(0, 0);
        A[5][4] = ddphidxidxj(0, 1);
        A[5][5] = ddphidxidxj(1, 1);

#else
        Eigen::Vector2d pp = XiNj * Nj;
        ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        A[1][0] = delta.x / scale.x - moment[1];
        A[1][1] = 1 / scale.x;
        A[1][2] = 0;
        A[2][0] = delta.y / scale.y - moment[2];
        A[2][1] = 0;
        A[2][2] = 1 / scale.y;

        A[1][3] = 0;
        A[1][4] = 0;
        A[1][5] = 0;
        A[2][3] = 0;
        A[2][4] = 0;
        A[2][5] = 0;

        A[3][0] = delta.x * delta.x / scale.x / scale.x - moment[3];
        A[4][0] = delta.x * delta.y / scale.x / scale.y - moment[4];
        A[5][0] = delta.y * delta.y / scale.y / scale.y - moment[5];

        A[3][1] = 2 * delta.x / scale.x / scale.x;
        A[4][1] = delta.y / scale.x / scale.y;
        A[5][1] = 0;

        A[3][2] = 0;
        A[4][2] = delta.x / scale.x / scale.y;
        A[5][2] = 2 * delta.y / scale.y / scale.y;

        A[3][3] = 2 / scale.x / scale.x;
        A[4][3] = 0;
        A[5][3] = 0;

        A[3][4] = 0;
        A[4][4] = 1 / scale.x / scale.y;
        A[5][4] = 0;

        A[3][5] = 0;
        A[4][5] = 0;
        A[5][5] = 2 / scale.y / scale.y;

#endif

        A[6][0] = rbfScale * RBF::RBF0(pc, crbf, RBF::F0) - moment[6];

        dphidetaj << RBF::RBF1x(pc, crbf, RBF::F0, RBF::F1), RBF::RBF1y(pc, crbf, RBF::F0, RBF::F1);
        dphidxj = iJacobi * dphidetaj;
        A[6][1] = rbfScale * dphidxj[0];
        A[6][2] = rbfScale * dphidxj[1];

        ScalarCfv::real rbfxx, rbfxy, rbfyy;
        rbfxx = RBF::RBF2xx(pc, crbf, RBF::F0, RBF::F1, RBF::F2);
        rbfxy = RBF::RBF2xy(pc, crbf, RBF::F0, RBF::F1, RBF::F2);
        rbfyy = RBF::RBF2yy(pc, crbf, RBF::F0, RBF::F1, RBF::F2);
        ddphidetaidetaj << rbfxx, rbfxy, rbfxy, rbfyy;
        ddphidxidxj = iJacobi * ddphidetaidetaj * iJacobi.transpose();
        A[6][3] = rbfScale * ddphidxidxj(0, 0);
        A[6][4] = rbfScale * ddphidxidxj(0, 1);
        A[6][5] = rbfScale * ddphidxidxj(1, 1);

        return true;
    }

    bool getMomentRBFB1_5_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 5> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO == 1
    {
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[1 + i] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0);
        }

        return true;
    }

    bool getBaseValueRBFB1_5_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 5> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 5> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO = 1
    {
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[1 + i] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0) - moment[1 + i];
        }

        return true;
    }

    bool getDiffBaseValueRBFB1_5_3(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 5> &moment,
        ScalarCfv::tensor2D<ScalarCfv::real, 5, 3> &A,
        ScalarCfv::cellFieldData &cell) // adding 4 for rO = 1
    {
        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
        Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
                                             {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
                                         {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        Eigen::Matrix2d Jacobi = dNjdetai * XiNj.transpose(); // = dxj/detai
        Eigen::Matrix2d iJacobi = Jacobi.inverse();           // = detaj/dxi // iJabcobi * dphideetaj = dphidxj

        for (int i = 0; i < 4; i++)
        {
            auto pcc = p - ScalarCfv::point(B1_QuadP4[i][0], B1_QuadP4[i][1]);
            A[1 + i][0] = rbfScale * RBF::RBF0(pcc, crbf, RBF::F0) - moment[1 + i];

            Eigen::Vector2d dphidetaj, dphidxj;
            dphidetaj << RBF::RBF1x(pcc, crbf, RBF::F0, RBF::F1), RBF::RBF1y(pcc, crbf, RBF::F0, RBF::F1);
            dphidxj = iJacobi * dphidetaj;
            A[1 + i][1] = rbfScale * dphidxj[0];
            A[1 + i][2] = rbfScale * dphidxj[1];
        }

        return true;
    }
}
