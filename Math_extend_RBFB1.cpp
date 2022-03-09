#pragma once

#include "Math.h"

#include "RBF.h"

#include "Point.h"

namespace CfvMath
{
    ScalarCfv::real rbfScale = 1.0;
    /*

    3-3

    */

    bool getMomentRBFB1(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 3> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {
        // Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        // Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
        // 									 {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        // Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
        // 								 {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        // Eigen::Vector2d pp = XiNj * Nj;
        // ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        // A[1] = delta.x / scale.x;
        // A[2] = delta.y / scale.y;

        A[1] = p.x - 0.5;
        A[2] = p.y - 0.5;
#ifdef TRIAL
        // std::cout << "Moment in " << p.x << p.y << std::endl;
#endif
        return true;
    }

    bool getBaseValueRBFB1(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 3> &moment,
        ScalarCfv::tensor1D<ScalarCfv::real, 3> &A,
        ScalarCfv::cellFieldData &cell) // adding 1
    {
        // Eigen::Vector4d Nj{(1 - p.x) * (1 - p.y), p.y * (1 - p.x), p.x * p.y, p.x * (1 - p.y)};
        // Eigen::Matrix<double, 2, 4> dNjdetai{{-(1 - p.y), -p.y, p.y, (1 - p.y)},
        // 									 {-(1 - p.x), (1 - p.x), p.x, -p.x}};
        // Eigen::Matrix<double, 2, 4> XiNj{{cell.cellNode[1].second.x, cell.cellNode[2].second.x, cell.cellNode[3].second.x, cell.cellNode[4].second.x},
        // 								 {cell.cellNode[1].second.y, cell.cellNode[2].second.y, cell.cellNode[3].second.y, cell.cellNode[4].second.y}};
        // Eigen::Vector2d pp = XiNj * Nj;
        // ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        // A[1] = delta.x / scale.x - moment[1];
        // A[2] = delta.y / scale.y - moment[2];

        A[1] = p.x - 0.5 - moment[1];
        A[2] = p.y - 0.5 - moment[2];
        return true;
    }

    bool getDiffBaseValueRBFB1(
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
        // Eigen::Vector2d pp = XiNj * Nj;
        // ScalarCfv::point delta = ScalarCfv::point(pp(0) - baryCenter.x, pp(1) - baryCenter.y);
        // A[1][0] = delta.x / scale.x - moment[1];
        // A[2][0] = delta.y / scale.y - moment[2];

        // A[1][1] = 1 / scale.x;
        // A[1][2] = 0 / scale.x;
        // A[2][1] = 0 / scale.y;
        // A[2][2] = 1 / scale.y;

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
        if (abs(baryCenter.x - 0.552) < 1e-2 && abs(baryCenter.y - 3.056) < 1e-1)
            return false;
        if (abs(baryCenter.x - 0.0833) < 1e-2 && abs(baryCenter.y - 3.038) < 1e-1)
            return false;

        return true;
    }

    /*

    4-6

    */

    bool getMomentRBFB1(
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

    bool getBaseValueRBFB1(
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

    bool getDiffBaseValueRBFB1(
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

    /*

    6-6

    */
#ifdef B1_6_6_Poly
    bool getMomentRBFB1(
        ScalarCfv::point p,          // parametric place
        ScalarCfv::point baryCenter, // dummy
        ScalarCfv::point scale,      // dummy
        ScalarCfv::tensor1D<ScalarCfv::real, 6> &A,
        ScalarCfv::cellFieldData &cell) // adding 0 for rO == 2
    {
        auto pc = p;
        pc.x -= 0.5;
        pc.y -= 0.5;
        A[1] = pc.x;
        A[2] = pc.y;
        A[3] = pc.x * pc.x;
        A[4] = pc.x * pc.y;
        A[5] = pc.y * pc.y;

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
        A[3] = pc.x * pc.x - moment[3];
        A[4] = pc.x * pc.y - moment[4];
        A[5] = pc.y * pc.y - moment[5];

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
        A[3][0] = pc.x * pc.x - moment[3];
        A[4][0] = pc.x * pc.y - moment[4];
        A[5][0] = pc.y * pc.y - moment[5];

        assert(cell.cellType_ == ScalarCfv::Quadrilateral);
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

        return true;
    }
#else

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
#endif

    ScalarCfv::real B1_QuadP4[4][2] = {{0.45, 0.45}, {0.55, 0.45}, {0.45, 0.55}, {0.55, 0.55}};
    /*

    7-6

    */

#if rO == 1

#define RBFB1P4_Diff2

#ifdef RBFB1P4_Diff2

    bool getMomentRBFB1(
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

    bool getBaseValueRBFB1(
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

    bool getDiffBaseValueRBFB1(
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

#else
    /*

    7-3

    */
    bool getMomentRBFB1(
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

    bool getBaseValueRBFB1(
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

    bool getDiffBaseValueRBFB1(
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

#endif

#elif rO == 2
    bool getMomentRBFB1(
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

    bool getBaseValueRBFB1(
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

    bool getDiffBaseValueRBFB1(
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
#endif
}