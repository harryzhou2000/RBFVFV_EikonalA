#ifndef _RBF_H
#define _RBF_H

#include "TypeDefine.h"
#include "Point.h"
#include <functional>

namespace RBF
{
    const ScalarCfv::real eps = 1e-7;

    ScalarCfv::real Gaussian0(ScalarCfv::real rr, ScalarCfv::real c)
    {
        return std::exp(-rr * rr / (c * c));
    }

    ScalarCfv::real Gaussian1(ScalarCfv::real rr, ScalarCfv::real c)
    {
        auto c2 = c * c;
        return -2 * rr / c2 * std::exp(-rr * rr / c2);
    }

    ScalarCfv::real Gaussian2(ScalarCfv::real rr, ScalarCfv::real c)
    {
        auto rr2 = rr * rr;
        auto c2 = c * c;
        auto f0 = exp(-rr2 / c2);
        return (4 * rr2 * f0) / (c2 * c2) - (2 * f0) / c2;
    }

    ScalarCfv::real MQ0(ScalarCfv::real rr, ScalarCfv::real c)
    {
        auto rdc = rr / c;
        return std::sqrt(1 + rdc * rdc);
    }

    ScalarCfv::real MQ1(ScalarCfv::real rr, ScalarCfv::real c)
    {
        auto c2 = c * c;
        auto rr2 = rr * rr;
        return rr / (c2 * std::sqrt(rr2 / c2 + 1));
    }

    ScalarCfv::real MQ2(ScalarCfv::real rr, ScalarCfv::real c)
    {
        auto rr2 = rr * rr;
        auto c2 = c * c;
        return 1.0 / (c2 * std::pow(rr2 / c2 + 1, 1.5));
    }

    ScalarCfv::real PHSpline3P0(ScalarCfv::real rr, ScalarCfv::real c)
    {
        auto rdc = rr / c;
        if (rdc <= eps)
            return 0.0;
        return std::log(rdc) * std::pow(rdc, 4.0);
    }

    ScalarCfv::real PHSpline3P1(ScalarCfv::real rr, ScalarCfv::real c)
    {
        auto rdc = rr / c;
        if (rdc <= eps)
            return 0.0;
        return std::pow(rdc, 3.0) / c * (1 + 4 * std::log(rdc));
    }

    ScalarCfv::real PHSpline3P2(ScalarCfv::real rr, ScalarCfv::real c)
    {
        auto rdc = rr / c;
        if (rdc <= eps)
            return 0.0;
        return (rdc * rdc) / (c * c) * (7 + 12 * std::log(rdc));
    }

    ScalarCfv::real PHSpline3P3(ScalarCfv::real rr, ScalarCfv::real c)
    {
        auto rdc = rr / c;
        if (rdc <= eps)
            return 0.0;
        return (rdc) / (c * c * c) * (26 + 24 * std::log(rdc));
    }

    // auto F0 = Gaussian0;
    // auto F1 = Gaussian1;
    // auto F2 = Gaussian2;
    // auto F0 = MQ0;
    // auto F1 = MQ1;
    // auto F2 = MQ2;
    auto F0 = PHSpline3P0;
    auto F1 = PHSpline3P1;
    auto F2 = PHSpline3P2;
    auto F3 = PHSpline3P3;

    typedef std::function<ScalarCfv::real(ScalarCfv::real, ScalarCfv::real)> functionRBF;

    ScalarCfv::real RBF0(ScalarCfv::point r, ScalarCfv::real c, const functionRBF &fRBF0)
    {
        auto rr = r.length();
        return fRBF0(rr, c);
    }

    ScalarCfv::real RBF1x(ScalarCfv::point r, ScalarCfv::real c, const functionRBF &fRBF0, const functionRBF &fRBF1)
    {
        auto rr = r.length();
        if (rr < eps)
            return 0.0;
        return fRBF1(rr, c) * r.x / rr;
    }

    ScalarCfv::real RBF1y(ScalarCfv::point r, ScalarCfv::real c, const functionRBF &fRBF0, const functionRBF &fRBF1)
    {
        auto rr = r.length();
        if (rr < eps)
            return 0.0;
        return fRBF1(rr, c) * r.y / rr;
    }

    ScalarCfv::real RBF2xx(ScalarCfv::point r, ScalarCfv::real c, const functionRBF &fRBF0, const functionRBF &fRBF1, const functionRBF &fRBF2)
    {
        auto rr = r.length();
        if (rr < eps)
            return fRBF2(0, c);
        auto rr2 = rr * rr;
        auto rirj = r.x * r.x;
        auto deltaij = 1.0;
        return fRBF2(rr, c) * rirj / rr2 + fRBF1(rr, c) * (deltaij / rr - rirj / (rr2 * rr));
    }

    ScalarCfv::real RBF2xy(ScalarCfv::point r, ScalarCfv::real c, const functionRBF &fRBF0, const functionRBF &fRBF1, const functionRBF &fRBF2)
    {
        auto rr = r.length();
        if (rr < eps)
            return 0.0;
        auto rr2 = rr * rr;
        auto rirj = r.x * r.y;
        auto deltaij = 0.0;
        return fRBF2(rr, c) * rirj / rr2 + fRBF1(rr, c) * (deltaij / rr - rirj / (rr2 * rr));
    }

    ScalarCfv::real RBF2yy(ScalarCfv::point r, ScalarCfv::real c, const functionRBF &fRBF0, const functionRBF &fRBF1, const functionRBF &fRBF2)
    {
        auto rr = r.length();
        if (rr < eps)
            return fRBF2(0, c);
        auto rr2 = rr * rr;
        auto rirj = r.y * r.y;
        auto deltaij = 1.0;
        return fRBF2(rr, c) * rirj / rr2 + fRBF1(rr, c) * (deltaij / rr - rirj / (rr2 * rr));
    }

    ScalarCfv::real RBF3xxx(ScalarCfv::point r,
                            ScalarCfv::real c, const functionRBF &fRBF0,
                            const functionRBF &fRBF1,
                            const functionRBF &fRBF2,
                            const functionRBF &fRBF3)
    {
        auto rr = r.length();
        auto rirjrk = r.x * r.x * r.x;
        auto r3 = rr * rr * rr;
        auto r2 = rr * rr;
        auto r4 = r2 * r2;


        if (rr < eps)
            return fRBF3(rr, c);
        return fRBF3(rr, c) * rirjrk / r3 + fRBF2(rr, c) * 3 * (r.x / r2 - rirjrk / r4) +
               fRBF1(rr, c) * 3 * (rirjrk / (r4 * rr) - r.x / r3);
    }

    ScalarCfv::real RBF3yyy(ScalarCfv::point r,
                            ScalarCfv::real c, const functionRBF &fRBF0,
                            const functionRBF &fRBF1,
                            const functionRBF &fRBF2,
                            const functionRBF &fRBF3)
    {
        auto rr = r.length();
        auto rirjrk = r.y * r.y * r.y;
        auto r3 = rr * rr * rr;
        auto r2 = rr * rr;
        auto r4 = r2 * r2;


        if (rr < eps)
            return fRBF3(rr, c);
        return fRBF3(rr, c) * rirjrk / r3 + fRBF2(rr, c) * 3 * (r.y / r2 - rirjrk / r4) +
               fRBF1(rr, c) * 3 * (rirjrk / (r4 * rr) - r.y / r3);
    }

    ScalarCfv::real RBF3xxy(ScalarCfv::point r,
                            ScalarCfv::real c, const functionRBF &fRBF0,
                            const functionRBF &fRBF1,
                            const functionRBF &fRBF2,
                            const functionRBF &fRBF3)
    {
        auto rr = r.length();
        auto rirjrk = r.x * r.x * r.y;
        auto r3 = rr * rr * rr;
        auto r2 = rr * rr;
        auto r4 = r2 * r2;


        if (rr < eps)
            return 0.0;
        return fRBF3(rr, c) * rirjrk / r3 + fRBF2(rr, c) * (r.y / r2 - 3 * rirjrk / r4) +
               fRBF1(rr, c) * (3 * rirjrk / (r4 * rr) - r.y / r3);
    }

    ScalarCfv::real RBF3xyy(ScalarCfv::point r,
                            ScalarCfv::real c, const functionRBF &fRBF0,
                            const functionRBF &fRBF1,
                            const functionRBF &fRBF2,
                            const functionRBF &fRBF3)
    {
        auto rr = r.length();
        auto rirjrk = r.x * r.y * r.y;
        auto r3 = rr * rr * rr;
        auto r2 = rr * rr;
        auto r4 = r2 * r2;

        if (rr < eps)
            return 0.0;
        return fRBF3(rr, c) * rirjrk / r3 + fRBF2(rr, c) * (r.x / r2 - 3 * rirjrk / r4) +
               fRBF1(rr, c) * (3 * rirjrk / (r4 * rr) - r.x / r3);
    }
}

#endif