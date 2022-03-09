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
            CfvMath::getDiffBaseValueRBFB1(P(ix * 1.0 / nx, iy * 1.0 / nx),
                                           P(0.5, 0.5),
                                           P(1, 1),
                                           moment, data, c);
            for (int a = 0; a < 7; a++)
                for (int b = 0; b < 6; b++)
                    std::cout << data[a][b] << '\t';
            std::cout << ";\n";
        }
    std::cout << ']';

    return 0;
}