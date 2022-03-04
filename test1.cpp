#include "RBF.h"
#include <iostream>

int main()
{
    using P=ScalarCfv::point;
    std::cout << RBF::RBF3xxx(P(0,0), 2, RBF::F0,RBF::F1,RBF::F2,RBF::F3) << std::endl;
    std::cout << RBF::RBF3xxy(P(1, 2), 2, RBF::F0, RBF::F1, RBF::F2, RBF::F3) << std::endl;
    std::cout << RBF::RBF3xyy(P(3, 2), 1, RBF::F0, RBF::F1, RBF::F2, RBF::F3) << std::endl;
    std::cout << RBF::RBF3yyy(P(1, 2), 1, RBF::F0, RBF::F1, RBF::F2, RBF::F3) << std::endl;
    std::cout << RBF::RBF2xx(P(1, 2), 3, RBF::F0, RBF::F1, RBF::F2) << std::endl;
    std::cout << RBF::RBF2xy(P(2, 2), 1, RBF::F0, RBF::F1, RBF::F2) << std::endl;
    std::cout << RBF::RBF2yy(P(3, 2), 1, RBF::F0, RBF::F1, RBF::F2) << std::endl;
    std::cout << RBF::RBF1x(P(1, 4), 4, RBF::F0, RBF::F1) << std::endl;
    std::cout << RBF::RBF1y(P(3, 3), 1, RBF::F0, RBF::F1) << std::endl;
    return 0;
}