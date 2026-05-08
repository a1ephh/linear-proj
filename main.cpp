#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "implementation.h"
#include "utils.h"

int main()
{
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "====================================================\n";
    std::cout << "   EXHIBITION: INNER PRODUCT SPACE VALIDATION       \n";
    std::cout << "====================================================\n"
              << std::endl;

    auto dot = [](const std::vector<double> &a, const std::vector<double> &b)
    {
        double s = 0;
        for (size_t i = 0; i < a.size(); ++i)
            s += a[i] * b[i];
        return s;
    };

    std::vector<double> u = {1.0, 2.0, 0.0};
    std::vector<double> v = {-1.0, 1.0, 2.0};
    std::vector<double> w = {0.0, 0.0, 1.0};
    std::vector<double> zero = {0.0, 0.0, 0.0};

    std::cout << "--- Phase 1: Validating R^3 Axioms ---" << std::endl;
    AxiomTester<std::vector<double>, decltype(dot)> testerRn(dot);
    testerRn.testPositivity(u, zero);
    testerRn.testSymmetry(u, v);
    testerRn.testLinearity(u, v, w, 0.5, -2.0);

    std::cout << "\n--- Phase 2: Gram-Schmidt Implementation (R^3) ---" << std::endl;
    std::vector<std::vector<double>> basisRn = {{1, 1, 0}, {1, 2, 0}, {0, 1, 2}};
    auto onbRn = gramSchmidt(basisRn, dot);
    for (size_t i = 0; i < onbRn.size(); ++i)
    {
        std::cout << "u" << i + 1 << ": (" << onbRn[i][0] << ", " << onbRn[i][1] << ", " << onbRn[i][2] << ")" << std::endl;
    }

    std::cout << "\n\n====================================================\n";
    std::cout << "   EXHIBITION: POLYNOMIAL SPACE VALIDATION          \n";
    std::cout << "====================================================\n"
              << std::endl;

    Polynomial p1 = {{1, 1}}; // 1 + x
    Polynomial p2 = {{0, 2}}; // 2x
    Polynomial p3 = {{1, 0}}; // 1
    Polynomial pZero = {{0, 0}};

    std::cout << "--- Phase 1: Validating P1 [0,1] Axioms ---" << std::endl;
    AxiomTester<Polynomial, decltype(polyInnerProduct)> testerPoly(polyInnerProduct);
    testerPoly.testPositivity(p1, pZero);
    testerPoly.testSymmetry(p1, p2);
    testerPoly.testLinearity(p1, p2, p3, 1.5, 3.0);

    std::cout << "\n--- Phase 2: Gram-Schmidt Implementation (P1) ---" << std::endl;
    std::vector<Polynomial> basisPoly = {{{1.0, 0.0}}, {{0.0, 1.0}}}; // {1, x}
    auto onbPoly = gramSchmidt(basisPoly, polyInnerProduct);
    for (size_t i = 0; i < onbPoly.size(); ++i)
    {
        std::cout << "u" << i + 1 << ": " << onbPoly[i].coeffs[0] << " + (" << onbPoly[i].coeffs[1] << ")x" << std::endl;
    }

    return 0;
}