#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <iomanip>

template <typename T, typename Func>
class AxiomTester
{
private:
    std::function<double(const T &, const T &)> innerProduct;
    double tolerance = 1e-9;

public:
    AxiomTester(Func ip) : innerProduct(ip) {}

    // positivity axiom
    void testPositivity(const T &v, const T &zeroV)
    {
        std::cout << "[Testing Positivity]" << std::endl;
        double selfIP = innerProduct(v, v);
        double zeroIP = innerProduct(zeroV, zeroV);

        bool cond1 = selfIP >= -tolerance;         // <v,v> >= 0
        bool cond2 = std::abs(zeroIP) < tolerance; // <0,0> = 0

        std::cout << " - <v, v> = " << selfIP << (cond1 ? " (>= 0 OK)" : " (Fail)") << std::endl;
        std::cout << " - <0, 0> = " << zeroIP << (cond2 ? " (Zero Vector OK)" : " (Fail)") << std::endl;
    }

    // b) symmetry axiom
    void testSymmetry(const T &u, const T &v)
    {
        std::cout << "[Testing Symmetry]" << std::endl;
        double uv = innerProduct(u, v);
        double vu = innerProduct(v, u);

        bool match = std::abs(uv - vu) < tolerance;
        std::cout << " - <u, v>: " << uv << " | <v, u>: " << vu
                  << (match ? " (Symmetry OK)" : " (Fail)") << std::endl;
    }

    // linearity axiom
    void testLinearity(const T &u, const T &v, const T &w, double a, double b)
    {
        std::cout << "[Testing Linearity]" << std::endl;

        // left: <au + bv, w>
        T combined = (u * a) + (v * b); // Needs operator+ and operator*
        double leftSide = innerProduct(combined, w);

        // right: a<u, w> + b<v, w>
        double rightSide = a * innerProduct(u, w) + b * innerProduct(v, w);

        bool match = std::abs(leftSide - rightSide) < tolerance;
        std::cout << " - Scalar a=" << a << ", b=" << b << std::endl;
        std::cout << " - <au+bv, w>: " << leftSide << " | a<u,w> + b<v,w>: " << rightSide
                  << (match ? " (Linearity OK)" : " (Fail)") << std::endl;
    }
};

template <typename T, typename Func>
T projection(const T &u, const T &v, Func innerProduct)
{
    double numerator = innerProduct(u, v);
    double denominator = innerProduct(v, v);

    if (std::abs(denominator) < 1e-9)
        return v * 0.0; // ensuring no division by zero so program doesn't throw an exception !!

    return v * (numerator / denominator);
}

template <typename T, typename Func>
std::vector<T> gramSchmidt(const std::vector<T> &basis, Func innerProduct)
{
    std::vector<T> orthogonalBasis;
    for (const auto &v : basis)
    {
        T w = v;
        for (const auto &prevW : orthogonalBasis)
        {
            // passing the innerProduct through to projection
            w = w - projection<T, Func>(v, prevW, innerProduct);
        }
        orthogonalBasis.push_back(w);
    }

    std::vector<T> orthonormalBasis;
    for (const auto &w : orthogonalBasis)
    {
        double norm = std::sqrt(innerProduct(w, w));
        if (norm > 1e-9)
            orthonormalBasis.push_back(w * (1.0 / norm));
    }
    return orthonormalBasis;
}
// verifies where the given basis is actually orthogonal or not
template <typename T, typename Func>
void verifyOrthonormality(const std::vector<T> &basis, Func innerProduct)
{
    std::cout << "\n--- Verification of Orthonormality ---" << std::endl;
    for (size_t i = 0; i < basis.size(); ++i)
    {
        for (size_t j = i; j < basis.size(); ++j)
        {
            double res = innerProduct(basis[i], basis[j]);
            std::cout << "<u" << i + 1 << ", u" << j + 1 << "> = "
                      << (std::abs(res) < 1e-9 ? 0.0 : res) << std::endl;
        }
    }
}
