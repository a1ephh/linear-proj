#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <iomanip>

// 1st implemention; vector in R^N: vector in R3

#include <vector>

// overload for vector - vector
std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b)
{
    std::vector<double> result;
    result.reserve(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        result.push_back(a[i] - b[i]);
    }
    return result;
}

// overload for vector * scalar
std::vector<double> operator*(const std::vector<double> &v, double scalar)
{
    std::vector<double> result;
    result.reserve(v.size());
    for (size_t i = 0; i < v.size(); ++i)
    {
        result.push_back(v[i] * scalar);
    }
    return result;
}

struct VectorR3
{
    double x, y, z;

    VectorR3 operator-(const VectorR3 &other) const { return {x - other.x, y - other.y, z - other.z}; }
    VectorR3 operator*(double s) const { return {x * s, y * s, z * s}; }

    static double dotProduct(const VectorR3 &a, const VectorR3 &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
};

// 2nd implementation: polynomials
struct Polynomial
{
    std::vector<double> coeffs; // where {a, b} corresponds to a + bx

    double eval(double x) const
    {
        double res = 0;
        for (size_t i = 0; i < coeffs.size(); ++i)
            res += coeffs[i] * std::pow(x, i);
        return res;
    }

    Polynomial operator-(const Polynomial &other) const
    {
        size_t n = std::max(coeffs.size(), other.coeffs.size());
        std::vector<double> res(n, 0.0);
        for (size_t i = 0; i < n; ++i)
        {
            if (i < coeffs.size())
                res[i] += coeffs[i];
            if (i < other.coeffs.size())
                res[i] -= other.coeffs[i];
        }
        return {res};
    }
};

// integral inner product for polynomials on [0, 1]
double polyInnerProduct(const Polynomial &p, const Polynomial &q)
{
    double sum = 0;
    int steps = 1000;
    double dx = 1.0 / steps;
    for (int i = 0; i < steps; ++i)
    {
        double x = (i + 0.5) * dx;
        sum += p.eval(x) * q.eval(x);
    }
    return sum * dx;
}

inline std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b)
{
    std::vector<double> res(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        res[i] = a[i] + b[i];
    return res;
}

inline Polynomial operator+(const Polynomial &a, const Polynomial &b)
{
    size_t n = std::max(a.coeffs.size(), b.coeffs.size());
    std::vector<double> res(n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        if (i < a.coeffs.size())
            res[i] += a.coeffs[i];
        if (i < b.coeffs.size())
            res[i] += b.coeffs[i];
    }
    return {res};
}

// 2. Also ensure your Polynomial operator* is defined similarly if it's outside:
inline Polynomial operator*(const Polynomial &p, double s)
{
    std::vector<double> res = p.coeffs;
    for (double &c : res)
        c *= s;
    return {res};
}