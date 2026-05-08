#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include "implementation.h"
#include "utils.h"

//Helper: parse CLI args
struct Config {
    int dimension = 3;
    std::vector<double> weights;
    double tolerance = 1e-9;
    std::string function = "dot"; // dot | weighted | custom
};

Config parseArgs(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "--dimension" || arg == "-d") && i + 1 < argc)
            cfg.dimension = std::stoi(argv[++i]);
        else if ((arg == "--tolerance" || arg == "-t") && i + 1 < argc)
            cfg.tolerance = std::stod(argv[++i]);
        else if ((arg == "--function" || arg == "-f") && i + 1 < argc)
            cfg.function = argv[++i];
        else if ((arg == "--weights" || arg == "-w") && i + 1 < argc) {
            ++i;
            while (i < argc && argv[i][0] != '-') {
                cfg.weights.push_back(std::stod(argv[i++]));
            }
            --i;
        }
    }
    // default weights = all 1s if not provided
    if (cfg.weights.empty())
        cfg.weights.assign(cfg.dimension, 1.0);
    return cfg;
}


// Standard dot product
auto makeDot() {
    return [](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0;
        for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
        return s;
    };
}

// Weighted inner product: <u,v> = sum(w_i * u_i * v_i)
auto makeWeighted(const std::vector<double>& w) {
    return [w](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0;
        for (size_t i = 0; i < a.size(); ++i) s += w[i] * a[i] * b[i];
        return s;
    };
}

// NON-inner-product: missing symmetry <u,v> = u1*v2 + u2*v1 (not symmetric in general)
auto makeAsymmetric() {
    return [](const std::vector<double>& a, const std::vector<double>& b) {
        if (a.size() < 2) return 0.0;
        // The 2.0 coefficient makes <u,v> different from <v,u>
        return a[0] * b[1] + 2.0 * a[1] * b[0] + (a.size() > 2 ? a[2]*b[2] : 0.0);
    };
}

// NON-inner-product: negative weight <u,v> = -u1v1 + u2v2 + ...
auto makeNegativeWeight() {
    return [](const std::vector<double>& a, const std::vector<double>& b) {
        double s = -a[0]*b[0];
        for (size_t i = 1; i < a.size(); ++i) s += a[i]*b[i];
        return s;
    };
}

// NON-inner-product: non-linear <u,v> = (u1v1)^2 + (u2v2)^2
auto makeNonLinear() {
    return [](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0;
        for (size_t i = 0; i < a.size(); ++i) s += std::pow(a[i]*b[i], 2);
        return s;
    };
}

// Extended tester
template<typename Func>
struct ExtendedTester {
    Func ip;
    double tol;
    int passCount = 0, failCount = 0;

    ExtendedTester(Func f, double t) : ip(f), tol(t) {}

    void testPositivity(const std::vector<std::vector<double>>& vecs,
                        const std::vector<double>& zero) {
        std::cout << "\n[POSITIVITY AXIOM: <v,v> >= 0, with equality iff v = 0]\n";
        bool allPass = true;

        for (size_t i = 0; i < vecs.size(); ++i) {
            double val = ip(vecs[i], vecs[i]);
            bool pass = val >= -tol;
            std::cout << "  Test " << i+1 << ": <v,v> = " << val
                      << (pass ? "  [PASS]" : "  [FAIL]") << "\n";
            if (!pass) allPass = false;
        }

        // zero vector test
        double zval = ip(zero, zero);
        bool zpass = std::abs(zval) < tol;
        std::cout << "  Zero vector: <0,0> = " << zval
                  << (zpass ? "  [PASS]" : "  [FAIL]") << "\n";
        if (!zpass) allPass = false;

        std::cout << "  >> Positivity: " << (allPass ? "PASS" : "FAIL") << "\n";
        allPass ? ++passCount : ++failCount;
    }

    void testSymmetry(const std::vector<std::pair<std::vector<double>,
                                                   std::vector<double>>>& pairs) {
        std::cout << "\n[SYMMETRY AXIOM: <u,v> = <v,u>]\n";
        bool allPass = true;

        for (size_t i = 0; i < pairs.size(); ++i) {
            double uv = ip(pairs[i].first, pairs[i].second);
            double vu = ip(pairs[i].second, pairs[i].first);
            double diff = std::abs(uv - vu);
            bool pass = diff < tol;
            std::cout << "  Test " << i+1 << ": <u,v>=" << uv
                      << "  <v,u>=" << vu
                      << "  diff=" << diff
                      << "  tol=" << tol
                      << (pass ? "  [PASS]" : "  [FAIL]") << "\n";
            if (!pass) allPass = false;
        }

        std::cout << "  >> Symmetry: " << (allPass ? "PASS" : "FAIL") << "\n";
        allPass ? ++passCount : ++failCount;
    }

    void testLinearity(const std::vector<std::tuple<std::vector<double>,
                                                     std::vector<double>,
                                                     std::vector<double>,
                                                     double, double>>& cases) {
        std::cout << "\n[LINEARITY AXIOM: <au+bv, w> = a<u,w> + b<v,w>]\n";
        bool allPass = true;

        for (size_t i = 0; i < cases.size(); ++i) {
            auto& currentCase = cases[i];
            const std::vector<double>& u = std::get<0>(currentCase);
            const std::vector<double>& v = std::get<1>(currentCase);
            const std::vector<double>& w = std::get<2>(currentCase);
            double a = std::get<3>(currentCase);
            double b = std::get<4>(currentCase);
            std::vector<double> combined(u.size());
            for (size_t j = 0; j < u.size(); ++j)
                combined[j] = a * u[j] + b * v[j];

            double left  = ip(combined, w);
            double right = a * ip(u, w) + b * ip(v, w);
            double diff  = std::abs(left - right);
            bool pass = diff < tol;

            std::cout << "  Test " << i+1 << " (a=" << a << ", b=" << b << "): "
                      << "<au+bv,w>=" << left
                      << "  a<u,w>+b<v,w>=" << right
                      << "  diff=" << diff
                      << (pass ? "  [PASS]" : "  [FAIL]") << "\n";
            if (!pass) allPass = false;
        }

        std::cout << "  >> Linearity: " << (allPass ? "PASS" : "FAIL") << "\n";
        allPass ? ++passCount : ++failCount;
    }

    void printSummary(const std::string& funcName) {
        std::cout << "\n================ SUMMARY REPORT ================\n";
        std::cout << "  Function: " << funcName << "\n";
        std::cout << "  Axioms Passed: " << passCount << " / 3\n";
        std::cout << "  Total Individual Tests: " << (passCount + failCount) * 4 << " (approx)\n"; // Estimated count

        if (failCount == 0) {
            std::cout << "  FINAL CONCLUSION: THIS IS A VALID INNER PRODUCT.\n";
        } else {
            std::cout << "  FINAL CONCLUSION: INVALID INNER PRODUCT.\n";
            std::cout << "  Reason: " << failCount << " axiom(s) failed the tolerance check.\n";
        }
        std::cout << "================================================\n\n";
    }
};

template<typename Func>
void runSuite(const std::string& label, Func ip, int dim, double tol) {
    std::cout << "\n####################################################\n";
    std::cout << "  FUNCTION: " << label << "\n";
    std::cout << "  Dimension: " << dim << "  |  Tolerance: " << tol << "\n";
    std::cout << "####################################################\n";

    std::vector<double> zero(dim, 0.0);
    std::vector<double> unit(dim, 0.0); unit[0] = 1.0;

    // build test vectors
    std::vector<double> u(dim), v(dim), w(dim), neg(dim), frac(dim);
    for (int i = 0; i < dim; ++i) {
        u[i]    =  i + 1.0;
        v[i]    = -(i + 1.0);
        w[i]    =  (i % 2 == 0) ? 0.5 : -0.5;
        neg[i]  = -(i * 2.0 + 1.0);
        frac[i] =  1.0 / (i + 1.0);
    }

    ExtendedTester<Func> tester(ip, tol);

    // ── Positivity: 4+ tests ──
    tester.testPositivity({u, v, unit, frac, neg}, zero);

    // ── Symmetry: 3+ tests ──
    tester.testSymmetry({
        {u, v},
        {unit, w},
        {frac, neg},
        {u, zero}
    });

    // ── Linearity: 3+ tests covering negative, zero, fractional scalars ──
    tester.testLinearity({
        {u, v, w,    2.0,  3.0},   // positive scalars
        {u, v, w,   -1.0,  4.0},   // negative scalar
        {u, v, w,    0.0,  1.0},   // zero scalar
        {u, v, w,    0.5, -0.5},   // fractional scalars
        {frac, neg, unit, 3.0, -2.0}
    });

    tester.printSummary(label);
}

// main
int main(int argc, char* argv[]) {
    std::cout << std::fixed << std::setprecision(6);

    Config cfg = parseArgs(argc, argv);

    std::cout << "====================================================\n";
    std::cout << "   INNER PRODUCT SPACE AXIOM CHECKER\n";
    std::cout << "====================================================\n";
    std::cout << "  Dimension : " << cfg.dimension << "\n";
    std::cout << "  Function  : " << cfg.function  << "\n";
    std::cout << "  Tolerance : " << cfg.tolerance << "\n";

    if (cfg.function == "dot" || cfg.function == "all") {
        runSuite("Standard Dot Product (VALID)", makeDot(), cfg.dimension, cfg.tolerance);
    }
    if (cfg.function == "weighted" || cfg.function == "all") {
        runSuite("Weighted Inner Product, weights=" + [&]{
            std::string s;
            for (double w : cfg.weights) s += std::to_string(w) + " ";
            return s;
        }() + " (VALID)", makeWeighted(cfg.weights), cfg.dimension, cfg.tolerance);
    }
    if (cfg.function == "custom" || cfg.function == "all") {
        runSuite("Asymmetric Function u1v2+u2v1 (INVALID - fails symmetry)",
                 makeAsymmetric(), cfg.dimension, cfg.tolerance);
        runSuite("Negative Weight -u1v1+u2v2+... (INVALID - fails positivity)",
                 makeNegativeWeight(), cfg.dimension, cfg.tolerance);
        runSuite("Non-linear (u1v1)^2+... (INVALID - fails linearity)",
                 makeNonLinear(), cfg.dimension, cfg.tolerance);
    }

    // Always run all when no specific function chosen
    if (cfg.function != "dot" && cfg.function != "weighted" && cfg.function != "custom" && cfg.function != "all") {
        std::cout << "\n[Unknown function type. Use: dot | weighted | custom | all]\n";
    }

    return 0;
}