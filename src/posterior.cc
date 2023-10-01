#include <iostream>
#include <cmath>
#include <vector>

double f(double r) {
    return std::sqrt((1.0 / (48.0 * M_PI)) * (3.0 * std::cosh(2.0 * r) + 1.0));
}

double f_prime(double r) {
    return std::sinh(2.0 * r) / (4.0 * std::sqrt(M_PI * std::cosh(2.0 * r) + M_PI / 3.0));
}

double pc_prior(double kappa, const std::vector<double>& v, double lambda, double lambda1) {
    double v_magnitude = std::sqrt(v[0] * v[0] + v[1] * v[1]);
    double f_val = f(v_magnitude);
    double f_prime_val = f_prime(v_magnitude);
    double exp_term = std::exp(-lambda1 * (f_val - f(0.0)) - lambda * f_val * kappa);

    return (lambda * lambda1 * f_prime_val * f_val) / (2.0 * M_PI * v_magnitude) * exp_term;
}

double log_pc_prior(const double kappa, const std::vector<double>& v, const double lambda, const double lambda1) {
    double v_magnitude = std::sqrt(v[0] * v[0] + v[1] * v[1]);
    double f_val = f(v_magnitude);
    double f_prime_val = f_prime(v_magnitude);

    return std::log(lambda) + std::log(lambda1) + std::log(f_prime_val) + std::log(f_val) - std::log(2.0 * M_PI * v_magnitude)
           - lambda1 * (f_val - f(0.0)) - lambda * f_val * kappa;
}