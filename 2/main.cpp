
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numbers>
#include <sstream>
// [cm]
using std::runtime_error;

const double D = 10.96;
const double d = 0.63;
const double l0 = 129.;
const double L = 115.;
// [kg]
const double m = 50. / 1000.;
// [m]
const double D_ = D / 100.;
// [m^4]
const double I_x = std::numbers::pi * D_ * D_ * D_ * D_ / 32.;

// [m]
const double r = d / 2. / 100.;

// Delta [cm]
double compute_phi(double Delta) { return Delta / (2. * L); }

// return [N]
double compute_P(size_t count) {
  // 2 * count * m * g
  // [1] * [kg] * [m / s^2]
  return count * m * 9.81;
}

// return [H * m]
double compute_M(double P) {
  // 2 * (r * m * g)
  // [m] * [kg] * [m / s^2]
  return 2 * r * P;
}

// M [H * m]
double compute_tau(double M) {
  // M * r / I
  // [H * m] * [m] / [m^4]
  return M * r / I_x;
}

double compute_gamma(double phi) {
  // phi * r / l0
  // [1] * cm / cm
  return phi * d / 2. / l0;
}

// M [H * m^2]
// return [H / m^2]
double compute_G(double M, double phi) {
  // M * l / (I * phi)
  // [H * m] * [m] / ([m^4] * [1])
  return (M * l0 / 100.) / (I_x * phi);
}

int main() {
  std::ifstream in_file("raw_data.csv");
  // std::ofstream out("out.csv");
  std::ostream &out = std::cout;
  if (!in_file || !out) {
    throw runtime_error("file not found!");
  }
  out << "N P Delta phi M tau gamma G" << std::endl;
  std::string buffer;
  while (std::getline(in_file, buffer)) {
    std::stringstream in(buffer);

    size_t N;
    double Delta;
    size_t mass_count;

    in >> N >> Delta >> mass_count;

    double P = compute_P(mass_count);
    double phi = compute_phi(Delta);
    double M = compute_M(P);
    double tau = compute_tau(M);
    double gamma = compute_gamma(phi);
    double G = compute_G(M, phi);

    out << N << " " << P << " " << Delta << " " << phi << " " << M << " " << tau
        << " " << gamma << " " << G << std::endl;
  }
  return 0;
}
