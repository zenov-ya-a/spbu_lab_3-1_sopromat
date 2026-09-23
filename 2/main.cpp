
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numbers>
#include <ranges>
#include <sstream>
#include <vector>

using std::runtime_error;

// [cm]
const double d = 0.63;
const double D = 10.96;
const double L = 115.;
const double l0 = 129.;
// [g]
const double m = 50.;

// [m^4]
const double I_x = std::numbers::pi * pow(d / 100., 4) / 32.;

// [cm]
const double d_d = 0.005;
const double d_D = 0.005;
const double d_L = 0.1;
const double d_l0 = 0.1;
const double d_Delta = 0.1;
// [g]
const double d_m = 1;

// Delta [cm]
double compute_phi(double Delta) { return Delta / (2. * L); }

// return [N]
double compute_P(size_t count) {
  // count * m * g
  // [1] * [kg] * [m / s^2]
  return count * m / 1000. * 9.81;
}

// return [H * m]
double compute_M(double P) {
  // D * P
  // [m] * [H]
  return D / 100. * P;
}

// M [H * m]
// return [Pa]
double compute_tau(double M) {
  // M * r / I
  // [H * m] * [m] / [m^4]
  return M * d / 2. / 100. / I_x;
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
  std::ofstream out_file("out.csv");
  std::ostream &out = out_file; // std::cout;
  if (!in_file || !out) {
    throw runtime_error("file not found!");
  }
  in_file.ignore(1000, '\n');
  std::string buffer;
  std::vector<std::tuple<size_t, double, size_t>> in_data;

  while (std::getline(in_file, buffer)) {
    std::stringstream in(buffer);

    size_t N;
    double Delta;
    size_t mass_count;

    in >> N >> Delta >> mass_count;
    in_data.push_back({N, Delta, mass_count});
  }

  out << "N" << "\t" << "P" << "\t" << "Delta" << "\t" << "phi" << "\t" << "M"
      << "\t" << "tau" << "\t" << "gamma"
      << "\t" << "G" << std::endl;

  double Delta0 = std::get<1>(in_data.front());

  for (auto [N, Delta, mass_count] : in_data) {
    Delta -= Delta0;
    double P = compute_P(mass_count);
    double phi = compute_phi(Delta);
    double M = compute_M(P);
    double tau = compute_tau(M);
    double gamma = compute_gamma(phi);
    double G = compute_G(M, phi);

    out << N << "\t" << P << "\t" << Delta << "\t" << phi << "\t" << M << "\t"
        << tau << "\t" << gamma << "\t" << G << std::endl;
  }

  double mean_G = 0.;
  double mean_Delta = 0.;
  double mean_mass_count = 0.;
  std::ranges::subrange in_data_computed{in_data.begin() + 1,
                                         in_data.end() - 1};
  for (auto [N, Delta, mass_count] : in_data_computed) {
    Delta -= Delta0;
    double P = compute_P(mass_count);
    double phi = compute_phi(Delta);
    double M = compute_M(P);
    double G = compute_G(M, phi);

    mean_G += G;
    mean_Delta += Delta;
    mean_mass_count += mass_count;
  }

  mean_G /= (double)in_data_computed.size();
  mean_Delta /= (double)in_data_computed.size();
  mean_mass_count /= (double)in_data_computed.size();
  std::printf("mean\n\tG : %f\n\tDelta : %f\n\tmass_count : %f\n", mean_G,
              mean_Delta, mean_mass_count);

  double epsilon_d = d_d / d;
  double epsilon_Delta = d_Delta / mean_Delta;
  double epsilon_m = d_m / m;
  double epsilon_D = d_D / D;
  double epsilon_l0 = d_l0 / l0;
  double epsilon_L = d_L / L;
  double epsilon_G =
      sqrt(pow(4 * epsilon_d, 2) + pow(epsilon_Delta, 2) + pow(epsilon_m, 2) +
           pow(epsilon_D, 2) + pow(epsilon_l0, 2) + pow(epsilon_L, 2));

  std::printf("%f %f %f %f %f %f\n\tepsilon_G :  %f\n", epsilon_d,
              epsilon_Delta, epsilon_m, epsilon_D, epsilon_l0, epsilon_L,
              epsilon_G);

  double Delta_G = mean_G * epsilon_G;

  std::cout << mean_G << " +- " << Delta_G << std::endl;

  std::cout << I_x << "\n";

  return 0;
}
