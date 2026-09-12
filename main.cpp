
#include <array>
#include <iostream>
#include <numeric>

#include <cmath>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <tuple>
#include <vector>

template <typename T> void from_file(T &data, const std::string &file_name) {
  std::ifstream in(file_name);
  if (in.is_open()) {
    char buf[100];
    in.getline(buf, 100);

    size_t number;
    char delimer;
    double measure1;
    double measure2;
    double measure3;

    for (; in >> number >> measure1 >> measure2 >> measure3;) {
      data.push_back({measure1, measure2, measure3});
    }
  } else {
    throw std::runtime_error("file!!!");
  }
}

template <typename D, typename... Ts>
void printer(std::ostream &out, D delimer, const std::tuple<Ts...> &t) {
  [&]<size_t... Is>(std::index_sequence<Is...>) {
    ((out << (Is == 0 ? "" : delimer) << std::get<Is>(t)), ...);
  }(std::make_index_sequence<sizeof...(Ts)>{});
}
template <typename T, typename S>
auto to_file(const std::vector<std::string> &description, const T &data,
             S &&out) {
  size_t i = 0;
  const char *delimer = "\t";
  for (const auto &d : description) {
    out << d << delimer;
  }
  out << std::endl;
  for (const auto &elem : data) {
    out << i + 1;
    for (const auto &t : elem) {
      out << delimer << t;
    }
    out << std::endl;
    ++i;
  }
  return std::forward<S>(out);
}

template <typename T> double mean(T &t) {
  return std::accumulate(t.begin(), t.end(), 0.) / (double)t.size();
}

double area(double d) { return std::numbers::pi * d * d / 4.; }

int main() {
  using data_type = std::vector<std::array<double, 11>>;

  // [mm]
  std::array l_array = {30.5, 30.62, 31.};
  // [mm]
  std::array d_array = {4.4, 4.4, 3.9};
  double l0 = mean(l_array);
  double d0 = mean(d_array);
  double F0 = area(d0);
  double d0_ocul = 6.925;
  double K = d0 / d0_ocul;

  std::cout << "l0 : " << l0 << " mm" << std::endl;
  std::cout << "d0 : " << d0 << " mm" << std::endl;
  std::cout << "F0 : " << F0 << " mm^2" << std::endl;
  std::cout << "K : " << K << std::endl;

  data_type control_points;
  from_file(control_points, "./control_points_raw.csv");

  for (auto &measure : control_points) {
    // F, [N]
    measure[0] *= 2 * 10 * 9.81;
    // dl, [mm]
    measure[1] *= 2 * 0.01;
  }
  for (auto &measure : control_points) {
    // d [mm]
    measure[3] = K * measure[2];
    // F - [mm^2]
    measure[4] = area(measure[3]);

    // sigma - [Pa / cm^2]
    measure[5] = measure[0] / F0 * 100.; // <- Pa / cm^2 | Pa / mm^2
    // S - [Pa / cm^2]
    measure[6] = measure[0] / measure[4] * 100.; // <- Pa / cm^2 | Pa / mm^2

    // epsilon
    measure[7] = measure[1] / l0 * 100.;
    // psi
    measure[8] = (F0 - measure[4]) / F0 * 100.;

    // e
    measure[9] = log(1 + measure[7]);
    // psi_
    measure[10] = log(F0 / measure[4]);
  }

  to_file({"N", "P", "dl", "d'", "d", "F", "sigma", "S", "epsilon", "psi", "e",
           "psi_"},
          control_points, std::ofstream("out.csv"));

  {
    using data_type = std::array<double, 4>;
    std::vector<data_type> many_p;
    std::ifstream in("points_raw.csv");
    if (!in.is_open())
      throw std::runtime_error("file!!");
    in.ignore(10000, '\n');
    for (double m1, m2; in >> m1 >> m2;) {
      many_p.push_back({m1, m2});
    };

    for (auto &measure : many_p) {
      // P, [N]
      measure[1] *= 2 * 10 * 9.81;
      // dl, [mm]
      measure[0] *= 2 * 0.01;
      // epsilon
      measure[2] = measure[0] / l0 * 100.;
      // sigma
      measure[3] = measure[1] / F0;
    }
    std::ofstream out("out_a.csv");
    out << "dl P epsilon sigma" << std::endl;
    for (auto &measure : many_p) {
      out << measure[0] << " " << measure[1] << " " << measure[2] << " "
          << measure[3] << std::endl;
    }
  }
  double limit_of_proportionality_1 = 25;
  double limit_of_proportionality_2 = 36;
  double yield_strength = 35;
  double temporary_tear_resistance = 38;
  auto from_grid_to_force = [](double n) { return n * 2 * 10 * 9.81; };

  std::cout << from_grid_to_force(limit_of_proportionality_1) << "\n"
            << from_grid_to_force(limit_of_proportionality_2) << "\n"
            << from_grid_to_force(yield_strength) << "\n"
            << from_grid_to_force(temporary_tear_resistance);

  // residual relative narrowing
  std::cout << "psi last " << (F0 - control_points.back()[4]) / F0 * 100.
            << std::endl;

  double def_plastic = 26;
  double def_elastic = 20;
  auto from_grid_to_len = [](double n) { return n * 0.01 * 2. / 1000.; };
  std::cout << "DEF" << std::endl
            << "\telastic: " << from_grid_to_len(def_elastic) << std::endl
            << "\tplastic: " << from_grid_to_len(def_plastic) << std::endl;

  return 0;
}
