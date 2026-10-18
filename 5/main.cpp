
#include "../utils.h"

#include <algorithm>
#include <ranges>
#include <tuple>
#include <vector>

namespace lab6 {

using utils::mm;
using utils::sm;
using utils::operator""_sm;
using utils::operator""_mm;

// H
constexpr double F = 0.1 * 9.8;

constexpr double K = 0.01_mm;

constexpr double d = 3.6_sm;
constexpr double h = 0.55_sm;

constexpr double l1 = 7_sm, l2 = (1 + 17.5 + 9.) * sm,
                 l3 = (8.5 + 28.5 + 1) * sm, l4 = 7_sm;
constexpr double L = l2 + l3;

// ---------|----------------|----------------|-----------------
//   l1          l2               l3              l4
//                           L
//                    |       p_dist   |

// m^4
constexpr double I_xx = d * h * h * h / 12.;
// Pa
constexpr double E = 2.5e11;

using data_type = std::tuple<double, double, double, double>;

struct aboba {
  static double compute_R1(double p_dist, double P) {
    return (L - p_dist) / L * P;
  }
  static double compute_R2(double p_dist, double P) { return p_dist / L * P; }
  static double compute_c1(double p_dist, double P) {
    return -1 *
           (std::pow(l2 + l3, 3) * compute_R1(p_dist, P) / 6. -
            std::pow(l2 + l3 - p_dist, 3) * P / 6.) /
           (l2 + l3);
  }
  static double compute_c2(double c1) { return -l1 * c1; }
  aboba(double p_dist, double P)
      : p_dist(p_dist), P(P), R1(compute_R1(p_dist, P)),
        R2(compute_R2(p_dist, P)), c1(compute_c1(p_dist, P)),
        c2(compute_c2(c1)) {}

public:
  double y(double z) { return y_private(z) / E / I_xx; }
  double tetta(double z) { return tetta_private(z) / E / I_xx; }

private:
  double y_private(double z) {
    if (z < l1) {
      return c2 + c1 * z;
    }
    if (z < l1 + p_dist) {
      return c2 + c1 * z + std::pow(z - l1, 3) * R1 / 6.;
    }
    if (z < l1 + l2 + l3) {
      return c2 + c1 * z + std::pow(z - l1, 3) * R1 / 6. -
             std::pow(z - p_dist - l1, 3) * P / 6.;
    }
    return c2 + c1 * z + std::pow(z - l1, 3) * R1 / 6. -
           std::pow(z - p_dist - l1, 3) * P / 6. +
           std::pow(z - l1 - l2 - l3, 3) * R2 / 6.;
  }
  double tetta_private(double z) {

    if (z < l1) {
      return c1;
    }
    if (z < l1 + p_dist) {
      return c1 + std::pow(z - l1, 2) * R1 / 2.;
    }
    if (z < l1 + l2 + l3) {
      return c1 + std::pow(z - l1, 2) * R1 / 2. -
             std::pow(z - p_dist - l1, 2) * P / 2.;
    }
    return c1 + std::pow(z - l1, 2) * R1 / 2. -
           std::pow(z - p_dist - l1, 2) * P / 2. +
           std::pow(z - l1 - l2 - l3, 2) * R2 / 2.;
  }
  double p_dist, P, R1, R2, c1, c2;
};

auto compute(double p_dist, const std::string &file_name) {
  return utils::read_file_csv<data_type>(file_name) |
         std::views::transform([](const auto &e) {
           auto [N, t1, t2, t3] = e;
           return std::make_tuple(N * F, t1 * K, t2 * K, t3 * K);
         }) |
         std::views::transform([p_dist](const auto &e) {
           double k = 1000.;
           auto [P, y1, y2, y3] = e;
           auto abo = aboba(p_dist, P);
           auto tetta_exp = std::make_tuple(-y1 / l1, y2 / l2, y3 / l3);
           auto y_comp = std::make_tuple(abo.y(0) * k, abo.y(l1 + l2) * k,
                                         abo.y(l1 + l2 + l3 + l4) * k);

           auto tetta_comp = std::make_tuple(abo.tetta(0), abo.tetta(l1 + l2),
                                             abo.tetta(l1 + l2 + l3 + l4));
           return std::tuple_cat(std::make_tuple(P, y1 * k, y2 * k, y3 * k),
                                 tetta_exp, y_comp, tetta_comp);
         });
}

void main() {
  double p_dist_1 = (1 + 17.5) * sm;
  double p_dist_2 = (1 + 17.5 + 9 + 8.5) * sm;

  auto descripsion =
      "P\ty1_exp\ty2_exp\ty3_exp\ttetta1_exp\ttetta2_exp\ttetta3_exp\ty1_"
      "comp\ty2_comp\ty3_comp\ttetta1_comp\ttetta2_comp\ttetta3_comp";
  utils::print_result(descripsion, compute(p_dist_1, "./raw_data_1.csv"),
                      std::ofstream("out1.csv"));

  utils::print_result(descripsion, compute(p_dist_2, "./raw_data_2.csv"),
                      std::ofstream("out2.csv"));
  auto f = [](double z, double p) {
    if (z < l1) {
      return 0.;
    }
    if (z < l1 + p) {
      return (z - l1) * (L - p) / L;
    }
    if (z < l1 + l2 + l3) {
      return (z - l1) * (L - p) / L - (z - p - l1);
    }
    return (z - l1) * (L - p) / L - (z - p - l1) + (z - l1 - l2 - l3) * p / L;
  };
  auto gen = [&f](double p) {
    std::vector<double> points({0., l1, p, L - p, l4});
    std::partial_sum(points.cbegin(), points.cend(), points.begin());

    return std::ranges::owning_view(std::move(points)) |
           std::views::transform([&f, p](const double &e) {
             return std::make_tuple(e, f(e, p));
           });
  };
  utils::print_result("z\tM", gen(p_dist_1), std::ofstream("./out_p1.csv"));
  utils::print_result("z\tM", gen(p_dist_2), std::ofstream("./out_p2.csv"));
}
} // namespace lab6

namespace lab7 {

void main() {}

} // namespace lab7

int main() {
  lab6::main();
  lab7::main();
  return 0;
}
