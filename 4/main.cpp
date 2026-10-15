
#include "../utils.h"

#include <algorithm>
#include <ranges>
#include <tuple>
#include <vector>

// koeff
constexpr double K = 6.4e-7;
using utils::operator""_mm;
using utils::operator""_sm;

namespace lab4 {

// mm
constexpr double a = 300._mm;
constexpr double l = 260._mm;
constexpr double d = 41._mm;
constexpr double h = 1._mm;

// kg
constexpr double m = 1.;
// m / c^2
constexpr double g = 9.81;
// Pa
constexpr double E = 0.7e11;
// 1
constexpr double nu = 0.33;

double compute_sigma_max(double F) {
  return 4. * F * l / std::numbers::pi / d / d / h;
}
double compute_tau_max(double F) {
  return 2. * F * a / std::numbers::pi / d / d / h;
}

std::tuple<double, double, double> compute_sigma_angle_from_force(double F) {
  double sigma_max = compute_sigma_max(F);
  double tau_max = compute_tau_max(F);
  double sqrt_from_sigma_tau =
      std::sqrt((sigma_max / 2.) * (sigma_max / 2.) + tau_max * tau_max);

  return std::make_tuple(sigma_max / 2. + sqrt_from_sigma_tau,
                         sigma_max / 2. - sqrt_from_sigma_tau,
                         std::atan(2. * tau_max / sigma_max) / 2. * 180. /
                             std::numbers::pi);
}

auto compute_sigma_angle_from_epsilon_uzv(double e_u, double e_z, double e_v) {
  double sqrt_from_sigma_tau =
      std::sqrt(((e_u - e_z) * (e_u - e_z) + (e_v - e_z) * (e_v - e_z)) / 2.);
  double e1 = (e_u + e_v) / 2. + sqrt_from_sigma_tau;
  double e3 = (e_u + e_v) / 2. - sqrt_from_sigma_tau;
  double betta = std::atan((2. * e_z - (e_u + e_v)) / (e_u - e_v)) / 2.;
  double coeff = E / (1. - nu * nu);

  return std::make_tuple(e1, e3, coeff * (e1 + nu * e3), coeff * (e3 + nu * e1),
                         (std::numbers::pi / 4. - betta) * 180. /
                             std::numbers::pi);
}

using data_type = std::tuple<size_t, double, double, double>;

auto load_and_prepare() {
  auto data1 = utils::read_file_csv<data_type>("./raw_data_1.csv");
  {
    auto first = data1.front();
    auto f = [&first](auto &element) {
      const auto &[N, td1, td2, td3] = element;
      auto [_, td1_0, td2_0, td3_0] = first;
      return std::make_tuple(N - 1., td1 - td1_0, td2 - td2_0, td3 - td3_0);
    };
    std::ranges::transform(data1, data1.begin(), f);
  }
  auto data2 = utils::read_file_csv<data_type>("./raw_data_2.csv");
  {
    auto first = data2.front();
    auto f = [first](auto &element) {
      const auto &[N, td1, td2, td3] = element;
      auto [_, td1_0, td2_0, td3_0] = first;
      return std::make_tuple(6 - N, td1 - td1_0, td2 - td2_0, td3 - td3_0);
    };
    std::ranges::transform(data2, data2.begin(), f);
  }
  std::vector<data_type> result;
  result.insert(result.end(), data1.begin() + 1, data1.end());
  result.insert(result.end(), data2.begin() + 1, data2.end() - 1);

  return result;
}
auto compute(std::ranges::view auto &&data) {
  auto f2 = [](const auto &element) {
    const auto &[n, d_td1, d_td2, d_td3] = element;

    return std::make_tuple(n * m * g, K * d_td1, K * d_td2, K * d_td3);
  };
  auto f3 = [](const auto &element) {
    const auto &[F, e_u, e_z, e_v] = element;

    return std::tuple_cat(element,
                          compute_sigma_angle_from_epsilon_uzv(e_u, e_z, e_v),
                          compute_sigma_angle_from_force(F));
  };
  return data | std::views::transform(f2) | std::views::transform(f3);
}

void main() {
  auto description = "F\te_u\te_z\te_v\te_1\te_"
                     "3\tsigma1\tsigma3\talpha\tsigma1\tsigma3\talpha";
  auto data_owner = load_and_prepare();

  auto computed = compute(data_owner | std::views::all);
  std::vector result(computed.begin(), computed.end());
  utils::print_result(description, result, std::ofstream("out1.csv"));
}
} // namespace lab4

namespace lab5 {
// N
constexpr double P = 50.;

// cm
constexpr auto d = 2.86_sm;
constexpr auto h = 0.246_sm;
// m^2
constexpr double S0 = d * h;

auto compute_sigma_xx(double P) { return P / S0; }
auto compute_E(double P, double e_xx) { return P / S0 / e_xx; }
auto compute_nu(double e_xx, double e_yy) { return -1 * e_yy / e_xx; }

using data_type = std::tuple<size_t, double, double, double, double>;

auto load_and_prepare(const std::string &in) {
  auto data = utils::read_file_csv<data_type>(in);
  auto first = data.front();
  auto f1 = [&first](auto &element) {
    const auto &[N, td1, td2, td3, td4] = element;
    auto [_, td1_0, td2_0, td3_0, td4_0] = first;
    return std::make_tuple(N, td1 - td1_0, td2 - td2_0, td3 - td3_0,
                           td4 - td4_0);
  };
  std::ranges::transform(data, data.begin(), f1);
  data.erase(data.begin());
  return data;
}
auto compute(std::ranges::view auto &&data) {
  auto f2 = [](const auto &element) {
    const auto &[n, d_td1, d_td2, d_td3, d_td4] = element;
    return std::make_tuple(n * P, (d_td1 + d_td4) / 2., (d_td2 + d_td3) / 2.);
  };
  auto f3 = [](const auto &element) {
    const auto &[P, td_xx, td_yy] = element;
    return std::tuple_cat(element, std::make_tuple(td_xx * K, td_yy * K));
  };
  auto f4 = [](const auto &element) {
    const auto &[P, td_xx, td_yy, e_xx, e_yy] = element;
    return std::tuple_cat(element, std::make_tuple(compute_sigma_xx(P),
                                                   compute_nu(e_xx, e_yy),
                                                   compute_E(P, e_xx)));
  };

  return data | std::views::transform(f2) | std::views::transform(f3) |
         std::views::transform(f4);
}

void main() {
  auto description = "F\ttd_xx\ttd_yy\tepsilon_xx\tepsilon_yy\tsigma_xx\tnu\tE";
  std::array all_data{load_and_prepare("./raw_data_3.csv"),
                      load_and_prepare("./raw_data_4.csv"),
                      load_and_prepare("./raw_data_5.csv")};
  auto data = compute(
      all_data | std::views::join | std::views::filter([](const auto &element) {
        return std::get<1>(element) >= 0. && std::get<4>(element) >= 0. &&
               std::get<2>(element) <= 0. && std::get<3>(element) <= 0.;
      }));
  std::vector result(data.begin(), data.end());
  utils::print_result(description, data, std::ofstream("out2.csv"));

  auto [nu, s_nu] = utils::compute_data_stat(result | std::views::elements<6>);
  auto [E, s_E] = utils::compute_data_stat(result | std::views::elements<7>);

  std::cout << std::format("nu (0.98): {:e} +- {:e}", nu,
                           s_nu * utils::student_coef(result.size(), 0.98))
            << "\n"
            << std::format("E (0.98): {:e} +- {:e}", E,
                           s_E * utils::student_coef(result.size(), 0.98))
            << "\n";
}

} // namespace lab5

int main() {
  lab4::main();
  lab5::main();
  return 0;
}
