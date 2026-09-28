
#include <array>
#include <boost/math/distributions/students_t.hpp>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numbers>
#include <ranges>
#include <sstream>
#include <utility>
#include <vector>

using std::tuple;
using std::ranges::views::transform;

static double student_coef(size_t n, double alpha) {
  if (n <= 1) {
    return 0.0;
  }
  boost::math::students_t dist(n - 1);
  /*
   *                     alpha
   <----|--------------------------------------|---->
      -t_stat              0                 +t_stat
  */
  auto t_stat = boost::math::quantile(dist, 0.5 + alpha / 2.0);
  return t_stat;
}

template <typename T, size_t N = 0> void read(std::istream &in, T &value) {
  in >> std::get<N>(value);
  if constexpr (N + 1 < std::tuple_size_v<T>) {
    read<T, N + 1>(in, value);
  }
}
template <typename T, typename D, size_t N = 0>
void print(std::ostream &out, const T &value, const D &delimer) {
  out << std::get<N>(value);
  if constexpr (N + 1 < std::tuple_size_v<T>) {
    out << delimer;
    print<T, D, N + 1>(out, value, delimer);
  }
}

template <typename T> auto read_file(const std::string &file_name) {
  std::ifstream in_file(file_name);

  if (!in_file) {
    throw std::runtime_error("file not found!");
  }
  in_file.ignore(1000, '\n');

  std::string buffer;
  std::vector<T> in_data;

  while (std::getline(in_file, buffer)) {
    std::stringstream in(buffer);

    T value;

    read(in, value);
    in_data.push_back(value);
  }

  return in_data;
}

// d [mkm]
// P [g]
double compute_microhardness(double d, double P) {
  return 2. * P * std::sin(std::numbers::pi * 68. / 180.) / (d * d) * 1000.;
}

const double P = 50;

std::pair<double, double> get_s(auto &&data) {
  double mean = 0., s = 0.;

  for (auto &&d : data) {
    mean += d;
  }
  mean /= data.size();
  for (auto &&d : data) {
    s += (d - mean) * (d - mean);
  }
  s /= data.size();
  s /= data.size() - 1;
  s = std::sqrt(s);

  return {mean, s};
}

void micro_hardness_1(const std::string &in_file_name,
                      const std::string &out_file_name) {
  using data_type = std::tuple<double, double>;
  auto data_range =
      read_file<data_type>(in_file_name) | transform([](const data_type &a) {
        auto h = compute_microhardness(std::get<0>(a), P);
        return std::tuple_cat(a, std::make_tuple(h));
      });
  std::vector data(data_range.begin(), data_range.end());

  std::ofstream out(out_file_name);
  out << "N\td\tH_lab\tH_computed\n";
  for (size_t i = 0; i < data.size(); ++i) {
    print(out, std::tuple_cat(std::make_tuple(i), data[i]), "\t");
    out << "\n";
  }
  std::cout << "FILE : " << in_file_name << std::endl;
  {
    auto [mean, s] = get_s(data | std::views::elements<1>);
    std::cout << mean << " +- " << s * student_coef(data.size(), 0.98)
              << std::endl;
  }
  {
    auto [mean, s] = get_s(data | std::views::elements<2>);
    std::cout << mean << " +- " << s * student_coef(data.size(), 0.98)
              << std::endl;
  }
}
void micro_hardness_time(const std::string &in_file_name,
                         const std::string &out_file_name) {
  using data_type = std::tuple<double, double, double>;
  auto data_range =
      read_file<data_type>(in_file_name) | transform([](const data_type &a) {
        auto [d, h_, t] = a;
        auto h = compute_microhardness(d, P);
        return std::tuple_cat(a, std::make_tuple(h));
      });
  std::vector data(data_range.begin(), data_range.end());

  std::ofstream out(out_file_name);
  for (size_t i = 0; i < data.size(); ++i) {
    print(out, std::tuple_cat(std::make_tuple(i), data[i]), "\t");
    out << "\n";
  }
  std::cout << "FILE : " << in_file_name << std::endl;
  {
    auto [mean, s] = get_s(data | std::views::elements<1>);
    std::cout << mean << " +- " << s * student_coef(data.size(), 0.98)
              << std::endl;
  }
  {
    auto [mean, s] = get_s(data | std::views::elements<2>);
    std::cout << mean << " +- " << s * student_coef(data.size(), 0.98)
              << std::endl;
  }
}

int main() {
  micro_hardness_1("./raw_data_1.csv", "./out1.csv");
  micro_hardness_1("./raw_data_3.csv", "./out2.csv");
  return 0;
}
