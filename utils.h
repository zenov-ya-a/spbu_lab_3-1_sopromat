
#ifndef __utils_h__
#define __utils_h__

#include <algorithm>
#include <boost/math/distributions/students_t.hpp>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <numbers>
#include <numeric>
#include <ranges>
#include <sstream>
#include <utility>
#include <vector>

namespace utils {

template <typename T>
concept TupleLike = requires(T t) {
  // Проверяем, что есть tuple_size
  {
    std::tuple_size<std::remove_cvref_t<T>>::value
  } -> std::convertible_to<std::size_t>;
  // Проверяем, что работает std::get<0> (если размер > 0)
  requires(std::tuple_size_v<std::remove_cvref_t<T>> == 0) ||
              requires { std::get<0>(t); };
};

template <typename T, typename D, size_t N = 0>
  requires TupleLike<T>
void print(std::ostream &out, const T &value, const D &delimer) {
  out << std::get<N>(value);
  if constexpr (N + 1 < std::tuple_size_v<T>) {
    out << delimer;
    print<T, D, N + 1>(out, value, delimer);
  }
}

template <typename T, size_t N = 0>
  requires TupleLike<T>
void read(std::istream &in, T &value) {
  in >> std::get<N>(value);
  if constexpr (N + 1 < std::tuple_size_v<T>) {
    read<T, N + 1>(in, value);
  }
}

template <typename T> auto read_file_csv(const std::string &file_name) {
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

void print_result(const std::string &description,
                  std::ranges::input_range auto &&data, std::ostream &out) {
  out << description;
  std::ranges::for_each(data, [&out](const auto &element) {
    out << std::endl;
    utils::print(out, element, '\t');
  });
}
void print_result(const std::string &description,
                  std::ranges::input_range auto &&data, std::ostream &&out) {
  print_result(description, data, out);
}

template <typename T>
concept DataLike = std::ranges::input_range<T> && std::ranges::sized_range<T>;

double compute_mean(DataLike auto &&data) {
  double mean = std::accumulate(data.begin(), data.end(), 0.) / data.size();
  return mean;
}

// sqrt((n - mean)^2 / n / (n - 1))
double compute_standart_deviation(DataLike auto &&data, double mean) {
  auto f = [mean](auto d) { return (mean - d) * (mean - d); };
  double s = std::ranges::fold_left(data | std::views::transform(f), 0.,
                                    std::plus<>{});
  s /= data.size();
  s /= data.size() - 1;
  s = std::sqrt(s);
  return s;
}

std::pair<double, double> compute_data_stat(DataLike auto &&data) {
  double mean = compute_mean(data);
  double s = compute_standart_deviation(data, mean);
  return std::make_pair(mean, s);
}

inline double student_coef(size_t n, double alpha) {
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

constexpr long double mm = 0.001;
constexpr long double sm = 0.01;

inline constexpr long double operator""_mm(unsigned long long x) {
  return x * mm;
}
inline constexpr long double operator""_mm(long double x) { return x * mm; }
inline constexpr long double operator""_sm(unsigned long long x) {
  return x * sm;
}
inline constexpr long double operator""_sm(long double x) { return x * sm; }

} // namespace utils

#endif // !__utils_h__
