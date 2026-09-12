
#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>

#include <cmath>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <ostream>
#include <ranges>
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
  }
}

template <typename D, typename... Ts>
void printer(std::ostream &out, D delimer, const std::tuple<Ts...> &t) {
  [&]<size_t... Is>(std::index_sequence<Is...>) {
    ((out << (Is == 0 ? "" : delimer) << std::get<Is>(t)), ...);
  }(std::make_index_sequence<sizeof...(Ts)>{});
}
template <typename T>
std::ostream &to_file(const std::string &description, const T &data,
                      std::ostream &out) {
  out << description << std::endl;
  size_t i = 0;
  const char *delimer = "\t";

  for (const auto &t : data) {
    out << i + 1 << delimer;
    printer(out, delimer, t);
    out << std::endl;
    ++i;
  }
  return out;
}

template <typename T> double mean(T &t) {
  return std::accumulate(t.begin(), t.end(), 0) / (double)t.size();
}

int main() {
  using data_type = std::vector<std::array<double, 11>>;

  std::array l_array = {3.05, 3.062, 3.1};
  std::array d_array = {0.44, 0.44, 0.39};
  double l = mean(l_array);
  double d = mean(d_array);
  double d0_ocul = 6.91;
  double K = d / d0_ocul;

  data_type control_points;
  from_file(control_points, "./control_points_raw.csv");
  to_file("N P dl d", control_points, std::cout);

  return 0;
}
