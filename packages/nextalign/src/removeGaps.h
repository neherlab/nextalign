#pragma once

#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <string>

inline void removeGapsInPlace(std::string& seq) {
  boost::remove_erase_if(seq, boost::is_any_of("-"));
}

inline std::string removeGaps(const std::string_view& seq) {
  auto copy = std::string(seq);
  removeGapsInPlace(copy);
  return copy;
}
