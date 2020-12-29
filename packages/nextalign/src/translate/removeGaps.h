#pragma once

#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <string>

template<typename Letter>
inline void removeGapsInPlace(Sequence<Letter>& seq) {
  seq.erase(std::remove(boost::begin(seq), boost::end(seq), Letter::GAP), boost::end(seq));
}

template<typename Letter>
inline Sequence<Letter> removeGaps(const SequenceView<Letter>& seq) {
  auto copy = Sequence<Letter>(seq);
  removeGapsInPlace(copy);
  return copy;
}
