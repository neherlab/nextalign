#pragma once

template<typename InputContainer, typename OutputContainer, typename UnaryOperation>
OutputContainer map(const InputContainer& input, UnaryOperation op, OutputContainer) {
  OutputContainer result = {};
  result.reserve(input.size());
  std::transform(input.cbegin(), input.cend(), std::back_inserter(result), op);
  return result;
}
