#pragma once


#include <string>
#include <type_safe/constrained_type.hpp>

namespace ts = type_safe;

using non_empty_string = ts::constrained_type<std::string, ts::constraints::non_empty>;
