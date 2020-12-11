#pragma once

#include <string>
#include <string_view>

using Aminoacid = char;
using Peptide = std::basic_string<Aminoacid>;

Peptide translate(const std::string_view& seq);
