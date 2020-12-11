#pragma once

#include <translate.h>

#include <string_view>

constexpr const Aminoacid AMINOACID_UNKNOWN = 'X';
constexpr const Aminoacid AMINOACID_GAP = '-';
constexpr const Aminoacid AMINOACID_STOP = '*';

Aminoacid decode(const std::string_view& codon);
