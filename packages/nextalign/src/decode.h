#pragma once


#include <string_view>

#include "translate.h"

constexpr const Aminoacid AMINOACID_UNKNOWN = 'X';
constexpr const Aminoacid AMINOACID_GAP = '-';
constexpr const Aminoacid AMINOACID_START = '#';
constexpr const Aminoacid AMINOACID_STOP = '*';

Aminoacid decode(const std::string_view& codon);

std::string encode(Aminoacid aa, const std::string_view& hint);
