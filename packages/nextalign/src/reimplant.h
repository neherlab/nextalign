#pragma once

struct Alignment;
struct CodonAlignmentResult;
struct Gene;

Alignment reimplant(Alignment& alignmentImproved, const CodonAlignmentResult& codonAlignmentResult, const Gene& gene);
