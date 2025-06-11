#ifndef COMMON_H
#define COMMON_H

#include <string>

#define PIZZLY_VERSION "0.45.0"

struct ProgramOptions {
  std::string gtf;
  std::string cache;
  std::string fusionFile;
  std::string fasta;
  std::string outprefix;
  bool ignoreProtein;
  bool longreads;
  bool paired;
  int k;
  int alignScore;
  int insertSize;
  int kmerScore;
  ProgramOptions() : kmerScore(2), insertSize(400), ignoreProtein(false), paired(true), longreads(false) {}
};

#endif // COMMON_H
