#ifndef COMMON_H
#define COMMON_H

#include <string>

#define PIZZLY_VERSION "0.37.3"

struct ProgramOptions {
  std::string gtf;
  std::string cache;
  std::string fusionFile;
  std::string fasta;
  std::string outprefix;
  bool ignoreProtein;
  int k;
  int alignScore;
  int insertSize;
  int kmerScore;
  ProgramOptions() : kmerScore(2), insertSize(400), ignoreProtein(false) {}
};

#endif // COMMON_H
