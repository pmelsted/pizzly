#ifndef COMMON_H
#define COMMON_H

#include <string>

#define PIZZLY_VERSION "0.37.1"

struct ProgramOptions {
  std::string gtf;
  std::string cache;
  std::string fusionFile;
  std::string fasta;
  std::string outprefix;
  int k;
  int alignScore;
  int insertSize;
  int kmerScore;
  ProgramOptions() : kmerScore(2), insertSize(400) {}
};

#endif // COMMON_H
