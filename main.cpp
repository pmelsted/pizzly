#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>


#include "common.h"

#include "GeneModel.hpp"
#include "FilterFusions.hpp"

seqan::ArgumentParser::ParseResult
parseCommandLine(ProgramOptions & options, int argc, char const ** argv) {

  seqan::ArgumentParser parser("pizzly");
  // We require one argument.
  seqan::addArgument(parser, seqan::ArgParseArgument(
      seqan::ArgParseArgument::STRING, "FUSION"));

  // Define Options
  seqan::addOption(parser, seqan::ArgParseOption(
      "k", "", "k-mer size used in kallisto",
      seqan::ArgParseArgument::INTEGER, "K"));
  seqan::addOption(parser, seqan::ArgParseOption(
      "a", "align-score", "Maximum number of mismatches allowed (default: 2)",
      seqan::ArgParseArgument::INTEGER, "ALIGN_SCORE"));
  seqan::addOption(parser, seqan::ArgParseOption(
      "i", "insert-size", "Maximum fragment size of library (default: 400)",
      seqan::ArgParseArgument::INTEGER, "INSERT_SIZE"));
  seqan::addOption(parser, seqan::ArgParseOption(
      "o", "output", "Prefix for output files",
      seqan::ArgParseArgument::STRING, "OUTPUT_PREFIX"));
  seqan::addOption(parser, seqan::ArgParseOption(
      "G", "gtf", "Annotation in GTF format",
      seqan::ArgParseArgument::STRING, "GTF"));
  seqan::addOption(parser, seqan::ArgParseOption(
      "C", "cache", "File for caching annotation (created if not present, otherwise reused from previous runs)",
      seqan::ArgParseArgument::STRING, "cache"));
  seqan::addOption(parser, seqan::ArgParseOption(
      "F", "fasta", "Fasta reference",
      seqan::ArgParseArgument::STRING, "FASTA"));

  seqan::setRequired(parser, "k");
  seqan::setRequired(parser, "o");
  seqan::setRequired(parser, "F");
  seqan::setVersion(parser, std::string(PIZZLY_VERSION));  

  // Parse command line.
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

  // Only extract  options if the program will continue after parseCommandLine()
  if (res != seqan::ArgumentParser::PARSE_OK)
      return res;

  // Extract option values.    
  getOptionValue(options.k, parser, "k");
  getOptionValue(options.alignScore, parser, "align-score");
  getOptionValue(options.insertSize, parser, "insert-size");
  getOptionValue(options.gtf, parser, "gtf");
  getOptionValue(options.fasta, parser, "fasta");
  getOptionValue(options.cache, parser, "cache");
  getOptionValue(options.outprefix, parser, "output");
  getArgumentValue(options.fusionFile, parser, 0);


  return seqan::ArgumentParser::PARSE_OK;

}

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ProgramOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // parse GTF file
    Transcriptome trx;
    if (options.cache.empty()) {
        parseGTF(trx, options.gtf);
    } else {
        std::ifstream in(options.cache);
        if (in.good()) {
            std::cerr << "Opening cached file ... ";
            loadTranscriptome(trx,in);
            std::cerr << "loaded " << trx.genes.size() << " genes and " << trx.trxToGeneId.size() << " transcripts" << std::endl;
            in.close();            
        } else {
            parseGTF(trx, options.gtf);
            in.close();
            std::ofstream out(options.cache);
            writeTranscriptome(trx,out);
            out.close();
        }
        in.close();
    }

    parseFasta(trx, options.fasta);

    

    std::cerr << "Read a total of " << trx.seqs.size() << " transcripts" << std::endl;

    // filter fusion edges
    processFusions(trx,options);


    return 0;
}
