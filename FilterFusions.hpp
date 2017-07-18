#ifndef FILTER_FUSIONS_HPP
#define FILTER_FUSIONS_HPP

#include "common.h"
#include "GeneModel.hpp"
#include <fstream>
#include <regex>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <string>
#include <functional>
#include <limits>

#include <seqan/align.h>
#include <seqan/align_split.h>
#include <seqan/index.h>
#include <seqan/store.h>


struct TranscriptPos {
  std::string tr;
  int pos;
  Strandedness strand;
  bool operator==(const TranscriptPos &o) const {
    return pos == o.pos && strand == o.strand && tr == o.tr;
  }
};

typedef std::unordered_map<std::string, std::unordered_map<std::string,std::vector<int>>> GeneGraph;



namespace std {
  template<> struct hash<TranscriptPos> {
    std::size_t operator()(const TranscriptPos &t) const {
      //return std::hash<std::string>{}(t.tr);
      return std::hash<int>{}(t.pos); // this is terrible
    }
  };
}

static std::regex trslist_re("\\(((\\w|\\.)+),(-?\\d+),(FW|RC)\\)");

void parseTrList(const std::string &trlist, std::vector<TranscriptPos> &trpos) {
  //std::regex tr_entry("([[:alnum:]]+,-?[[:digit:]]+,[FW|RC]\)")
  std::smatch m;
  auto sre_begin = std::sregex_iterator (trlist.begin(), trlist.end(), trslist_re);
  auto sre_end = std::sregex_iterator();
  for (std::sregex_iterator it = sre_begin; it != sre_end; ++it) {
    auto m = *it;
    TranscriptPos tp;
    tp.tr = m[1];
    tp.pos = std::stoi(m[3]);
    tp.strand = (m[4] == "FW") ? Strandedness::FORWARD : Strandedness::REVERSE;
    trpos.push_back(std::move(tp));
  }
}

// half open intervals [a,b[ and [c,d[
// returns 0 if overlap, otherwise distance between them
int interval_dist(int a, int b, int c, int d) {
  if (a <= c) {
    if (c < b) {
      return 0;
    } else {
      return c-b+1; // half open
    }
  } else {
    if (a < d) {
      return 0;
    } else {
      return a-d+1;
    }
  }
}

struct FusionRecord {
  enum TYPE {PAIR,SPLIT} type;
  std::string n1,n2,seq1,seq2;
  std::pair<int,int> kpos1,kpos2;
  std::pair<bool,int> splitat; // true if read 1 is split, k-mer based position
  std::vector<TranscriptPos> tr1,tr2;  
};

struct SplitRecord {
  enum TYPE {PAIR, SPLIT} type;
  std::string n1,n2;
  std::vector<TranscriptPos> tr1,tr2;
  int splitPos1, splitPos2;
  bool dir1,dir2; // direction of tr1, tr2 
  std::vector<std::pair<TranscriptPos, TranscriptPos>> sp1, sp2;
  std::unordered_set<std::string> g1, g2;
};


struct SplitInfo {
  std::string tr1,tr2;
  int pos1,pos2; // position of split
  bool dir1, dir2; // true for [:pos1+1], false for [pos1:] for tr1 tr2 
  bool strand1, strand2; // true for forward
  bool operator==(const SplitInfo &o) const {
    return pos1 == o.pos1 && pos2 == o.pos2 &&
           strand1 == o.strand1 && strand2 == o.strand2 &&
           tr1 == o.tr1 && tr2 == o.tr2 &&
           dir1 == o.dir1 && dir2 == o.dir2; 
  }
};

namespace std {
  template<> struct hash< SplitInfo> {
      std::size_t operator()(const SplitInfo &t) const {        
        return std::hash<int>{}(t.pos1) ^ std::hash<int>{}(t.pos2); // this is terrible
      }
    };
};

std::istream& operator>>(std::istream& stream, FusionRecord &rec)
{
  std::string tmp;
  if (stream >> tmp) {
    if (tmp == "PAIR") {
      rec.type = FusionRecord::PAIR;
    } else if (tmp == "SPLIT") {
      rec.type = FusionRecord::SPLIT;
    } 
    stream >> rec.n1 >> rec.seq1 >> tmp; // need to parse kpos1,2
    int cpos = tmp.find(',');
    rec.kpos1.first = std::stoi(tmp.substr(0,cpos));
    rec.kpos1.second = std::stoi(tmp.substr(cpos+1));
    stream >> rec.n2 >> rec.seq2 >> tmp;
    cpos = tmp.find(',');
    rec.kpos2.first = std::stoi(tmp.substr(0,cpos));
    rec.kpos2.second = std::stoi(tmp.substr(cpos+1));
    stream >> tmp; // parse splitat or NA
    if (tmp != "NA" && tmp.size() >= 13) {
      //splitat=(1,X)
      rec.splitat.first = (tmp[9] == '0');
      rec.splitat.second = std::stoi(tmp.substr(11,tmp.size()-12));
    }

    std::string tr1,tr2;
    stream >> tr1 >> tr2;
    parseTrList(tr1, rec.tr1);
    parseTrList(tr2, rec.tr2);
  } 
  return stream;
}

// returns true if the record should be filtered
bool filterFusionRecord(const Transcriptome &trx, ProgramOptions& options, const FusionRecord& rec, SplitRecord &splitrec) {

  std::unordered_set<std::string> g1,g2;

  auto get_gene_set = [&](const std::vector<TranscriptPos> &trlist, std::unordered_set<std::string> &gene_list) {
    for (auto &tr : trlist) {
      auto it = trx.trxToGeneId.find(tr.tr);
      if (it != trx.trxToGeneId.end()) {      
        auto git = trx.genes.find(it->second)->second;
        if (git.chr == "chrM" || git.chr == "M" || git.chr == "MT") {
          return true; // discard anything from the mitochondria, TODO: exclude unplaced regions
        }
        gene_list.insert(it->second);
      }    
    }
    return false;
  };

  if (get_gene_set(rec.tr1, g1)) {
    return true;
  }
  if (get_gene_set(rec.tr2, g2)) {
    return true;
  }
  
  for (auto &g : g1) {
    if (g2.count(g)) {
      return true;
    }
  }


  
  auto find_gene_range = [&](const std::unordered_set<std::string> glist) -> ExonModel {
    int limit = 1000;
    ExonModel m;
    m.start = std::numeric_limits<int>::max();
    m.stop = -1;
    


    for (auto &g : glist) {
      // check if first
      const GeneModel &gene = trx.genes.find(g)->second;
      if (gene.type != BioType::PROTEIN) {
        continue; // strict type checking
      }
      if (m.stop < 0) {
        m.chr = gene.chr;
        m.start = gene.start;
        m.stop = gene.stop;
      } else {
        if (m.chr != gene.chr) {
          m.stop = -1; m.start =0;
          return std::move(m);
        }
        int dist = interval_dist(m.start,m.stop, gene.start, gene.stop);
        if (dist <= limit) {
          m.start = std::min(m.start,gene.start);
          m.stop = std::max(m.stop, gene.stop);
        }
      } 
    }
    return std::move(m);
  };
  
  auto i1 = find_gene_range(g1);
  auto i2 = find_gene_range(g2);

  if (i1.stop < 0 || i2.stop < 0) {
    return true;
  }

  if (i1.chr == i2.chr) {
    int dist = interval_dist(i1.start, i1.stop, i2.start, i2.stop);
    if (dist < 10000) {
      return true;
    }
  }

  
  // see if the k-mer matching region of read 1 matches any of the transcripts of read 2
  typedef seqan::CharString TSeq;
  TSeq rs1(rec.seq1.substr(rec.kpos1.first, options.k + rec.kpos1.second-rec.kpos1.first).c_str());
  TSeq rs2(rec.seq2.substr(rec.kpos2.first, options.k + rec.kpos2.second-rec.kpos2.first).c_str());
  

  seqan::Align<TSeq, seqan::ArrayGaps> align;
  seqan::resize(rows(align), 2);

  auto mapKmersWithErrors = [&](TSeq& krseq, const std::vector<TranscriptPos> &trlist) -> bool {
    typedef seqan::Index<seqan::StringSet<TSeq>, seqan::IndexQGram<seqan::Shape<seqan::Dna, seqan::UngappedShape<8> >, seqan::OpenAddressing> > TIndex;
    typedef seqan::Pattern<TIndex, seqan::Swift<seqan::SwiftSemiGlobal> > TPattern;
    typedef seqan::Finder<const TSeq, seqan::Swift<seqan::SwiftSemiGlobal> > TFinder;
    seqan::StringSet<TSeq> reads;
    appendValue(reads,krseq);
    seqan::reverseComplement(krseq);
    appendValue(reads,krseq);
    TIndex index(reads);
    TPattern pattern(index);
    int matchlen = seqan::length(krseq);

    double epsilon = std::min(((double) options.kmerScore) / options.k,0.1);

    for (auto &tr : trlist) {
      const auto &it = trx.seqs.find(tr.tr);
      if (it != trx.seqs.end()) {      
        const TSeq &trseq = it->second;
        TFinder finder(trseq);
        while (find(finder, pattern, epsilon, options.k)) {
          seqan::Finder<const TSeq> verifyFinder(trseq);
          seqan::setPosition(verifyFinder, seqan::beginPosition(finder));
          seqan::Pattern<const TSeq, seqan::Myers<>> verifyPattern(reads[seqan::position(pattern).i1]);
          while (seqan::find(verifyFinder, verifyPattern, -((int) (epsilon*options.k))) 
              && seqan::position(verifyFinder) < seqan::endPosition(seqan::infix(finder))) {
            //std::cout << seqan::infix(finder)  << " -- " << tr.tr << std::endl;
            //std::cout << reads[seqan::position(pattern).i1] << " -- " << rec.n1 << ((rec.type == FusionRecord::PAIR) ? " PAIR" : " SPLIT") << std::endl;
            return true;
          }
        }
      }
    }
    return false;
  };
  
  // does any of the matcing k-mers map to one of the other transcript with errors?
  if (mapKmersWithErrors(rs1,rec.tr2)) {
    return true;
  }
  if (mapKmersWithErrors(rs2,rec.tr1)) {
    return true;
  }
  
  
  // does the entire read map to the respective transcripts?
  std::vector<TranscriptPos> filtered_tr1, filtered_tr2;

  auto getSequencesFromTrs = [&](const std::vector<TranscriptPos> &trs, int length, bool orient = true) -> std::unordered_map<std::string, std::vector<TranscriptPos>> {
    std::unordered_map<std::string, std::vector<TranscriptPos>> tlist;
    for (const TranscriptPos &trpos : trs) {
      auto git = trx.trxToGeneId.find(trpos.tr);
      if (git == trx.trxToGeneId.end()) {
        continue;
      }
      auto gene_id = git->second; //TODO refactor this out
      auto &genemodel = trx.genes.find(gene_id)->second;
      auto &trmodel = genemodel.transcripts.find(trpos.tr)->second;
      if (trmodel.type != BioType::PROTEIN) {
        continue; // modify
      }
      const auto& trseq = trx.seqs.find(trmodel.id)->second;
      int pos = trpos.pos-1; // 1-based like SAM
      int tr_start = (trpos.strand == Strandedness::FORWARD) ? pos : pos - length+1;
      int tr_stop = (trpos.strand == Strandedness::FORWARD) ? pos + length : pos+1; // off by one?

      // discard off bounds things
      if (tr_start < 0) {
        tr_start = 0;
      }

      if (tr_stop > seqan::length(trseq)) {
        tr_stop = seqan::length(trseq);
      }

      auto read_match = seqan::CharString(seqan::infix(trseq, tr_start, tr_stop));
      if ((trpos.strand == Strandedness::REVERSE) == orient) {
        seqan::reverseComplement(read_match);
      }

      std::string tseq(seqan::toCString(read_match));
      tlist[tseq].push_back(trpos);
    }
    return std::move(tlist);
  };

  auto mapReadsToTranscriptPos = [&](const std::string &seq, const std::vector<TranscriptPos> &trs, std::vector<TranscriptPos> &out_tr) {
    std::unordered_map<std::string, std::vector<TranscriptPos>> tlist;
    int readlen = seq.size();
    tlist = getSequencesFromTrs(trs, readlen);

    for (auto &it : tlist) {
      auto& tseq = it.first;
      auto& tposv = it.second;

      seqan::CharString read(seq.c_str());
      seqan::CharString ref(tseq.c_str()); // inefficient, refactor at some point

      int ascore = seqan::globalAlignmentScore(read, ref, seqan::MyersBitVector());
      if ( ascore >= -options.alignScore) {
        out_tr.insert(out_tr.end(), tposv.begin(), tposv.end());
        // debuggingstd::cout <<
        /*
        typedef seqan::Align<seqan::CharString, seqan::ArrayGaps> TAlign;
        TAlign align;
        seqan::resize(seqan::rows(align), 2);
        seqan::assignSource(seqan::row(align, 0), read);
        seqan::assignSource(seqan::row(align, 1), ref); 
        int score = seqan::globalAlignment(align, seqan::Score<int, seqan::Simple>(0, -1, -1));
        std::cout << "Score: " << score << std::endl;
        std::cout << align + std::endl;
        */
      }

    }
  };

  mapReadsToTranscriptPos(rec.seq1, rec.tr1, filtered_tr1);
  mapReadsToTranscriptPos(rec.seq2, rec.tr2, filtered_tr2);


  auto allSameStrand = [&](const std::vector<TranscriptPos> trp, Strandedness& s) -> bool {
    if (trp.empty()) {
      return false;
    }
    s = trp[0].strand;    
    for (const auto &tp : trp) {
      if (tp.strand != s) {
        return false;
      }
    }
    return true;
  };

  auto splitAlignReads = [&](const std::string &seq, const std::vector<TranscriptPos> &trsplit, const std::vector<TranscriptPos> &trsafe, std::vector<std::pair<TranscriptPos, TranscriptPos>> &out, int& splitPos) {
    // 
    splitPos = -1;
    std::unordered_map<std::string, std::vector<TranscriptPos>> trsplitlist, trsafelist;
    int readlen = seq.size();
    int fraglen = options.insertSize;
    if (fraglen < 2*readlen) {
      fraglen = 2*readlen;
    }
    trsplitlist = getSequencesFromTrs(trsplit, readlen);
    trsafelist = getSequencesFromTrs(trsafe, fraglen, false); // use false here?

    for (auto &tsplit : trsplitlist) {      
      for (auto &tsafe : trsafelist) {

        Strandedness splitStrand, safeStrand;
        bool allSameSplit = allSameStrand(tsplit.second, splitStrand);
        bool allSameSafe = allSameStrand(tsafe.second, safeStrand);

        if (!allSameSafe || !allSameSplit) {          
          //std::cout << "not all same direction " + std::endl;
          continue; // weird mapping, fix later
        }
        if (splitStrand == safeStrand) {
          //std::cout << "not same strand " + std::endl;
          continue; // weird potential fusion, TODO look into this
        }
        bool forwardRead = (splitStrand == Strandedness::FORWARD);
        // see if we can align across the sequences        
        seqan::CharString read(seq.c_str());
        seqan::CharString splitref(tsplit.first.c_str()); // inefficient, refactor at some point
        seqan::CharString saferef(tsafe.first.c_str());

        seqan::Gaps<seqan::CharString> readRowLeft, readRowRight, splitRowLeft, safeRowRight;
        seqan::setSource(readRowLeft, read);
        seqan::setSource(readRowRight, read);
        seqan::setSource(splitRowLeft, splitref);
        seqan::setSource(safeRowRight, saferef);

        seqan::Score<int, seqan::Simple> scoringScheme(0,-1,-1);
        seqan::AlignConfig<false, true, true, false> config;
        int splitScore = seqan::splitAlignment(readRowLeft, splitRowLeft, readRowRight, safeRowRight, scoringScheme, config);
        // debug
        /*
        std::cout << std::endl + "Score : " + splitScore + std::endl; 
        std::cout << "splitRowLEFT: " + splitRowLeft + std::endl;
        std::cout << "READLEFT    : " + readRowLeft + std::endl;
        std::cout << "safeRowRight: " + safeRowRight + std::endl;
        std::cout << "READRIGHT   : " + readRowRight + std::endl;
        */
        if (splitScore >= -options.alignScore) {
          if (splitPos < 0) {
            splitPos = seqan::toSourcePosition(readRowLeft, seqan::clippedEndPosition(readRowLeft));
          }
          for (const auto& t1 : tsplit.second) {
            for (const auto& t2 : tsafe.second) {
              TranscriptPos t2copy;
              t2copy.tr = t2.tr;
              t2copy.strand = (t2.strand == Strandedness::FORWARD) ? Strandedness::REVERSE : Strandedness::FORWARD; 
              int t2pos = seqan::toSourcePosition(safeRowRight, 0);
              int t2clip= seqan::toSourcePosition(readRowRight, 0); //seqan::clippedBeginPosition(readRowRight));
              if (forwardRead) {
                t2copy.pos = t2.pos - seqan::length(saferef) + t2pos + (seqan::length(splitref) - t2clip);
              } else {
                t2copy.pos = t2.pos + seqan::length(saferef) - t2pos -1;
              }
              out.push_back({t1, t2copy});
            }
          }
        }
      }
    }    
  };



  std::vector<std::pair<TranscriptPos, TranscriptPos>> split_reads1, split_reads2;

  
  if (filtered_tr1.empty() && filtered_tr2.empty()) {
    int pos1,pos2;
    splitAlignReads(rec.seq1, rec.tr1, rec.tr2, split_reads1, pos1);
    if (split_reads1.empty()) {
      return true;
    } 
    splitAlignReads(rec.seq2, rec.tr2, rec.tr1, split_reads2, pos2);
    if (split_reads2.empty()) {
      return true;
    }
    // both records are split
    splitrec.type = SplitRecord::SPLIT;
    splitrec.n1 = rec.n1;
    splitrec.n2 = rec.n2;

    {
      std::unordered_set<TranscriptPos> tmp;
      for (const auto &p : split_reads1) {
        tmp.insert(p.first);
      }
      for (auto &t : tmp) {
        splitrec.tr1.push_back(std::move(t));
      }
      tmp.clear();
      for (const auto &p : split_reads2) {
        tmp.insert(p.first);
      }
      for (auto &t : tmp) {
        splitrec.tr2.push_back(std::move(t));
      }
    }    

    get_gene_set(splitrec.tr1, splitrec.g1);
    get_gene_set(splitrec.tr2, splitrec.g2);
    
    splitrec.splitPos1 = pos1;
    splitrec.splitPos2 = pos2;
    splitrec.dir1 = split_reads1[0].first.strand == Strandedness::FORWARD;
    splitrec.dir2 = !splitrec.dir2; 
    if ((split_reads2[0].first.strand == Strandedness::FORWARD) != splitrec.dir1) {
      //std::cerr << "reads do not agree on direction!" << std::endl;
    }
    splitrec.sp1 = std::move(split_reads1);
    splitrec.sp2 = std::move(split_reads2); 
  } else if (filtered_tr1.empty() && !filtered_tr2.empty()) {
    int pos1;
    splitAlignReads(rec.seq1, rec.tr1, rec.tr2, split_reads1, pos1);
    if (split_reads1.empty()) {
      return true;
    }
    splitrec.type = SplitRecord::SPLIT;
    splitrec.n1 = rec.n1;
    splitrec.n2 = rec.n2;    
    {
      std::unordered_set<TranscriptPos> tmp;
      for (const auto &p : split_reads1) {
        tmp.insert(p.first);
      }
      for (auto &t : tmp) {
        splitrec.tr1.push_back(std::move(t));
      }
    }
    splitrec.tr2 = std::move(filtered_tr2);
    get_gene_set(splitrec.tr1, splitrec.g1);
    get_gene_set(splitrec.tr2, splitrec.g2);
    splitrec.splitPos1 = pos1;
    splitrec.dir1 = split_reads1[0].first.strand == Strandedness::FORWARD;
    splitrec.dir2 = !splitrec.dir1; 
    splitrec.splitPos2 = -1;
    splitrec.sp1 = std::move(split_reads1); 
  } else if (!filtered_tr1.empty() && filtered_tr2.empty()) {
    // splitalign read 2
    int pos2;
    splitAlignReads(rec.seq2, rec.tr2, rec.tr1, split_reads2, pos2);
    if (split_reads2.empty()) {
      return true;
    }
    splitrec.type = SplitRecord::SPLIT;
    splitrec.n1 = rec.n1;
    splitrec.n2 = rec.n2;
    splitrec.tr1 = std::move(filtered_tr1);
    {
      std::unordered_set<TranscriptPos> tmp;
      for (const auto &p : split_reads2) {
        tmp.insert(p.first);
      }
      for (auto &t : tmp) {
        splitrec.tr2.push_back(std::move(t));
      }
    }
    get_gene_set(splitrec.tr1, splitrec.g1);
    get_gene_set(splitrec.tr2, splitrec.g2);
    splitrec.splitPos1 = -1;
    splitrec.splitPos2 = pos2;
    splitrec.dir1 = split_reads2[0].first.strand == Strandedness::FORWARD;
    splitrec.dir2 = !splitrec.dir2;
    // check if other read agrees 
    splitrec.sp2 = std::move(split_reads2); 
  } else {
    // both are paired.

    Strandedness st1,st2;
    bool allSameRead1 = allSameStrand(filtered_tr1, st1);
    bool allSameRead2 = allSameStrand(filtered_tr2, st2);

    if (!allSameRead1 && !allSameRead2) {
      //std::cout << "not all same direction, can't fix'" << std::endl;
      return true; // weird mapping
    } else if (!allSameRead1 && allSameRead2) {
      // fix read 1, discard discordant pairs      
      std::copy_if(filtered_tr1.begin(), filtered_tr1.end(), 
                   std::back_inserter(splitrec.tr1), 
                   [&](const auto &x) {
                     return x.strand != st2;
                   });      
      if (splitrec.tr1.empty()) {

        return true;
      }
      splitrec.tr2 = std::move(filtered_tr2);

    } else if (allSameRead1 && !allSameRead2) {
      std::copy_if(filtered_tr2.begin(), filtered_tr2.end(), 
                        std::back_inserter(splitrec.tr2), 
                        [&](const auto &x) {
                          return x.strand != st1;
                        });      
      if (splitrec.tr2.empty()) {
        return true;
      }
      splitrec.tr1 = std::move(filtered_tr1);
    } else {
      if (st1 == st2) {
        //std::cout << "not same strand" << std::endl;
        return true; // weird potential fusion, TODO look into this
      }
      splitrec.tr1 = std::move(filtered_tr1);
      splitrec.tr2 = std::move(filtered_tr2);   
    }
    

    splitrec.type = SplitRecord::PAIR;
    splitrec.n1 = rec.n1;
    splitrec.n2 = rec.n2;
    
    get_gene_set(splitrec.tr1, splitrec.g1);
    get_gene_set(splitrec.tr2, splitrec.g2);
    splitrec.splitPos1 = -1;
    splitrec.splitPos2 = -1;
    splitrec.dir1 = st1 == Strandedness::FORWARD; // check read 1
    splitrec.dir2 = st2 == Strandedness::FORWARD; // check read 2
  }

  return false;
}

void processFusions(const Transcriptome &trx, ProgramOptions& options) {
  // open fusion file and parse
  std::ifstream in(options.fusionFile);

  std::string line;
  std::vector<std::pair<FusionRecord, SplitRecord>> fusions;
  int n = 0;
  int count=0;

  while (std::getline(in,line)) {
    if (line.size() >= 4 && line.substr(0,4) == "TYPE") {
      //std::cout << line << "\n";
      continue;
    }
    
    count++;
    std::istringstream iss(line);
    FusionRecord rec;
    SplitRecord splitrec;
    iss >> rec;
    if (!filterFusionRecord(trx,options,rec, splitrec)) {
      //std::cout << line << "\n";
      n++;
      fusions.push_back({std::move(rec), std::move(splitrec)});
    }
  }
  std::cout.flush();

  // create gene graph
  GeneGraph G;
  int fi = 0;
  for (const auto& p : fusions) {
    for (const auto &g1 : p.second.g1) {
      for (const auto &g2 : p.second.g2) {
        /*G[g1].insert({g2,fi});
        G[g2].insert({g1,fi});*/
        G[g1][g2].push_back(fi);
        G[g2][g1].push_back(fi);
      }
    }
    ++fi;
  }
  

  using SplitInfoMap = std::unordered_map<SplitInfo, int>;

  auto findTranscriptSplit = [&](const std::pair<TranscriptPos, TranscriptPos> &sp, int splitPos, bool dir, int readlen, SplitInfoMap &out) {
    if (readlen - splitPos < std::max(8, 2*options.alignScore)) {
      return;
    }    
    SplitInfo spi;

    spi.tr1 = sp.first.tr;
    spi.tr2 = sp.second.tr;    
    spi.strand1 = sp.first.strand == Strandedness::FORWARD;
    spi.strand2 = sp.second.strand == Strandedness::FORWARD;
    

    // set direction as well.
    if (spi.strand1) {
      spi.pos1 = sp.first.pos + (splitPos-1);
    } else {
      spi.pos1 = sp.first.pos - splitPos;
    }
    if (spi.strand2) {
      spi.pos2 = sp.second.pos - (readlen - splitPos);
    } else {
      spi.pos2 = sp.second.pos; //+ (readlen - splitPos);
    }

    // reverse complement of 
    if (!spi.strand1 && !spi.strand2) {
      std::swap(spi.tr1,spi.tr2);
      std::swap(spi.pos1,spi.pos2);
      spi.strand1 = !spi.strand1;
      spi.strand2 = !spi.strand2;
    }

    spi.dir1 = true;
    spi.dir2 = false;

    auto it = out.find(spi);
    if (it == out.end()) {
      out.insert({std::move(spi),1});
    } else {
      it->second++;
    }
  };



  auto mapToClosestExon = [&](const std::string &tr, int pos, bool forward) -> int {
    const auto git = trx.trxToGeneId.find(tr);
    int mindist = pos;
    if (git != trx.trxToGeneId.end()) {
      const std::string gid = git->second;
      const auto &gm = (trx.genes.find(gid))->second;
      const auto &trm = (gm.transcripts.find(tr))->second;
      int a = 0,b=0;
      
      for (int i = 0; i < trm.exons.size(); i++) {        
        b += (trm.exons[i].stop - trm.exons[i].start);
        // exon is from [a,b)
        if (a <= pos  && pos < b) {
          if ((pos-a) < abs(mindist)) {
            mindist = pos-a;
          }  
          if ((b - pos) < abs(mindist)) {
            mindist = pos - b;
          }
          //mindist = std::min(std::min(pos-a, b - pos), mindist);
        }
        a = b;
      }
    }
    return mindist; 
  };

  auto mapToForwardExon = [&](const std::string &tr, int pos) -> std::pair<int,int> {
    const auto git = trx.trxToGeneId.find(tr);

    int minPosDist = std::numeric_limits<int>::max();
    int minNegDist = 0;
    if (git != trx.trxToGeneId.end()) {
      const std::string gid = git->second;
      const auto &gm = (trx.genes.find(gid))->second;
      const auto &trm = (gm.transcripts.find(tr))->second;
      int a = 0,b=0;
      
      for (int i = 0; i < trm.exons.size(); i++) {        
        b += (trm.exons[i].stop - trm.exons[i].start);
        // exon is from [a,b)
        if (a <= pos  && pos < b) {
          if ((a-pos) < minNegDist) {
            minNegDist = a-pos;
          }  
          if ((b - pos) < minPosDist) {
            minPosDist = b - pos;
          }          
        }
        a = b;
      }
    }
    return std::make_pair(minNegDist, minPosDist); 
  };

  
  auto snapToJunction = [&](const SplitInfoMap &splits, int snapDist) -> SplitInfoMap {
    SplitInfoMap ret;
    for (const auto &spp : splits) {
      auto sp = spp.first;      
      int split_reads = spp.second;

      int d1 = mapToClosestExon(sp.tr1,sp.pos1, sp.strand1);
      int d2 = mapToClosestExon(sp.tr2,sp.pos2, sp.strand2);

      int gap = d2 - d1;
      
      if (std::abs(d1) <= snapDist && std::abs(d2) <= snapDist) {
        SplitInfo spcopy = sp;
        spcopy.pos1 = sp.pos1 - d1;
        spcopy.pos2 = sp.pos2 - d2;
        auto it = ret.find(spcopy);
        if (it != ret.end()) {
          it->second += split_reads;
        } else {
          ret.insert({spcopy,split_reads});
        }
      } else {
        ret.insert(spp);
      }
    }
    return std::move(ret);
  };


  auto findTranscriptInVector = [&](const std::vector<TranscriptPos> v, std::string t) {
    for (const auto &x : v) {
      if (x.tr == t) {
        return true;
      }
    }
    return false;
  };


  auto writeFusions = [&](std::ofstream& jsonOut, std::ofstream& fastaOut, bool filter) {
    // write json headers and info
    jsonOut << "{\n  \"genes\" : [\n";
    bool firstJsonCommaGenes = true;
    std::vector<int> readsInTr, readsInGene;
    readsInTr.reserve(1000);
    readsInGene.reserve(1000);
    // write output
    for (const auto &gp1 : G) {
      const auto &g1 = gp1.first;
      const auto &gm1 = trx.genes.find(g1); 
      const auto &Gg1 = gp1.second;
      for (const auto &gg2 : Gg1) {
        const auto& g2 = gg2.first;
        if (g1 < g2) {

          readsInGene.clear();
          const auto &gm2 = trx.genes.find(g2);
          const auto &v = gg2.second;
          int paircount = 0;
          int splitcount = 0;
          SplitInfoMap splits;
          for (int i : v) {
            const FusionRecord& fr = fusions[i].first;
            const SplitRecord& sr = fusions[i].second;
            if (sr.type == SplitRecord::PAIR) {
              paircount++;
            } else if (sr.type == SplitRecord::SPLIT) {
              splitcount++;
            }

            if (sr.splitPos1 >= 0) {
              int readlen = fr.seq1.size();
              for (const auto &sp : sr.sp1) {
                findTranscriptSplit(sp, sr.splitPos1, sr.dir1, readlen, splits);
              }
            }
            if (sr.splitPos2 >= 0) {
              int readlen = fr.seq2.size();
              for (const auto &sp : sr.sp2) {
                findTranscriptSplit(sp, sr.splitPos2, sr.dir2, readlen, splits);
              }
            } 
          }

          bool hasAnyTr = false;
          auto writeGeneInfoToJsonStream = [&](bool swap) {
            if (hasAnyTr) {
              return;
            }
            if (!firstJsonCommaGenes) {
              jsonOut <<  ",\n"; 
            }
            firstJsonCommaGenes = false;
            if (!swap) {
              jsonOut << "    {\n      \"geneA\" : { \"id\" : \"" << g1 << "\", \"name\" : \""<< gm1->second.name <<"\"},\n"
                      <<        "      \"geneB\" : { \"id\" : \"" << g2 << "\", \"name\" : \""<< gm2->second.name << "\"},\n";
            } else {
              jsonOut << "    {\n      \"geneA\" : { \"id\" : \"" << g2 << "\", \"name\" : \""<< gm2->second.name <<"\"},\n"
                      <<        "      \"geneB\" : { \"id\" : \"" << g1 << "\", \"name\" : \""<< gm1->second.name << "\"},\n";
            }
            jsonOut <<        "      \"paircount\" : " << std::to_string(paircount) << ",\n      \"splitcount\" : " << std::to_string(splitcount) << ",\n" 
                    <<        "      \"transcripts\" : [\n";
            hasAnyTr = true;
          };

          std::unordered_set<std::string> fnames;

          bool firstJsonCommaTrans = true;
          SplitInfoMap snappedSplits = snapToJunction(splits, 4);      
          for (auto & spp : snappedSplits) {
            auto &sp = spp.first;          
            int split_reads = spp.second;
            int ed1 = mapToClosestExon(sp.tr1,sp.pos1, sp.strand1);
            int ed2 = mapToClosestExon(sp.tr2,sp.pos2, sp.strand2);
            int tp1 = sp.pos1 - ed1;
            int tp2 = sp.pos2 - ed2;
            const auto & seq1 = trx.seqs.find(sp.tr1)->second;
            const auto & seq2 = trx.seqs.find(sp.tr2)->second;
            
            
            readsInTr.clear();
            int t_count = 0;

            if ((gm1->second.transcripts.count(sp.tr1) > 0 && gm2->second.transcripts.count(sp.tr2) > 0)
            ||  (gm1->second.transcripts.count(sp.tr2) > 0 && gm2->second.transcripts.count(sp.tr1) > 0)) {
              
              
              for (int i = 0; i < v.size(); i++) {              
                if ((findTranscriptInVector(fusions[v[i]].second.tr1, sp.tr1) && findTranscriptInVector(fusions[v[i]].second.tr2, sp.tr2))
                 || (findTranscriptInVector(fusions[v[i]].second.tr1, sp.tr2) && findTranscriptInVector(fusions[v[i]].second.tr2, sp.tr1))) {
                  t_count++;
                  readsInTr.push_back(i);
                }
              }
              
              if (!filter || ((t_count >= 1) && std::abs(ed1) <= 10 && std::abs(ed2) <= 10)) {
                bool swap = false;
                auto tg1 = trx.trxToGeneId.find(sp.tr1);
                auto tg2 = trx.trxToGeneId.find(sp.tr2);
                assert(tg1 != trx.trxToGeneId.end());
                assert(tg2 != trx.trxToGeneId.end());
                if (g1 != tg1->second) {
                  swap = true;
                  assert(g1 == tg2->second);
                  assert(g2 == tg1->second);
                } else {
                  assert(g1 == tg1->second);
                  assert(g2 == tg2->second);
                }
                writeGeneInfoToJsonStream(swap);
                
                std::string fasta_name = sp.tr1 + "_0:" + std::to_string(tp1) + "_" + sp.tr2 + "_" + std::to_string(tp2) + ":" + std::to_string(seqan::length(seq2));
                if (fnames.count(fasta_name) > 0) {
                  continue;
                } else {
                  fnames.insert(fasta_name);
                }
               
                if (!firstJsonCommaTrans) {
                  jsonOut << ",\n";
                }
                firstJsonCommaTrans = false;
                jsonOut <<"        {\n          \"fasta_record\": \"" << fasta_name <<  "\",\n"
                        <<           "          \"transcriptA\": {\"id\" : \"" << sp.tr1 << "\", \"startPos\" : " << std::to_string((sp.dir1) ? 0 : tp1)
                        << ", \"endPos\" : " << std::to_string((sp.dir1) ? tp1 : seqan::length(seq1)) << ", \"edit\" : " << std::to_string(ed1) 
                        << ", \"strand\" : " << ((sp.strand1) ? "true" : "false") << "},\n";
                        
                jsonOut <<           "          \"transcriptB\": {\"id\" : \"" << sp.tr2 << "\", \"startPos\" : " << std::to_string((sp.dir2) ? 0 : tp2)
                        << ", \"endPos\" : " << std::to_string((sp.dir2) ? tp2 : seqan::length(seq2)) << ", \"edit\" : " << std::to_string(ed2)
                        << ", \"strand\" : " << ((sp.strand2) ? "true" : "false") << "},\n"
                        << "          \"support\" : " << std::to_string(t_count) << ",\n"
                        << "          \"reads\" : [";
                
                for (int i = 0; i < readsInTr.size(); i++) {
                  if (i > 0) {
                    jsonOut << ", ";
                  }
                  jsonOut << i;
                }                
                jsonOut <<  "]\n";
                jsonOut << "        }";
                std::copy(readsInTr.begin(), readsInTr.end(), std::back_inserter(readsInGene));

                fastaOut << ">" << fasta_name << "\n"
                         << seqan::prefix(seq1, tp1)
                         << seqan::suffix(seq2, tp2)
                         << "\n"; 


              }    
                        
            }

            //std::cout << std::endl; 
          }
          if (!hasAnyTr) {
            
            std::unordered_map<SplitInfo, std::vector<std::pair<int,int>>> trfmap;

            for (int i : v) {
              const FusionRecord& fr = fusions[i].first;
              const SplitRecord& sr = fusions[i].second;
              if (sr.type == SplitRecord::PAIR) {
                for (const auto &t1 : sr.tr1) {
                  for (const auto &t2 : sr.tr2) {
                    int rlen1 = seqan::length(fr.seq1);
                    int rlen2 = seqan::length(fr.seq2); 

                    const auto & seq1 = trx.seqs.find(t1.tr)->second;
                    const auto & seq2 = trx.seqs.find(t2.tr)->second;

                    bool t1strand = t1.strand == Strandedness::FORWARD;
                    bool t2strand = t2.strand == Strandedness::FORWARD;

                    // position of split
                    int tp1 = (t1strand) ? std::min(t1.pos+rlen1, (int)seqan::length(seq1)-1) : std::max(t1.pos-rlen1,0);
                    int tp2 = (t2strand) ? std::min(t2.pos+rlen2, (int)seqan::length(seq2)-1) : std::max(t2.pos-rlen2,0);

                    // distances to exon boundaries
                    auto ped1 = mapToForwardExon(t1.tr, tp1);
                    auto ped2 = mapToForwardExon(t2.tr, tp2);

                    // pick most reliable
                    int ed1 = 0;
                    int ed2 = 0;
                    int insLen = 0;
                    if (t1strand) {
                      ed1 = (ped1.first >= -10 && ped1.second > abs(ped1.first)) ? ped1.first : ped1.second;
                      insLen += rlen1 + ed1;
                    } else {
                      ed1 = (ped1.second <= 10 && abs(ped1.first) > ped1.second) ? ped1.second : ped1.first;
                      insLen += rlen1 - ed1;
                    }
                    if (t2strand) {
                      ed2 = (ped2.first >= -10 && ped2.second > abs(ped2.first)) ? ped2.first : ped2.second;
                      insLen += rlen2 + ed2;
                    } else {
                      ed2 = (ped2.second <= 10 && abs(ped2.first) > ped2.second) ? ped2.second : ped2.first;
                      insLen += rlen2 - ed2;
                    }

                    tp1 += ed1;
                    tp2 += ed2;    

                    assert(t1strand != t2strand);                

                    if ((gm1->second.transcripts.count(t1.tr) > 0 && gm2->second.transcripts.count(t2.tr) > 0)
                    ||  (gm1->second.transcripts.count(t2.tr) > 0 && gm2->second.transcripts.count(t1.tr) > 0)) {

                      SplitInfo spi;
                      if (t1strand) {
                        spi.tr1 = t1.tr;
                        spi.pos1 = tp1;
                        spi.dir1 = true;
                        spi.strand1 = Strandedness::FORWARD;
                        spi.tr2 = t2.tr;
                        spi.pos2 = tp2;
                        spi.dir2 = false;
                        spi.strand2 = Strandedness::FORWARD;
                      } else {
                        spi.tr1 = t2.tr;
                        spi.pos1 = tp2;
                        spi.dir1 = true;
                        spi.strand1 = Strandedness::FORWARD;
                        spi.tr2 = t1.tr;
                        spi.pos2 = tp1;
                        spi.dir2 = false;
                        spi.strand2 = Strandedness::FORWARD;
                      }

                      trfmap[spi].push_back({i,insLen});
                    }
                  }
                }
              }
            }

            for (const auto p : trfmap) {
              const SplitInfo& spi = p.first;
              const auto rv = p.second;
              
              const auto & seq1 = trx.seqs.find(spi.tr1)->second;
              const auto & seq2 = trx.seqs.find(spi.tr2)->second;
              
              readsInTr.clear();

              for (auto x : rv) {
                int vi = x.first;
                int insLen = x.second;
                if ((findTranscriptInVector(fusions[vi].second.tr1, spi.tr1) && findTranscriptInVector(fusions[vi].second.tr2, spi.tr2))
                         || (findTranscriptInVector(fusions[vi].second.tr1, spi.tr2) && findTranscriptInVector(fusions[vi].second.tr2, spi.tr1))) {                          
                  if (!filter || (insLen <= options.insertSize)) {
                    readsInTr.push_back(vi);
                  }
                }
              }
              int t_count = readsInTr.size();
              if (t_count == 0) {
                continue;
              }

              if (!filter || (t_count >= 2)) {
                bool swap = false;
                auto tg1 = trx.trxToGeneId.find(spi.tr1);
                auto tg2 = trx.trxToGeneId.find(spi.tr2);
                assert(tg1 != trx.trxToGeneId.end());
                assert(tg2 != trx.trxToGeneId.end());
                if (g1 != tg1->second) {
                  swap = true;
                  assert(g1 == tg2->second);
                  assert(g2 == tg1->second);
                } else {
                  assert(g1 == tg1->second);
                  assert(g2 == tg2->second);
                }
                writeGeneInfoToJsonStream(swap);

                if (!firstJsonCommaTrans) {
                  jsonOut << ",\n";
                }
                firstJsonCommaTrans = false;

                std::string fasta_name;
                
                fasta_name = spi.tr1 + "_0:" + std::to_string(spi.pos1) + "_" + spi.tr2 + "_" + std::to_string(spi.pos2) + ":" + std::to_string(seqan::length(seq2));
                jsonOut << "        {\n          \"fasta_record\": \"" << fasta_name << "\",\n"
                        <<            "          \"transcriptA\": {\"id\" : \"" << spi.tr1 << "\", \"startPos\" : " << std::to_string((spi.dir1) ? 0 : spi.pos1)
                        << ", \"endPos\" : " << std::to_string((spi.dir1) ? spi.pos1 : seqan::length(seq1))
                        << ", \"strand\" : " << ((spi.strand1 == Strandedness::FORWARD) ? "true" : "false") << "},\n";
                jsonOut <<           "          \"transcriptB\": {\"id\" : \"" << spi.tr2 << "\", \"startPos\" : " << std::to_string((spi.dir2) ? 0 : spi.pos2)
                        << ", \"endPos\" : " << std::to_string((spi.dir2) ? spi.pos2 : seqan::length(seq2))
                        << ", \"strand\" : " << ((spi.strand2 == Strandedness::FORWARD) ? "true" : "false") << "},\n"
                        << "          \"support\" : " << std::to_string(t_count) << ",\n"
                        << "          \"reads\" : [";
              
                for (int i = 0; i < readsInTr.size(); i++) {
                  if (i > 0) {
                    jsonOut << ", ";
                  }
                  jsonOut << i;
                }                
                jsonOut <<  "]\n        }";

                std::copy(readsInTr.begin(), readsInTr.end(), std::back_inserter(readsInGene));
                
                fastaOut << ">" << fasta_name << "\n"
                        << seqan::prefix(seq1, spi.pos1)
                        << seqan::suffix(seq2, spi.pos2)
                        << "\n"; 
              }
            }               
          }

          // write out supporting reads
          if (hasAnyTr) {
            jsonOut << "\n      ],\n"; // closing transcripts
            bool firstJsonRead = true;
            jsonOut <<  "      \"readpairs\" : [\n";
            /*   */
            for (int i = 0; i < v.size(); i++) {
              const auto& fr = fusions[v[i]].first;
              const auto& sr = fusions[v[i]].second;
              if (i > 0) {
                jsonOut <<",\n";
              }
              jsonOut << "        {\n" 
                      << "          \"type\" : \"" << ((sr.type == SplitRecord::SPLIT) ? "SPLIT" : "PAIR" ) << "\",\n"
                      << "          \"read1\" : { \"name\" : \"" << fr.n1 << "\", \"seq\" : \"" << fr.seq1 
                      << "\", \"splitpos\" : " << std::to_string(sr.splitPos1) << ", \"direction\" : \"" << ((sr.dir1)? "true" : "false") << "\", " 
                      << "\"kmerpos\" : { \"start\" : " << std::to_string(fr.kpos1.first) << ", \"stop\" : " << std::to_string(fr.kpos1.second) << "}},\n"
                      << "          \"read2\" : { \"name\" : \"" << fr.n2 << "\", \"seq\" : \"" << fr.seq2 
                      << "\", \"splitpos\" : " << std::to_string(sr.splitPos2) << ", \"direction\" : \"" << ((sr.dir2)? "true" : "false") << "\", "
                      << "\"kmerpos\" : { \"start\" : " << std::to_string(fr.kpos2.first) << ", \"stop\" : " << std::to_string(fr.kpos2.second) << "}}\n"
                      << "        }";
            }
            jsonOut << "\n      ]\n    }";
          }
        }
      }
    }
    jsonOut << "\n  ]\n}";
  };
  
  
  // write unfiltered files
  std::ofstream unFilteredFastaOut(options.outprefix + ".unfiltered.fusions.fasta");
  std::ofstream unFilteredJsonOut(options.outprefix + ".unfiltered.json");

  writeFusions(unFilteredJsonOut, unFilteredFastaOut, false);
  unFilteredFastaOut.close();
  unFilteredJsonOut.close();

  // write filtered fiels
  std::ofstream fastaOut(options.outprefix + ".fusions.fasta");
  std::ofstream jsonOut(options.outprefix + ".json");

  writeFusions(jsonOut, fastaOut, true);
  fastaOut.close();
  jsonOut.close();

  std::cerr << "Number of kept records " << n  << " out of " << count << std::endl; 

}

#endif // FILTER_FUSIONS_HPP