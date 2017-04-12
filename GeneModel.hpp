#ifndef GENE_MODEL_HPP
#define GENE_MODEL_HPP

#include <vector>
#include <string>
#include <unordered_map>
#include <seqan/gff_io.h>

enum Strandedness {FORWARD, REVERSE, UNKNOWN};
enum BioType {PROTEIN, PSEUDO, OTHER};

Strandedness operator~(const Strandedness &s) {
  switch (s) {
    case FORWARD: return REVERSE; break;
    case REVERSE: return FORWARD; break;
    case UNKNOWN: return REVERSE; break;
  }
}

std::string typeToString(BioType t) {
  switch(t) {
    case BioType::PROTEIN:
      return "PROTEIN"; break;
    case BioType::PSEUDO:
      return "PSEUDO"; break;
    default: 
      return "OTHER"; break;
  }
}

char strandToChar(Strandedness s) {
  switch(s) {
    case Strandedness::FORWARD:
      return '+'; break;
    case Strandedness::REVERSE:
      return '-'; break;
    default:
      return '?'; break;
  }
}

BioType stringToType(const std::string &s) {
  if (s == "PROTEIN") {
    return BioType::PROTEIN;
  } else if (s == "PSEUDO") {
    return BioType::PSEUDO;
  } else {
    return BioType::OTHER;
  }
}

Strandedness charToStrand(char c) {
  switch (c) {
    case '+':
      return Strandedness::FORWARD; break;
    case '-':
      return Strandedness::REVERSE; break;
    default:
      return Strandedness::UNKNOWN; break;
  }
}

struct ExonModel {
  std::string chr;
  int start,stop;
  Strandedness strand;
};

struct TranscriptModel {
  std::string id; // like ENST
  std::string chr;
  int start,stop;
  std::vector<ExonModel> exons;
  Strandedness strand;
  BioType type;
};


struct GeneModel {
  std::string id; // like ENSG
  std::string name; // like ALB
  std::unordered_map<std::string,TranscriptModel> transcripts;
  std::string chr;
  int start,stop; // start 0-based, 1 bp past end for stop
  Strandedness strand;
  BioType type;
};

struct Transcriptome {
  std::unordered_map<std::string,GeneModel> genes;
  std::unordered_map<std::string,std::string> trxToGeneId;
  std::unordered_map<std::string, seqan::CharString> seqs;
  
  // maps transcript tr and 0-based position trpos
  //      to chr:chrpos in genome mapping to gene gene_id
  bool translateTrPosition(const std::string &tr, const int _trpos, std::string &chr, int& chrpos, std::string &gene_id) {
    auto gid_it = trxToGeneId.find(tr);
    if (gid_it == trxToGeneId.end()) {
      return false;
    }
    gene_id = gid_it->second;
    auto g_it = genes.find(gene_id);
    if (g_it == genes.end()) {
      return false;
    }

    auto t_it = g_it->second.transcripts.find(tr);
    if (t_it == g_it->second.transcripts.end()) {
      return false;
    }

    auto &trxmodel = t_it->second;
    chr = trxmodel.chr;    
    if (trxmodel.strand == Strandedness::UNKNOWN) {
      return false;
    } 
    bool fw = (trxmodel.strand == Strandedness::FORWARD);
    int trpos = _trpos;
    chrpos = trxmodel.start;
    for (auto& exon : trxmodel.exons) {
      int len = (exon.strand == Strandedness::FORWARD) ? (exon.stop - exon.start) : (exon.start - exon.stop);
      assert(len > 0);
      if (trpos < len) {
        // maps to this exon
        chrpos = exon.start;
        if (exon.strand == Strandedness::FORWARD) {
          chrpos += trpos;
        } else {
          chrpos -= trpos;
        }
        break;
      } else {
        trpos -= len;
      }
    }
    if (trpos > 0) {
      // goes beyond last exon, map to end
      chrpos = trxmodel.stop + trpos - 1;
    }
    return true;
  }
};

void writeTranscriptome(Transcriptome &transcriptome, std::ostream &out) {
  for(const auto &gene : transcriptome.genes) {
    const auto &gene_id = gene.first;
    const auto &model = gene.second;
    
    out << "GENE" << "\t" << model.id  << "\t"
        << model.name << "\t"
        << model.chr << "\t"
        << typeToString(model.type) << "\t"
        << strandToChar(model.strand) << "\t"
        << model.start << "\t"
        << model.stop << "\t";
    
    bool firsttr = true;
    for (auto &trlist : model.transcripts) {
      if (!firsttr) {
        out << ";";
      } else {
        firsttr = false;
      }
      out << trlist.first;
    }
    out << "\n";
  }
  for(const auto &gene : transcriptome.genes) {
    const auto &gene_id = gene.first;
    const auto &model = gene.second;
    for (auto &trlist : model.transcripts) {
      const auto &tr_id = trlist.first;
      const auto &tr = trlist.second;
      
      out << "TRANSCRIPT" << "\t" 
          << tr.id << "\t"
          << gene_id << "\t"
          << tr.chr << "\t"
          << typeToString(tr.type) << "\t"
          << strandToChar(tr.strand) << "\t"
          << tr.start << "\t"
          << tr.stop << "\t";
      
      bool firstex = true;
      for(const auto &exon : tr.exons) {
        if (!firstex) {
          out << ";";
        } else {
          firstex = false;
        }
        out << exon.start << "," << exon.stop;
      }
      out << "\n";
    }
  }
}

void loadTranscriptome(Transcriptome& transcriptome, std::istream &in) {
  std::string line;
  std::string segment;
  std::string type;
  std::string gtype, ttype;
  std::string gene_id, tr_id;
  while (in.good()) {
    in >> type;
    if (type == "GENE") {
      GeneModel model;
      char strand = '?';
      in >> gene_id>> model.name >> model.chr >> gtype >> strand >> model.start >> model.stop;
      model.id = gene_id;
      model.strand = charToStrand(strand);
      model.type = stringToType(gtype);
      in >> line; // rest of transcripts, ignore for now
      transcriptome.genes.insert({gene_id, std::move(model)});
    } else if (type == "TRANSCRIPT") {
      TranscriptModel model;
      char strand = '?';
      in >> tr_id >> gene_id >> model.chr >> ttype >> strand >> model.start >> model.stop;
      model.id = tr_id;
      model.strand = charToStrand(strand);
      model.type = stringToType(ttype);
      in >> line;
      std::stringstream strline(line);
      while (std::getline(strline, segment, ';')) {
        ExonModel ex;
        ex.chr = model.chr;
        ex.strand = model.strand;
        int p = segment.find(',');
        ex.start = std::stoi(segment.substr(0,p));
        ex.stop = std::stoi(segment.substr(p+1));
        model.exons.push_back(std::move(ex));
      }
      auto it = transcriptome.genes.find(gene_id);
      if (it != transcriptome.genes.end()) {
        GeneModel &gene_model = it->second;
        gene_model.transcripts.insert({tr_id, std::move(model)});
        transcriptome.trxToGeneId.insert({tr_id, gene_id});
      } else {
        std::cerr << "unknown gene id " << gene_id << std::endl;
      }
    } else {
      std::getline(in, line); // read until end of line
    }
  }
}


void parseFasta(Transcriptome &transcriptome, const std::string &fasta_fn) {
  seqan::SeqFileIn seqFileIn(seqan::toCString(fasta_fn));
  seqan::CharString id;
  seqan::CharString seq;

  while(!atEnd(seqFileIn)) {
    seqan::readRecord(id,seq,seqFileIn);
    std::string name = std::string(seqan::toCString(id));
    size_t sp = name.find(' ');
    size_t dot = name.find('.');
    transcriptome.seqs.insert({name.substr(0,std::min(sp,dot)), std::move(seq)});
  }  
}

void parseGTF(Transcriptome &transcriptome, const std::string &gtf_fn) {
  seqan::GffFileIn gtf(gtf_fn.c_str());
  seqan::GffRecord record;

  int n = 0;
  while(!seqan::atEnd(gtf)) {
    n++;
    seqan::readRecord(record,gtf);
    if (record.type == "gene") {
      GeneModel model;
      for (int i = 0; i < length(record.tagNames); i++) {
        if (record.tagNames[i] == "gene_id") {
          model.id = seqan::toCString(record.tagValues[i]);
        }
        if (record.tagNames[i] == "gene_name") {
          model.name = seqan::toCString(record.tagValues[i]);
        }
        if (record.tagNames[i] == "gene_biotype") {
          std::string val = std::string(seqan::toCString(record.tagValues[i]));
          //auto &val = record.tagValues[i];
          if (val == "protein_coding") {
            model.type = BioType::PROTEIN;
          } else if (val.find("pseudogene") != std::string::npos) {
            model.type = BioType::PSEUDO;
          } else {
            model.type = BioType::OTHER;
          }
        }

      }
      model.chr = seqan::toCString(record.ref);
      model.start = record.beginPos;
      model.stop = record.endPos;
      if (record.strand == '+') {
        model.strand = Strandedness::FORWARD;
      } else if (record.strand == '-') {
        model.strand = Strandedness::REVERSE;
      } else {
        model.strand = Strandedness::UNKNOWN;
      }
      transcriptome.genes.insert({model.id,std::move(model)});
    } else if (record.type == "transcript") {
      TranscriptModel model;
      std::string gene_id;
      for (int i = 0; i < length(record.tagNames); i++) {
        if (record.tagNames[i] == "gene_id") {
          gene_id = seqan::toCString(record.tagValues[i]);
        } else if (record.tagNames[i] == "transcript_id") {
          model.id = seqan::toCString(record.tagValues[i]);
        }
        if (record.tagNames[i] == "transcript_biotype") {
          std::string val = std::string(seqan::toCString(record.tagValues[i]));
          //auto &val = record.tagValues[i];
          if (val == "protein_coding") {
            model.type = BioType::PROTEIN;
          } else if (val.find("pseudogene") != std::string::npos) {
            model.type = BioType::PSEUDO;
          } else {
            model.type = BioType::OTHER;
          }
        }
      }
      model.chr = seqan::toCString(record.ref);
      model.start = record.beginPos;
      model.stop = record.endPos;
      if (record.strand == '+') {
        model.strand = Strandedness::FORWARD;
      } else if (record.strand == '-') {
        model.strand = Strandedness::REVERSE;
      } else {
        model.strand = Strandedness::UNKNOWN;
      }
      assert(!gene_id.empty());
      auto it = transcriptome.genes.find(gene_id);
      assert(transcriptome.trxToGeneId.find(model.id) == transcriptome.trxToGeneId.end());
      transcriptome.trxToGeneId.insert({model.id,gene_id});
      assert(it != transcriptome.genes.end());
      it->second.transcripts.insert({model.id,model});      
    } else if (record.type == "exon") {
      ExonModel model;
      std::string trx_id;
      std::string gene_id;
      for (int i = 0; i < length(record.tagNames); i++) {
         if (record.tagNames[i] == "gene_id") {
          gene_id = seqan::toCString(record.tagValues[i]);
        } else if (record.tagNames[i] == "transcript_id") {
          trx_id = seqan::toCString(record.tagValues[i]);
        }
      }
      model.chr = seqan::toCString(record.ref);
      model.start = record.beginPos;
      model.stop = record.endPos;
      if (record.strand == '+') {
        model.strand = Strandedness::FORWARD;
      } else if (record.strand == '-') {
        model.strand = Strandedness::REVERSE;
      } else {
        model.strand = Strandedness::UNKNOWN;
      }

      auto g_it = transcriptome.genes.find(gene_id);
      assert(g_it != transcriptome.genes.end());
      auto t_it = g_it->second.transcripts.find(trx_id);
      assert(t_it != g_it->second.transcripts.end());
      t_it->second.exons.push_back(std::move(model));
    }
  }
  std::cerr << "Found " << transcriptome.genes.size() << " genes out of " << n << " records" << std::endl;

  for (auto it : transcriptome.genes) {
    const auto &gene = it.second;
    //std::cout << "Gene = " << gene.id << ", name = " << gene.name << ", pos = " << gene.chr << ":" << (gene.start+1) << "-" << gene.stop << std::endl;
    
  }
}






#endif // GENE_MODEL_HPP