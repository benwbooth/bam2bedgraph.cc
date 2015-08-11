#include <cstdio>
#include <cstdlib>
#include <sys/wait.h>
#include <algorithm>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <tuple>
#include <cmath>

#include <api/BamReader.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include "IntervalTree.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::pair;
using std::unique_ptr;
using std::map;
using std::tuple;
using std::to_string;
using std::ifstream;
using std::ofstream;
using std::abs;
using BamTools::BamReader;
using BamTools::SamHeader;
using BamTools::RefData;
using BamTools::BamAlignment;
using BamTools::CigarOp;
using boost::filesystem::path;
using boost::algorithm::to_lower;
using boost::format;
using boost::numeric_cast;
using boost::regex;


static bool split_exons=true, split_read, zero, paired_only, proper_only, primary_only, trackline=true, bigwig, uniq, fixchr;
static string bamfile, trackname, out, autostrand, split_strand="uu";

void cigar2exons(vector<pair<size_t,size_t>> &exons, const vector<CigarOp> &cigar, size_t pos) {
  for (auto &op : cigar) {
    if (op.Type == 'M') {
      pos += op.Length;
      exons.push_back({pos - op.Length, pos});
    }
    else if (op.Type == 'N' || op.Type == 'D') {
      pos += op.Length;
    }
    else if (op.Type == 'I' || op.Type == 'S' || op.Type == 'H' || op.Type == 'P') {}
    else {
      throw string("Bad CIGAR string: ")+op.Type;
    }
  }
}

string
open_file(int32_t read_number,
          string strand,
          string split_strand,
          map<string,unique_ptr<ofstream>> &fhs)
{
  string track_name = (format("%s%s%s") %
                       (!trackname.empty()? trackname : path{bamfile}.filename().replace_extension(path{}).native()) %
                       (split_read && read_number? (format(".r%s") % read_number).str() : "") %
                       (split_strand != "uu" && !strand.empty()? (format(".%s") % strand).str() : "")).str();

  string filename = (format("%s%s%s%s") %
                     (!out.empty()? out : path{bamfile}.replace_extension(path{}).native()) %
                     (split_read && read_number? (format(".r%s") % read_number).str() : "") %
                     (split_strand != "uu" && !strand.empty()? (format(".%s") % strand).str() : "") %
                     ".bedgraph").str();

  // initialize the file if needed
  if (!fhs.count(filename)) {
    fhs[filename] = unique_ptr<ofstream>{new ofstream{}};
    auto &fh = *fhs[filename];
    fh.exceptions(std::ios::badbit | std::ios::failbit);
    fh.open(filename);
    if (trackline) fh <<"track type=bedGraph name=\""<<track_name<<"\" description=\""<<track_name<<"\" visibility=full"<<endl;
  }
  return filename;
}

void write_chr(
    const RefData &chr,
    const map<tuple<int32_t,string>,vector<int32_t>> &histogram,
    map<string,unique_ptr<ofstream>> &fhs,
    string split_strand)
{
  for (auto& tuple : histogram) {
    auto read_number = std::get<0>(tuple.first);
    auto strand = std::get<1>(tuple.first);
    auto& histo = tuple.second;

    string filename = open_file(read_number, strand, split_strand, fhs);
    auto &fh = *fhs[filename];
    // scan the histogram to produce the bedgraph data
    size_t start = 0;
    size_t end = 0;
    size_t ref_length = numeric_cast<size_t>(chr.RefLength);
    while (start < ref_length) {
      while ((end < histo.size()? histo[end] : 0) == (start < histo.size()? histo[start] : 0) && end < ref_length) ++end;
      if (zero || (start < histo.size()? histo[start] : 0)) {
        fh << chr.RefName << "\t" <<
              start << "\t" <<
              end << "\t" <<
              (strand == "-"? -histo[start] : histo[start]) << endl;
      }
      start = end;
    }
  }
}

static void
analyzeBam(string split_strand,
           bool autostrandPass,
           map<string,IntervalTree<int8_t>>& intervals)
{
  BamReader bam;
  if (!bam.Open(bamfile)) {
    throw (format("Could not open input BAM file ") % bamfile).str();
  }
  auto header = bam.GetHeader();
  auto refs = bam.GetReferenceData();
  if (fixchr) {
    for (auto& ref : refs) {
      ref.RefName = !boost::regex_search(ref.RefName, boost::regex{R"(^(chr|Zv9_))"})?
        string("chr")+ref.RefName :
        ref.RefName;
    }
  }

  if (autostrandPass) {
    cerr << "Running strand detection phase on " << bamfile << endl;
  }
  else {
    cerr << "Building histograms for " << bamfile << endl;
  }

  // build a lookup map for the refseqs
  map<string,size_t> refmap;
  for (size_t i=0; i<refs.size(); ++i) {
    refmap[refs[i].RefName] = i;
  }

  int32_t lastchr = -1;
  map<string,unique_ptr<ofstream>> fhs;
  map<tuple<int32_t,string>,vector<int32_t>> histogram;

  map<char,int64_t> autostrandTotals{};
  map<char,int64_t> autostrandTotals2{};
  BamAlignment read;
  while (bam.GetNextAlignment(read)) {
    // if we've hit a new chr, write out the bedgraph data and clear the histogram
    if (lastchr == -1 || read.RefID != lastchr) {
      if (!autostrandPass && !histogram.empty() && lastchr != -1) {
        write_chr(refs[lastchr], histogram, fhs, split_strand);
      }
      histogram.clear();
      lastchr = read.RefID;
    }

    // skip this read if it's no good
    bool paired = read.IsPaired();
    bool proper = read.IsProperPair();
    bool primary  = read.IsPrimaryAlignment();
    if ((paired_only && !paired) || (primary_only && !primary) || (proper_only && !proper)) continue;

    // skip if it's not unique and we want unique alignments
    if (uniq) {
      int32_t hits=0;
      if (!read.GetTag("NH",hits) || hits != 1) continue;
    }

    vector<pair<size_t,size_t>> exons;
    if (split_exons) cigar2exons(exons, read.CigarData, read.Position);
    else exons.push_back({read.Position, read.GetEndPosition()});
    int32_t read_number = read.IsSecondMate()? 2 : 1;

    // attempt to determine the strandedness of the transcript
    uint8_t xs = 0;
    // read numbers match, is not reverse, is not flipped
    string strand =
        read.GetTag("XS",xs) && xs? string{static_cast<char>(xs)} :
          read_number == 1? split_strand[0] == 'r'? read.IsReverseStrand()? "+" : "-" :
                            split_strand[0] == 's'? read.IsReverseStrand()? "-" : "+" :
                            "" :
          read_number == 2? split_strand[1] == 's'? read.IsReverseStrand()? "-" : "+" :
                            split_strand[1] == 'r'? read.IsReverseStrand()? "+" : "-" :
                            "" :
        "";

    int32_t read_num = split_read? read_number : 0;

    size_t ref_length = numeric_cast<size_t>(refs[read.RefID].RefLength);
    // add the read to the histogram
    for (auto &exon : exons) {
      // try to determine the strandedness of the data
      if (autostrandPass) {
        if (intervals.count(refs[lastchr].RefName)) {
          vector<Interval<int8_t>> overlappingAnnot;
          intervals[refs[lastchr].RefName].findOverlapping(exon.first+1, exon.second, overlappingAnnot);
          for (const auto& interval : overlappingAnnot) {
            size_t overlap_length = std::min(exon.second, numeric_cast<size_t>(interval.stop)) -
              std::max(exon.first, numeric_cast<size_t>(interval.start-1));

            char strandtype = read.IsReverseStrand() == (interval.value == '-')? 's' : 'r';
            if (read_number == 1) autostrandTotals[strandtype] += overlap_length;
            else if (read_number == 2) autostrandTotals2[strandtype] += overlap_length;
          }
        }
      }
      else {
        auto tuple = std::make_tuple(read_num, strand);
        if (!histogram.count(tuple)) histogram.insert({tuple, vector<int32_t>{}});
        // keep track of chromosome sizes
        if (ref_length < exon.second) refs[read.RefID].RefLength = numeric_cast<int32_t>(exon.second);
        if (histogram[tuple].size() < ref_length) histogram[tuple].resize(ref_length);

        for (size_t pos=exon.first; pos < exon.second; ++pos) {
          histogram[tuple][pos]++;
        }
      }
    }
  }
  bam.Close();

  if (!autostrandPass && !histogram.empty() && lastchr != -1) {
    write_chr(refs[lastchr], histogram, fhs, split_strand);
  }

  // make sure empty files were created
  if (histogram.empty() && !autostrandPass) {
    for (auto &read_number : split_read? vector<int32_t>{1,2} : vector<int32_t>{0}) {
      for (auto &s : split_strand != "uu"? vector<string>{"+","-"} : vector<string>{""}) {
        open_file(read_number, s, split_strand, fhs);
      }
    }
  }

  // close the filehandles
  for (auto &fh : fhs) {
    fh.second->close();
  }
  if (autostrandPass) {
    // get the read 1 and read2 totals
    int64_t total1 = 0;
    int64_t total2 = 0;
    for (auto &i : autostrandTotals) {
      total1 += i.second;
    }
    for (auto &i : autostrandTotals2) {
      total2 += i.second;
    }
    // figure out the best and second-best strand types for reads 1 and 2
    pair<int32_t,int64_t> best1;
    pair<int32_t,int64_t> second_best1;
    pair<int32_t,int64_t> best2;
    pair<int32_t,int64_t> second_best2;
    for (auto &i : autostrandTotals) {
      cerr << "Total evidence for read 1 strand type " << i.first << ": " << i.second << endl;
      if (best1.second < i.second) {
        second_best1 = best1;
        best1 = i;
      }
      else if (second_best1.second < i.second) {
        second_best1 = i; 
      }
    }
    for (auto &i : autostrandTotals2) {
      cerr << "Total evidence for read 2 strand type " << i.first << ": " << i.second << endl;
      if (best2.second < i.second) {
        second_best2 = best2;
        best2 = i;
      }
      else if (second_best2.second < i.second) {
        second_best2 = i;
      }
    }
    const double threshold = 0.0; //threshold isn't working, set to zero
    char strand1 = (total1 > 0.0)?
      threshold < numeric_cast<double>(best1.second - second_best1.second) / numeric_cast<double>(total1)? best1.first : 'u'
                  : 'u';
    char strand2 = (total2 > 0.0)?
      threshold < numeric_cast<double>(best2.second - second_best2.second) / numeric_cast<double>(total2)? best2.first : 'u'
                  : 'u';
    string best_strand {strand1, strand2};
    cerr << "autostrandPass found best strand type: " << best_strand << endl;

    // re-run analyzeBam with the strand type indicated
    analyzeBam(best_strand, false, intervals);
  }
  if (!autostrandPass && bigwig) {
    // Convert bedgraph file to bigwig file
    for (auto &fh : fhs) {
      // write the genome file for bigwigs
      string genome_filename = fh.first+".genome";
      ofstream genome_fh {};
      genome_fh.exceptions(std::ios::badbit | std::ios::failbit);
      genome_fh.open(genome_filename);
      for (auto& ref : refs) {
        genome_fh << ref.RefName << "\t" << ref.RefLength << endl;
      }
      genome_fh.close();

      string cmd = (format("bedGraphToBigWig '%s' '%s' '%s'") %
          boost::regex_replace(fh.first, regex(R"(')"), "'\\''") %
          boost::regex_replace(genome_filename, regex(R"(')"), "'\\''") %
          boost::regex_replace(boost::regex_replace(fh.first, regex(R"(\.bedgraph$)"),"")+".bw", regex(R"(')"), "'\\''")).str();
      int result = system(cmd.c_str());
      if (WEXITSTATUS(result) != 0) {
        throw (format("Command \"%s\" returned bad exit status: %d") % cmd % WEXITSTATUS(result)).str();
      }

      // remove the bedgraph file
      boost::filesystem::remove(fh.first);
      // remove the genome file
      boost::filesystem::remove(genome_filename);
    }
  }
}

int main(int argc, char **argv) {
  try {
    using boost::program_options::options_description;
    using boost::program_options::value;
    using boost::program_options::bool_switch;
    using boost::program_options::positional_options_description;
    using boost::program_options::variables_map;
    using boost::program_options::command_line_parser;
    namespace cls = boost::program_options::command_line_style;
    using boost::program_options::notify;

    bool nosplit=false, noread=false, nozero=false, nofixchr=false,
      nopaired=false, noproper=false, noprimary=false,nobigwig=false,
      nouniq=false, notrackline=false;

    options_description desc;
    desc.add_options()
      ("help", "produce help")
      ("bamfile", value<string>(&bamfile)->value_name("FILE"),
       "Input BAM filename")
      ("split", bool_switch(&split_exons)->default_value(true),
       "Use CIGAR string to split alignment into separate exons (default)")
      ("nosplit", bool_switch(&nosplit))
      // if your protocol is from Illumina (fr-secondstrand), read_one is on the correct strand
      // if your protocol is dUTP (ala CSH) (fr-firststrand), then read_two is on the correct strand
      ("autostrand", value<string>(&autostrand)->value_name("ANNOT_FILE"),
       "Attempt to determine the strandedness of the input data using an annotation file. Must be a BAM file.")
      ("strand", value<string>(&split_strand)->value_name("[TYPE]"),
       "Split output bedgraph by strand: Possible values: u s r uu us ur su ss sr ru rs rr, first char is read1, second is read2, u=unstranded, s=stranded, r=reverse")
      ("read", bool_switch(&split_read), "Split output bedgraph by read number")
      ("noread", bool_switch(&noread), "(default)")
      ("zero", bool_switch(&zero), "Pad output bedgraph with zeroes")
      ("nozero", bool_switch(&nozero), "(default)")
      ("fixchr", bool_switch(&fixchr), "Transform chromosome names to be UCSC-compatible")
      ("nofixchr", bool_switch(&nofixchr), "(default)")
      ("paired", bool_switch(&paired_only), "Only output paired read alignments")
      ("nopaired", bool_switch(&nopaired), "(default)")
      ("proper", bool_switch(&proper_only), "Only output proper-paired read alignments")
      ("noproper", bool_switch(&noproper), "(default)")
      ("primary", bool_switch(&primary_only), "Only output primary alignments")
      ("noprimary", bool_switch(&noprimary), "(default)")
      ("bigwig", bool_switch(&bigwig), "Output bigwig files (requires bedGraphToBigWig in $PATH)")
      ("nobigwig", bool_switch(&nobigwig), "(default)")
      ("uniq", bool_switch(&uniq), "Keep only unique alignments (NH:i:1)")
      ("nouniq", bool_switch(&nouniq), "(default)")
      ("out", value<string>(&out)->value_name("FILE"), "Output file prefix")
      ("trackline", bool_switch(&trackline)->default_value(true), "Output a UCSC track line (default)")
      ("notrackline", bool_switch(&notrackline))
      ("trackname", value<string>(&trackname)->value_name("TRACKNAME"), "Name of track for the track line") ;
    positional_options_description pod;
    pod.add("bamfile", 1);
    variables_map vm;
    store(command_line_parser(argc, argv).options(desc).
      style(cls::default_style | cls::allow_long_disguise).
      positional(pod).run(), vm);

    if (vm.count("help") || !vm.count("bamfile")) {
      if (!vm.count("bamfile")) cerr << "No input bam file was given" << endl;
      cerr << "Usage: " << argv[0] << " [ options ] bamfile" << endl << desc << endl;
      return 1;
    }
    notify(vm);

    if (nosplit) split_exons = false;
    if (noread) split_read = false;
    if (nozero) zero = false;
    if (nofixchr) fixchr = false;
    if (nopaired) paired_only = false;
    if (noproper) proper_only = false;
    if (noprimary) primary_only = false;
    if (nobigwig) bigwig = false;
    if (nouniq) uniq = false;
    if (notrackline) trackline = false;

    boost::algorithm::to_lower(split_strand);
    if (split_strand.size() == 1) split_strand += "u";
    if (!boost::regex_search(split_strand, boost::regex{R"(^[usr][usr]$)"})) {
      cerr << "Invalid value for split_strand: \"" << split_strand << "\": values must be one of: u s r uu us ur su ss sr ru rs rr" << endl;
      return 1;
    }

    // read in the annotation file
    map<string,vector<Interval<int8_t>>> interval_lists;
    BamReader bam;
    if (!autostrand.empty()) {
      if (!bam.Open(autostrand)) {
        throw (format("Could not open strand annotation BAM file %s") % autostrand).str();
      }
      auto refs = bam.GetReferenceData();
      if (fixchr) {
        for (auto& ref : refs) {
          ref.RefName = !boost::regex_search(ref.RefName, boost::regex{R"(^(chr|Zv9_))"})?
            string("chr")+ref.RefName :
            ref.RefName;
        }
      }
      BamAlignment read;
      while (bam.GetNextAlignment(read)) {
        int32_t tid = read.RefID;
        if (!interval_lists.count(refs[tid].RefName)) {
          interval_lists[refs[tid].RefName] = vector<Interval<int8_t>>{};
        }
        interval_lists[refs[tid].RefName].push_back(Interval<int8_t>{
            read.Position+1,
            read.GetEndPosition(),
            read.IsReverseStrand()? '-' : '+'});
      }
    }
    map<string,IntervalTree<int8_t>> intervals;
    for (auto &kv : interval_lists) {
      intervals[kv.first] = IntervalTree<int8_t>{kv.second};
    }

    // analyze the bam file and produce histograms
    if (!autostrand.empty()) {
      // make both stranded and unstranded files
      analyzeBam(split_strand, !autostrand.empty(), intervals);
      analyzeBam(string("uu"), false, intervals);
    }
    else {
      analyzeBam(split_strand, !autostrand.empty(), intervals);
    }
  }
  catch(const std::exception& e) {
    cerr << "Caught exception: " << e.what() << endl;
    return 1;
  }
  catch(const string& e) {
    cerr << "Caught exception: " << e << endl;
    return 1;
  }
  catch(...) {
    cerr << "Caught unknown exception" << endl;
    return 1;
  }
  return 0;
}
