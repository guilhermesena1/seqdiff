#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <unistd.h>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "htslib_wrapper.hpp"
#include "chromosome_utils.hpp"
#include "cigar_utils.hpp"

using std::runtime_error;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::begin;
using std::end;
using std::unordered_map;
using std::unordered_set;
using std::to_string;

static bool
file_exists(const string &file) {
  return access(file.c_str(), R_OK) == 0;
}

struct pos_stats {
  uint16_t cnt_A;
  uint16_t cnt_C;
  uint16_t cnt_G;
  uint16_t cnt_T;
  uint16_t cnt_DEL;

  uint16_t pos_bases;
  uint16_t neg_bases;
  static const char DEL_BASE = '#';
  void reset() { 
    cnt_A = cnt_C = cnt_G = cnt_T = cnt_DEL = pos_bases = neg_bases = 0;
  }

  pos_stats() {
    reset();
  }

  void add_base(const char c, const bool is_neg) {
    if (c == 'A') ++cnt_A;
    else if (c == 'C') ++cnt_C;
    else if (c == 'G') ++cnt_G;
    else if (c == 'T') ++cnt_T;
    else if (c == DEL_BASE) ++cnt_DEL;
    else return;
    pos_bases += !is_neg;
  }

  string tostring() const {
    const size_t _depth = depth();
    const double denom = static_cast<double>(_depth);
    return to_string(_depth) + '\t' +
           to_string(pos_bases) + '\t' +
           to_string(_depth - pos_bases) + '\t' +
           to_string(cnt_A/denom) + '\t' +
           to_string(cnt_C/denom) + '\t' +
           to_string(cnt_G/denom) + '\t' +
           to_string(cnt_T/denom) + '\t' +
           to_string(cnt_DEL/denom);
  }

  uint16_t depth() const {
    return cnt_A + cnt_C + cnt_G + cnt_T + cnt_DEL;
  }
};

static void
get_chrom(const string &chrom_name, const vector<string> &all_chroms,
          const unordered_map<string, size_t> &chrom_lookup,
          string &chrom, vector<pos_stats> &stats) {
  const auto the_chrom = chrom_lookup.find(chrom_name);
  if (the_chrom == end(chrom_lookup))
    throw runtime_error("could not find chrom " + chrom_name);

  chrom = all_chroms[the_chrom->second];

  const size_t chrom_sz = chrom.size();
  stats.resize(chrom_sz);
  for (size_t i = 0; i < chrom_sz; ++i)
    stats[i].reset();

  if (chrom.empty())
    throw runtime_error("problem with chrom sequence " + chrom_name);
}

static void
load_chroms(const bool VERBOSE,
            const string &genome_file,
            vector<string> &all_chroms,
            vector<string> &chrom_names,
            unordered_map<string, size_t> &chrom_lookup) {
  if (VERBOSE)
    cerr << "[loading FASTA file: " << genome_file << "]\n";

  vector<string> tmp_chroms;
  vector<string> tmp_names;
  read_fasta_file(genome_file, tmp_names, tmp_chroms);
  for (size_t i = 0; i < tmp_chroms.size(); ++i) {
    chrom_names.push_back(tmp_names[i]);
    chrom_lookup[chrom_names.back()] = all_chroms.size();
    all_chroms.push_back("");
    all_chroms.back().swap(tmp_chroms[i]);
    if (VERBOSE)
      cerr << "  chrom: " << chrom_names.back() << "\tsize: "
           << all_chroms.back().size() << "\n";
  }

  if (VERBOSE)
    cerr << "[read " << all_chroms.size() << " chromosomes]\n";
}

inline void
swap_bisulfite_base(sam_rec &aln) {
  // GS: pass as command args?
  static const size_t the_cv_tag_pos = 1;
  const char base = aln.tags[the_cv_tag_pos].back();
  if (base == 'A')
    aln.tags[the_cv_tag_pos] = 'T';
  else if (base == 'T')
    aln.tags[the_cv_tag_pos] = 'A';
  else
    throw runtime_error("bad bisulfite base: " + base);
}

template<const bool bisulfite_bases>
void
adjust_alignment(sam_rec &aln, bool &is_neg) {
  // (1) make pos zero-based
  --aln.pos;

  // (2) revcomp to match the available sequence
  if (check_flag(aln, samflags::read_rc)) {
    is_neg = true;
    unset_flag(aln, samflags::read_rc);
    revcomp_inplace(aln.seq);
    if (bisulfite_bases)
      swap_bisulfite_base(aln);
  }

  // (3) apply CIGAR with deletion symbol
  inflate_with_cigar(aln, aln.seq, pos_stats::DEL_BASE);
}

inline char
get_bisulfite_base(const sam_rec &aln) {
  static const size_t the_cv_tag_pos = 1;
  return aln.tags[the_cv_tag_pos].back();
}

template<const bool bisulfite_bases>
void
process_alignment(const bool is_neg, const sam_rec &aln,
                  vector<pos_stats> &stats) {
  const size_t lim = aln.seq.size();
  const size_t start = aln.pos;
  const char bisulfite_base = get_bisulfite_base(aln);
  for (size_t i = 0; i < lim; ++i) {
    if (!bisulfite_bases || aln.seq[i] != bisulfite_base) {
      stats[start + i].add_base(aln.seq[i], is_neg);
    }
  }
}

static void
summarize_chrom_stats(const string & cur_chrom, const vector<pos_stats> &stats) {
  const size_t lim = cur_chrom.size();
  cout << "pos\tbase\tdepth\t+\t-\tA\tC\tG\tT\tDEL\n";
  for (size_t i = 0; i < lim; ++i) {
    if (stats[i].depth() > 0)
      cout << i << '\t' << cur_chrom[i] << '\t' << stats[i].tostring() << '\n';
  }
}

template<const bool bisulfite_bases>
static void
process_reads(const bool VERBOSE, const string &mapped_reads_file,
              const vector<string> &all_chroms,
              const vector<string> &chrom_names,
              const unordered_map<string, size_t> &chrom_lookup) {
  SAMReader sam_reader(mapped_reads_file);
  sam_rec aln;

  unordered_set<string> chroms_seen;
  string cur_chrom, cur_chrom_name = "";

  vector<pos_stats> stats;
  size_t num_reads_for_chrom = 0;
  while (sam_reader >> aln) {
    if (aln.rname != cur_chrom_name) {
      // new chrom to process
      if (chroms_seen.find(aln.rname) != end(chroms_seen))
        throw runtime_error("chroms out of order in mapped reads file\n");

      if (num_reads_for_chrom > 0) {
        if (VERBOSE)
          cerr << "[writing stats for " << cur_chrom_name << "]\n";

        cout << ">" << cur_chrom_name <<"\n";
        summarize_chrom_stats(cur_chrom, stats);
      }
      if (VERBOSE)
        cerr << "[processing chromosome " << aln.rname << "]\n";

      cur_chrom_name = aln.rname;
      chroms_seen.insert(cur_chrom_name);
      get_chrom(cur_chrom_name, all_chroms, chrom_lookup, cur_chrom, stats);
    }
    ++num_reads_for_chrom;
    bool is_neg = false;
    adjust_alignment<bisulfite_bases>(aln, is_neg);
    process_alignment<bisulfite_bases>(is_neg, aln, stats);
  }

  if (num_reads_for_chrom > 0) {
    if (VERBOSE)
      cerr << "[writing stats for " << cur_chrom_name << "]\n";

    cout << ">" << cur_chrom_name <<"\n";
    summarize_chrom_stats(cur_chrom, stats);
  }
}

int
main(const int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool bisulfite = false;
    OptionParser opt_parse(strip_path(argv[0]),
                           "estimate edit distance between reference and mapped reads",
                           "<reference.fa> <seq-file.sam>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("bisulfite", 'b', "reads come from WGBS", false, bisulfite);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }

    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }

    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }

    cout.precision(3);

    // program starts here
    const string genome_file = leftover_args.front();
    const string mapped_reads_file = leftover_args.back();

    if (!file_exists(genome_file))
      throw runtime_error("cannot open genome file: " + genome_file);

    if (!file_exists(mapped_reads_file))
      throw runtime_error("cannot open mapped reads file: " +
          mapped_reads_file);

    vector<string> all_chroms;
    vector<string> chrom_names;
    unordered_map<string, size_t> chrom_lookup;
    load_chroms(VERBOSE, genome_file, all_chroms, chrom_names, chrom_lookup);
    if (bisulfite)
      process_reads<true>(VERBOSE, mapped_reads_file,
                          all_chroms, chrom_names, chrom_lookup);
    else
      process_reads<false>(VERBOSE, mapped_reads_file,
                          all_chroms, chrom_names, chrom_lookup);

  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
