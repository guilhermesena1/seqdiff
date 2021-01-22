#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <parallel/algorithm>
#include <unistd.h>
#include <cmath>
#include <stdexcept>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "htslib_wrapper.hpp"
#include "chromosome_utils.hpp"
#include "cigar_utils.hpp"

#include <omp.h>

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
using std::min;
using std::max;
using std::ostream;
using std::ostringstream;

typedef uint16_t count_t;

static bool
file_exists(const string &file) {
  return access(file.c_str(), R_OK) == 0;
}

//sizeof(pos_stats) = 64
struct pos_stats {
  count_t bases_ref;
  count_t bases_nref;
  count_t bases_pos;
  count_t bases_neg;
  
  void reset() {
    bases_ref = bases_nref = bases_pos = bases_neg = 0;
  }

  pos_stats() { reset(); }

  double get_geno_likelihood(const size_t geno_count) const;
  void add_base(const char seq_base, const char ref_base, const bool is_pos);
  string tostring() const;

  bool sure_homozygous() const {
    return ((depth() >= 5) && (10*max(bases_ref, bases_nref) >= 9*depth()));
  }

  bool sure_ref() const {
    return ((depth() >= 5) && (10*bases_ref >= 9*depth()));
  }

  inline count_t depth() const {
    return bases_ref + bases_nref;
  }

  double calc_psi(const double err) const;
  double loglik(const size_t g, const double err) const;
  static char DEL_BASE;
  static size_t num_bases;
  static size_t num_alleles;
  static double psi_not_calculated;
};

size_t pos_stats::num_bases = 5;
size_t pos_stats::num_alleles = 2;
char pos_stats::DEL_BASE = '#';
double pos_stats::psi_not_calculated = 2.0;

double
pos_stats::get_geno_likelihood(const size_t geno_count) const {
  return 1.0;
}

void
pos_stats::add_base(const char seq_base, const char ref_base,
                    const bool is_pos) {
  if (seq_base == 'N') return;
  const bool is_ref = (seq_base == ref_base);
  //assert(bases_ref + bases_nref == bases_pos + bases_neg);
  bases_ref += is_ref; bases_nref += !is_ref;
  bases_pos += is_pos; bases_neg += !is_pos;
}

double
pos_stats::loglik(const size_t g, const double err) const {
  //assert(g >= 0 && g <= num_alleles)
  static const double m = static_cast<double>(num_alleles);
  const double gg = static_cast<double>(g);
  return (static_cast<double>(bases_ref)*
           log((m - gg)*err + gg*(1 - err)) +
         static_cast<double>(bases_nref)*
           log((m - gg)*(1 - err) + gg*err) -
         static_cast<double>(depth())*log(2.0));
}

double
pos_stats::calc_psi(const double err) const {
  // optimization: no nref bases -> psi = 1
  if (bases_nref == 0) return 1.0;
  double psi = 0.5;
  double psi_prev = 0.0;
  double l0 = exp(loglik(0, err));
  double l1 = exp(loglik(1, err));
  double l2 = exp(loglik(2, err));

  if (isnan(l0) || isnan(l1) || isnan(l2))
    throw runtime_error("bad loglik: " + to_string(l0) + " " + to_string(l1) +
                        " " + to_string(l2));
  // assert(eps < fabs(psi - psi_prev))
  const double eps = 1e-10;
  size_t niter = 0;
  while (fabs(psi - psi_prev) > eps) {
    psi_prev = psi;

    const double num = 2*l1*psi*(1 - psi) + 2*l2*psi*psi;
    const double denom = l0*(1 - psi)*(1 - psi) + 2*l1*psi*(1 - psi) + l2*psi*psi;

    if (denom == 0.0) {
      // GS: too much data. fall back to naive conversion?
      // this is a hack. Need to implement log sum log exp
      return bases_ref / static_cast<double>(depth());
    }

    psi = num/(num_alleles*denom);
    if (++niter > 10000) {
      cerr << "[WARNING: psi calculation did not converge. Value:  " << psi << "]\n";
      return psi;
    }
  }
  
  // GS: not sure if needed, but it makes sense to truncate to zero when
  // psi is very small because it went under epsilon
  if (psi < eps) psi = 0;
  return psi;
}

string
pos_stats::tostring() const {
  ostringstream oss;
  oss << bases_ref << '\t' << bases_nref << '\t'
      << bases_pos << '\t' << bases_neg;

  return oss.str();
}

ostream &
operator<<(ostream &the_stream, const pos_stats &stats) {
  the_stream << stats.tostring();
  return the_stream;
}

struct total_counts {
  size_t total_a;
  size_t total_c;
  size_t total_g;
  size_t total_t;
  size_t total_del;

  // expected from the genome
  size_t exp_a;
  size_t exp_c;
  size_t exp_g;
  size_t exp_t;

  // wrong bases
  size_t err_a;
  size_t err_c;
  size_t err_g;
  size_t err_t;
  void reset() {
    total_a = total_c = total_g = total_t = total_del = 0;
    exp_a = exp_c = exp_g = exp_t = 0;
    err_a = err_c = err_g = err_t = 0;
  }
  total_counts() {
    reset();
  }

  void add_exp(const char ref_base) {
    exp_a += (ref_base == 'A');
    exp_c += (ref_base == 'C');
    exp_g += (ref_base == 'G');
    exp_t += (ref_base == 'T');
  }

  void add_base(const char seq_base, const char ref_base) {
    total_a += (seq_base == 'A');
    total_c += (seq_base == 'C');
    total_g += (seq_base == 'G');
    total_t += (seq_base == 'T');
    total_del += (seq_base == pos_stats::DEL_BASE);

    err_a += (seq_base == 'A' && ref_base != 'A');
    err_c += (seq_base == 'C' && ref_base != 'C');
    err_g += (seq_base == 'G' && ref_base != 'G');
    err_t += (seq_base == 'T' && ref_base != 'T');
  }

  string tostring() const {
    static const string line_start = "    ";
    ostringstream oss;
    oss << line_start << "A_obs_exp: " << total_a / static_cast<double>(exp_a)
        << '(' << err_a / static_cast<double>(exp_a) << ")\n";
    oss << line_start << "C_obs_exp: " << total_c / static_cast<double>(exp_c)
        << '(' << err_c / static_cast<double>(exp_c) << ")\n";
    oss << line_start << "G_obs_exp: " << total_g / static_cast<double>(exp_g)
        << '(' << err_g / static_cast<double>(exp_g) << ")\n";
    oss << line_start << "T_obs_exp: " << total_t / static_cast<double>(exp_t)
        << '(' << err_t / static_cast<double>(exp_t) << ")\n";
    oss << line_start << pos_stats::DEL_BASE << ": " << total_del;
    return oss.str();
  }
};

ostream &
operator<<(ostream &the_stream, const total_counts &cnts) {
  the_stream << cnts.tostring();
  return the_stream;
}

static void
get_chrom(const string &chrom_name, const vector<string> &all_chroms,
          const unordered_map<string, size_t> &chrom_lookup,
          string &chrom, vector<pos_stats> &stats,
          total_counts &cnts) {
  const auto the_chrom = chrom_lookup.find(chrom_name);
  if (the_chrom == end(chrom_lookup))
    throw runtime_error("could not find chrom " + chrom_name);

  chrom = all_chroms[the_chrom->second];
  transform(begin(chrom), end(chrom), begin(chrom), ::toupper);

  const size_t chrom_sz = chrom.size();
  stats.resize(chrom_sz);

  cnts.reset();
  for (size_t i = 0; i < chrom_sz; ++i) {
    cnts.add_exp(chrom[i]);
    stats[i].reset();
  }

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
adjust_alignment(sam_rec &aln) {
  --aln.pos; // 0-based pos
  if (check_flag(aln, samflags::read_rc)) {
    revcomp_inplace(aln.seq);
    if (bisulfite_bases) swap_bisulfite_base(aln);
  }

  inflate_with_cigar(aln, aln.seq, pos_stats::DEL_BASE);
  // GS: (not sure if needed) uppercase seq
}

inline char
get_bisulfite_base(const sam_rec &aln) {
  static const size_t the_cv_tag_pos = 1;
  const char c = aln.tags[the_cv_tag_pos].back();
  if (c != 'A' && c != 'T')
    throw runtime_error("bad bisulfite base: " + aln.tostring());
  return c;
}

inline bool
is_bisulfite_base(const char c, const char b) {
  return (b == 'T') ?
         (c == 'T' || c == 'C') :
         (c == 'A' || c == 'G');
}

template<const bool do_bisulfite_bases>
void
process_alignment(const sam_rec &aln,
                  const string &chrom, vector<pos_stats> &stats,
                  total_counts &cnts) {
  char bisulfite_base = 0;
  if (do_bisulfite_bases)
    bisulfite_base = get_bisulfite_base(aln);
  const size_t lim = aln.seq.size();
  const size_t start = aln.pos;
  for (size_t i = 0; i < lim; ++i) {
    const char the_aln_base = aln.seq[i];
    if (!do_bisulfite_bases || !is_bisulfite_base(the_aln_base, bisulfite_base)) {
      stats[start + i].add_base(
        the_aln_base, chrom[start + i], // ref base
        !check_flag(aln, samflags::read_rc) // +/- strand
      );
      cnts.add_base(the_aln_base, chrom[start+i]);
    }
  }
}

static double
calc_error_freq(const vector<pos_stats> &stats) {
  size_t err = 0;
  size_t tot = 0;
  const size_t lim = stats.size();
  for (size_t i = 0; i < lim; ++i) {
    if (stats[i].sure_homozygous()) {
      err += min(stats[i].bases_ref, stats[i].bases_nref);
      tot += stats[i].depth();
    }
  }

  return err / static_cast<double>(tot);
}

static void
summarize_chrom_stats(const bool VERBOSE,
                      const string &cur_chrom, const string &cur_chrom_name,
                      vector<double> &psi,
                      const vector<pos_stats> &stats,
                      const total_counts &cnts,
                      double &similarity) {
  const double error_freq = calc_error_freq(stats);
  const size_t lim = cur_chrom.size();

  if (VERBOSE)
    cerr << "[calculating reference allele frequencies]\n";

#pragma omp parallel for
  for (size_t i = 0; i < lim; ++i) {
    if (stats[i].depth() > 0)
      psi[i] = stats[i].calc_psi(error_freq);
  }

  if (VERBOSE)
    cerr << "[calculating similarity to the reference]\n";

  double sim = 0.0;
  double tot = 0.0;
  for (size_t i = 0; i < lim; ++i) {
    if (stats[i].depth() > 0) {
      sim += psi[i];
      tot += 1.0;
    }
  }
  similarity = sim/tot;
}

static void
write_psi_histogram(vector<double> &psi, ostream &out) {
  // GS: histogram up to 3 sig digits

  static const string line_start = "    ";
  const size_t mult = 1000;
  const size_t lim = psi.size();
  if (lim == 0) return;

  __gnu_parallel::sort(begin(psi), end(psi));
  size_t cur = static_cast<size_t>(mult*psi[0]);
  const size_t undef = static_cast<size_t>(pos_stats::psi_not_calculated*mult);
  size_t cnt = 1;

  for (size_t i = 1; i != lim; ++i) {
    const size_t v = static_cast<size_t>(mult*psi[i]);
    if (v != cur) {
      if (cur == undef)
        out << line_start << "not_calculated";
      else out << line_start <<  cur/static_cast<double>(mult);

      out << ": " << cnt << '\n';
      cur = v;
      cnt = 1;
    } else ++cnt;
  }
  if (cur == undef) out << line_start << "not_calculated";
  else out << line_start << cur/static_cast<double>(mult);
  out << ": " << cnt << '\n';
}


static void
write_chrom_stats(const string &chrom_name,
                  const size_t &chrom_size,
                  const double &similarity,
                  vector<double> &psi,
                  const total_counts &cnts, ostream &out) {
  static const string sep = "  ";
  out << chrom_name << ":\n";
  out << sep << "size: " << chrom_size << "\n";
  out << sep << "similarity: " << similarity << "\n";
  out << sep << "global_base_stats:\n";
  out << cnts << "\n";
  out << sep << "psi_histogram:\n";
  write_psi_histogram(psi, out);
}

template<const bool bisulfite_bases>
static void
process_reads(const bool VERBOSE, const string &mapped_reads_file,
              const vector<string> &all_chroms,
              const vector<string> &chrom_names,
              const unordered_map<string, size_t> &chrom_lookup,
              ostream &out) {
  total_counts cnts;
  vector<double> psi;
  vector<pos_stats> stats;
  double similarity;

  bool chrom_has_reads = false;
  unordered_set<string> chroms_seen;
  string cur_chrom, cur_chrom_name = "";

  SAMReader sam_reader(mapped_reads_file);
  sam_rec aln;
  while (sam_reader >> aln) {
    if (aln.rname != cur_chrom_name) {
      // new chrom to process
      if (chroms_seen.find(aln.rname) != end(chroms_seen))
        throw runtime_error("chroms out of order in mapped reads file\n");

      if (chrom_has_reads) {
        summarize_chrom_stats(VERBOSE, cur_chrom, cur_chrom_name, 
                              psi, stats, cnts, similarity);

        write_chrom_stats(cur_chrom_name, cur_chrom.size(), 
                          similarity, psi, cnts, out);
      }

      if (VERBOSE)
        cerr << "[processing SAM reads in chromosome " << aln.rname << "]\n";

      cur_chrom_name = aln.rname;
      chroms_seen.insert(cur_chrom_name);
      get_chrom(cur_chrom_name, all_chroms, chrom_lookup, cur_chrom, stats, cnts);

      // reset statistics
      const size_t chrom_sz = cur_chrom.size();
      psi.resize(chrom_sz);
      for (size_t i = 0; i != chrom_sz; ++i)
        psi[i] = pos_stats::psi_not_calculated;
    }

    chrom_has_reads = true;
    adjust_alignment<bisulfite_bases>(aln);
    process_alignment<bisulfite_bases>(aln, cur_chrom, stats, cnts);
  }

  if (chrom_has_reads) {
    summarize_chrom_stats(VERBOSE, cur_chrom, cur_chrom_name,
                          psi, stats, cnts, similarity);
    write_chrom_stats(cur_chrom_name, cur_chrom.size(), 
                      similarity, psi, cnts, out);
  }
}

int
main(const int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    bool bisulfite = false;
    size_t num_threads = 1;
    string outfile;
    OptionParser opt_parse(strip_path(argv[0]),
                           "estimate edit distance between reference and mapped reads",
                           "<reference.fa> <seq-file.sam>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("threads", 't', "number of threads", false, num_threads);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("bisulfite", 'b', "reads come from WGBS", false, bisulfite);
    opt_parse.add_opt("output", 'o', "output YAML file", false, outfile);
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

    omp_set_num_threads(num_threads);
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

    std::ofstream of;
    if (!outfile.empty()) 
      of.open(outfile.c_str(), std::ios::binary);

    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    if (bisulfite)
      process_reads<true>(VERBOSE, mapped_reads_file,
                          all_chroms, chrom_names, chrom_lookup, out);
    else
      process_reads<false>(VERBOSE, mapped_reads_file,
                          all_chroms, chrom_names, chrom_lookup, out);

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
