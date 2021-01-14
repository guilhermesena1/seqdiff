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

using std::runtime_error;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::begin;
using std::end;
using std::unordered_map;
using std::unordered_set;

static bool file_exists(const string &file) {
  return access(file.c_str(), R_OK) == 0;
}

int
main(const int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    OptionParser opt_parse(strip_path(argv[0]),
                           "estimate edit distance between reference and mapped reads",
                           "<reference.fa> <seq-file.sam>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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

    // program starts here
    const string genome_file = leftover_args.front();
    const string mapped_reads_file = leftover_args.back();

    if (!file_exists(genome_file))
      throw runtime_error("cannot open genome file: " + genome_file);

    if (!file_exists(mapped_reads_file))
      throw runtime_error("cannot open mapped reads file: " +
          mapped_reads_file);

    SAMReader sam_reader(mapped_reads_file);
    sam_rec aln;

    // read chroms into strings
    vector<string> all_chroms;
    vector<string> chrom_names;
    unordered_map<string, size_t> chrom_lookup;
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
        cerr << "  chrom: " << chrom_names.back() << "\tsize: " <<
                  all_chroms.back().size() << "\n";
    }

    if (VERBOSE)
      cerr << "[read " << all_chroms.size() << " chromosomes]\n";

  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }

    return 0;
}
