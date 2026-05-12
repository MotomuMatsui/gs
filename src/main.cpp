/******************************************\
| Graph Splitting Method v2.5 (2026/04/23) |
|                                          |
| Copyright (c) Motomu Matsui              |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2019)        |
|     Systematic Biology, 69(2):265-279.   |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <random>
#include <unistd.h>

#include "blastn.h"
#include "ep.h"
#include "gs.h"
#include "format.h"
#include "messages.h"
#include "mmseqs.h"
#include "transitivity.h"

using namespace std;

int main(int argc, char* argv[]){

  // Pin BLAS to one thread per process before any BLAS call.  gs2's
  // parallelism lives at the EP / OMP level; allowing multi-threaded
  // OpenBLAS underneath an OMP-parallel call site oversubscribes the
  // CPU (N_omp × N_blas threads) and tanks throughput on 64-core hosts.
  // setenv with overwrite=0 respects an explicit user override.
  setenv("OPENBLAS_NUM_THREADS", "1", 0);

  int silence   = 0;
  int ep_num    = 0;
  int seed      = 0;
  int label     = 0;
  int dist_mode = 0;
  int sim_mode  = 0;
  string threads = "1";
  string sensitivity = "7.5";
  string bs_method = "";

  opterr = 0;
  int opt;
  // Integer-only options (-e, -r, -t) must reject decimals; -m accepts floats.
  regex reint(R"(^\d+$)");
  regex reflt(R"(^\d+(\.\d+)?$)");
  while ((opt = getopt(argc, argv, "shlve:r:t:m:b:dw")) != -1){
    if(opt == 'e'){
      if(regex_match(optarg, reint)){
	ep_num = atoi(optarg);
      }
      else{
	print_banner();
	cerr << "Option -e requires an integer argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'r'){
      if(regex_match(optarg, reint)){
	seed = atoi(optarg);
      }
      else{
	print_banner();
	cerr << "Option -r requires an integer argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 't'){
      if(regex_match(optarg, reint)){
	auto th = atoi(optarg);
	if(th >0){
	  threads = string(optarg);
	}
	else{
	  print_banner();
	  cerr << "Option -t requires an integer argument (>=1).\n" << endl;
	  print_usage(argv[0]);
	  return -1;
	}
      }
      else{
	print_banner();
	cerr << "Option -t requires an integer argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'm'){
      if(regex_match(optarg, reflt)){
	auto sen = atof(optarg);
	if(1<=sen && sen<=7.5){
	  sensitivity = string(optarg);
	}
	else{
	  print_banner();
	  cerr << "Option -m requires a double number argument [1, 7.5].\n" << endl;
	  print_usage(argv[0]);
	  return -1;
	}
      }
      else{
	print_banner();
	cerr << "Option -m requires a double number argument [1, 7.5].\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'h'){
      print_banner();
      print_usage(argv[0]);
      return 0;
    }
    else if(opt == 'v'){
      print_banner();
      return 0;
    }
    else if(opt == 's'){
      silence = 1;
    }
    else if(opt == 'l'){
      label = 1;
    }
    else if(opt == 'd'){
      dist_mode = 1;
    }
    else if(opt == 'w'){
      sim_mode = 1;
    }
    else if(opt == 'b'){
      bs_method = optarg;
      if(bs_method != "fbs" && bs_method != "tbe"){
        print_banner();
        cerr << "Option -b requires a string that equals 'fbs' or 'tbe'.\n" << endl;
        print_usage(argv[0]);
        return -1;
      }
    }
    else if (opt == '?'){
      if(optopt == 'e'){
	print_banner();
	cerr << "Option -e requires an integer argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
      else if(optopt == 'r'){
	print_banner();
	cerr << "Option -r requires an integer argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
      else if(optopt == 't'){
	print_banner();
	cerr << "Option -t requires an integer argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
      else if(optopt == 'm'){
	print_banner();
	cerr << "Option -m requires a flout argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
      else if(optopt == 'b'){
	print_banner();
	cerr << "Option -b requires a string argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
      else{
	print_banner();
	cerr << argv[0] << ": invalid option\n" <<  endl;
	print_usage(argv[0]);
	return -1;
      }
    }
  }

  string input = "";
  if(optind < argc){
    input = argv[optind];
  }
  else{
    print_banner();
    cerr << argv[0] << " requires an input file.\n" << endl;
    print_usage(argv[0]);
    return -1;
  }

  double* W_raw = nullptr;
  int     size  = 0;
  int*    gs_raw = nullptr;
  string  newick;

  regex re(R"(\.[^\.]+$)");
  auto original_input = string(input);
  auto annotation_txt = regex_replace(original_input, re, "_annotation.txt");
  string simple_fasta;
  string search_result;

  if(dist_mode && sim_mode){
    print_banner();
    cerr << "Options -d and -w cannot be used together.\n" << endl;
    print_usage(argv[0]);
    return -1;
  }
  else if(dist_mode || sim_mode){

    ifstream ifs(original_input);
    ofstream ofs_ann(annotation_txt);

    if(ifs.fail()){
      cerr << "\nCannot access " << original_input << "!" << endl;
      return -1;
    }
    if(ofs_ann.fail()){
      cerr << "\nCannot create " << annotation_txt << "!" << endl;
      return -1;
    }

    string input_type;
    int parse_result;
    if(dist_mode){
      input_type   = "distance matrix";
      parse_result = distmat2mat(ifs, ofs_ann, W_raw, size);
    }
    else{
      input_type   = "similarity matrix";
      parse_result = simmat2mat(ifs, ofs_ann, W_raw, size);
    }

    if(parse_result != 0){
      print_banner();
      cerr << "Failed to parse " << input_type << " from " << original_input << ".\n" << endl;
      print_usage(argv[0]);
      return -1;
    }

    if(!silence){
      print_banner();
      cerr << "Settings:" << endl;
      cerr << "-Input (" << input_type << ")" << endl;
      cerr << "  file = " << input << endl;
      cerr << "  # of sequences = " << size << endl << endl;

      cerr << "-EP method" << endl;
      if(ep_num>0){
        if(bs_method == "fbs"){
          cerr << "  method = Felsenstein's bootstrap proportion" << endl;
        }
        else{
          cerr << "  method = Transfer bootstrap expectation (F. Lemoine, et al., Nature, 2018)" << endl;
        }
      }

      if(seed>0){
        cerr << "  random seed = " << seed << endl;
      }
      else{
        cerr << "  random seed = " << "a random number (default)" << endl;
      }
      cerr << "  # of iterations = " << ep_num << endl << endl;
      cerr << "Progress:" << endl;
    }
  }
  else{

    simple_fasta  = regex_replace(original_input, re, "_simple.fst");
    search_result = regex_replace(original_input, re, "_mmseqs.txt");

    ifstream ifs1(original_input);
    ofstream ofs1(annotation_txt);
    ofstream ofs2(simple_fasta);

    if(ifs1.fail()){
      cerr << "\nCannot access " << original_input << "!" << endl;
      return -1;
    }
    if(ofs1.fail()){
      cerr << "\nCannot create " << annotation_txt << "!" << endl;
      return -1;
    }
    if(ofs2.fail()){
      cerr << "\nCannot create " << simple_fasta << "!" << endl;
      return -1;
    }

    bool is_nucleotide = false;
    auto file_condition = readFASTA(ifs1, ofs1, ofs2, size, is_nucleotide);

    ofs1.close();
    ofs2.close();

    // Pick the aligner based on the detected sequence type: blastn for
    // nucleotide input, mmseqs2 for protein input. Verify the selected
    // tool is on PATH before proceeding.
    auto required_tool = is_nucleotide ? string("blastn") : string("mmseqs");
    // `&>` is a bash extension; under POSIX /bin/sh (dash on Debian/Ubuntu)
    // it parses as `&` then `>`, leaking the `which` stdout into our output.
    auto which_cmd     = "which " + required_tool + " > /dev/null 2>&1";
    auto status        = system(which_cmd.c_str());
    if(WEXITSTATUS(status) != 0){
      print_banner();
      cerr << required_tool << " not found!\nCheck our web page (https://github.com/MotomuMatsui/gs) for more information" << endl;
      return -1;
    }

    if(file_condition == 1){
      print_banner();
      cerr << "Input file had an empty entry.\n" << endl;
      print_usage(argv[0]);
      return -1;
    }
    else if(file_condition == 2){
      print_banner();
      cerr << "Input file had an empty entry (1st entry lacked '>').\n" << endl;
      print_usage(argv[0]);
      return -1;
    }
    else if(file_condition == 3){
      print_banner();
      cerr << "Input file had an empty entry (last entry lacked '>').\n" << endl;
      print_usage(argv[0]);
      return -1;
    }
    else if(file_condition == 4){
      print_banner();
      cerr << "Input file have to contain more than two sequences.\n" << endl;
      print_usage(argv[0]);
      return -1;
    }

    if(!silence){
      print_banner();
      cerr << "Settings:" << endl;
      cerr << "-Input" << endl;
      cerr << "  file = " << input << endl;
      cerr << "  type = " << (is_nucleotide ? "nucleotide" : "amino acid") << " (auto-detected)" << endl;
      cerr << "  # of sequences = " << size << endl << endl;

      if(is_nucleotide){
        cerr << "-BLASTN" << endl;
        cerr << "  # of threads = " << threads << endl << endl;
      }
      else{
        cerr << "-MMseqs" << endl;
        cerr << "  sentitivity = " << sensitivity << endl;
        cerr << "  # of threads = " << threads << endl << endl;
      }

      cerr << "-EP method" << endl;
      if(ep_num>0){
        if(bs_method == "fbs"){
          cerr << "  method = Felsenstein's bootstrap proportion" << endl;
        }
        else{
          cerr << "  method = Transfer bootstrap expectation (F. Lemoine, et al., Nature, 2018)" << endl;
        }
      }

      if(seed>0){
        cerr << "  random seed = " << seed << endl;
      }
      else{
        cerr << "  random seed = " << "a random number (default)" << endl;
      }
      cerr << "  # of iterations = " << ep_num << endl << endl;
      cerr << "Progress:" << endl;
    }

    if(is_nucleotide){
      if(!silence) cerr << "-BLASTN\n  " << size << "x" << size << " pairwise alignment\n" << "  searching...\r" << flush;
      blastn(simple_fasta, search_result, threads);
    }
    else{
      if(!silence) cerr << "-MMseqs\n  " << size << "x" << size << " pairwise alignment\n" << "  searching...\r" << flush;
      mmseqs(simple_fasta, search_result, threads, sensitivity);
    }

    ifstream ifs2(search_result);
    if(ifs2.fail()){
      cerr << "  Cannot access " << search_result << "!" << endl;
      return -1;
    }
    else{
      if(!silence) cerr << "  done.        " << endl << endl;
    }

    int same_sequence = bl2mat(ifs2, W_raw, size);
    if(!silence) if(same_sequence>0) cerr << "  <WARNING> This dataset has " << same_sequence << " duplicated sequence pair(s)" << endl << endl;
  }

  auto W = unique_ptr<double[]>(W_raw); // adopt ownership

  auto transitivity_score = transitivity(W.get(), size);
  if(!silence) cerr << "  transitivity = "<< transitivity_score << endl << endl;

  if(!silence) cerr << "-GS method\n" << "  executing...\r" << flush;
  GS(W.get(), gs_raw, size);
  auto gs = unique_ptr<int[]>(gs_raw); // adopt ownership from GS
  if(!silence) cerr << "  done.         " << endl << endl;

  sc2nwk(gs.get(), newick, size);

  if(ep_num>0){
    if(!silence) cerr << "-EP method\n  executing " << ep_num << " iterations (parallelised)...\r" << flush;

    // Derive a concrete seed when the user did not specify one
    int effective_seed = (seed > 0) ? seed : static_cast<int>(random_device{}());

    unordered_map<std::uint64_t, double> ep;
    string newick_EP;

    if(bs_method == "fbs"){
      EP_fbs(W.get(), ep, effective_seed, ep_num, size);
    }
    else{
      int* list_ori_raw = nullptr;
      sc2list(gs.get(), list_ori_raw, size);
      auto list_ori = unique_ptr<int[]>(list_ori_raw); // adopt ownership from sc2list

      EP_tbe(W.get(), list_ori_raw, ep, effective_seed, ep_num, size);
    }

    if(!silence) cerr << "\n  done." << endl << endl;
    if(!silence) cerr << "------------------------------------------\n" << endl;

    addEP(newick, newick_EP, ep, ep_num, size);

    if(label==1){
      string newick_ann;
      addLABEL(newick_EP, newick_ann, annotation_txt, size);
      cout << newick_ann << endl;
    }
    else{
      cout << newick_EP << endl;
    }
  }
  else{
    if(!silence) cerr << "------------------------------------------\n" << endl;

    if(label==1){
      string newick_ann;
      addLABEL(newick, newick_ann, annotation_txt, size);
      cout << newick_ann << endl;
    }
    else{
      cout << newick << endl;
    }
  }

  std::remove(annotation_txt.c_str());
  if(!simple_fasta.empty())  std::remove(simple_fasta.c_str());
  if(!search_result.empty()) std::remove(search_result.c_str());

  return 0; // W and gs are released automatically by unique_ptr
}
