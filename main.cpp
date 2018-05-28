/******************************************\
| Graph Splitting Method v1.0 (2018/06/11) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>
#include <unistd.h>
#include <random>   
#include <functional>

using namespace std;

/// format.cpp (File/Data-handling)
extern void readFASTA(ifstream&, ofstream&, ofstream&, int&);
extern void bl2mat(ifstream&, double*&, int const&);
extern void sc2nwk(int* const&, string&, int const&);
extern void addEP(string const&, string&, unordered_map<string, double>&, int const&, int const&);

/// mmseqs.cpp (Wrapper function of MMseqs)
extern void mmseqs(string const&, string const&, string const&);

/// gs.cpp (Core functions of GS method)
extern void GS(double* const&, int*&, int const&);
extern void EP(double* const&, unordered_map<string, double>&, function<double()>&, int const&);

/// messages.cpp
extern void print_banner();
extern void print_usage(char*&);

int main(int argc, char* argv[]){

  /*/ Getopt /*/
  int silence = 0;
  int ep_num  = 0;
  int seed    = 0;
  string threads = "1";
  int opt;
  regex renum(R"(^\d+$)"); // -e/-r/-t option requires an integer number
  opterr = 0;              // default error messages -> OFF

  while ((opt = getopt(argc, argv, "shve:r:t:")) != -1){
    if(opt == 'e'){ // OK! (./gs -e 100 IN.fst)
      if(regex_match(optarg, renum)){
	ep_num = atoi(optarg);
      }
      else{ // NG! (./gs -e hundred IN.fst)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'r'){ // OK! (./gs -r 12345 IN.fst)
      if(regex_match(optarg, renum)){
	seed = atoi(optarg);
      }
      else{ // NG! (./gs -r one_two_three IN.fst)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -r requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 't'){ // OK! (./gs -t 4 IN.fst)
      if(regex_match(optarg, renum)){
	auto th = atoi(optarg);
	if(th >0){
	  threads = string(optarg);
	}
	else{
	  /*PRINT*/ print_banner();
	  /*PRINT*/ cerr << "Option -t requires an integer argument (>=1).\n" << endl;
	  /*PRINT*/ print_usage(argv[0]);
	  return -1;
	}
      }
      else{ // NG! (./gs -t four IN.fst)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -t requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'h'){ // HELP message (./gs -h)
      /*PRINT*/ print_banner();
      /*PRINT*/ print_usage(argv[0]);
      return 0;
    }
    else if(opt == 'v'){ // Version (./gs -v)
      /*PRINT*/ print_banner();
      return 0;
    }
    else if(opt == 's'){ // SILENT mode (./gs -s -e 100 IN.fst)
      silence = 1;
    }
    else if (opt == '?'){
      if(optopt == 'e'){ // NG! (./gs IN.fst -e)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
      else{ // NG! (./gs -Z)
	/*PRINT*/ print_banner();
	/*PRINT*/ cerr << argv[0] << ": invalid option\n" <<  endl;
	/*PRINT*/ print_usage(argv[0]);
	return -1;
      }
    }
  }

  /*/ Input file name /*/
  string input = "";
  if(optind < argc){ // OK! (./gs -e 100 IN.fst)
    input = argv[optind];

    if(!silence){
      /*PRINT*/ print_banner();
      /*PRINT*/ cerr << "Number of threads used in MMseqs:\n  " << threads << endl;
      /*PRINT*/ cerr << "Number of iterations in EP method:\n  " << ep_num << endl;
      if(seed>0){
	/*PRINT*/ cerr << "Random seed number for EP method:\n  " << seed << endl;
      }
      else{
	/*PRINT*/ cerr << "Seed for the random number generator in EP method:\n  " << "a random number (default)" << endl;
      }
      /*PRINT*/ cerr << "Input file:\n  " << input << endl;
    }
  }
  else{ // NG! (./gs -e 100)
    /*PRINT*/ cerr << argv[0] << " requires an input file (fasta format).\n" << endl;
    /*PRINT*/ print_banner();
    /*PRINT*/ print_usage(argv[0]);
    return -1;
  }
  
  /*/ Variables /*/
  double* W;     // Sequence similarity matrix
  int size;      // Row size of W (W is a synmetry matrix)
  int* gs;       // Result of GS method
  string newick; // GS tree (without EP values)

  /*/ File I/O /*/
  regex re(R"(\.[^\.]+$)"); // Extention of the input file
  auto original_fasta = string(input);
  auto annotation_txt = regex_replace(original_fasta, re, "_annotation.txt");
  auto simple_fasta   = regex_replace(original_fasta, re, "_simple.fst");
  auto mmseqs_result  = regex_replace(original_fasta, re, "_mmseqs.txt");
  
  ifstream ifs1(original_fasta); // Fasta file (original)
  ofstream ofs1(annotation_txt); // Annotation file
  ofstream ofs2(simple_fasta);   // Fasta file (simple)
  
  if(ifs1.fail()){
    /*PRINT*/ cerr << "\nCannot access " << original_fasta << "!" << endl;
    return -1;
  }
  if(ofs1.fail()){
    /*PRINT*/ cerr << "\nCannot access " << annotation_txt << "!" << endl;
    return -1;
  }
  if(ofs2.fail()){
    /*PRINT*/ cerr << "\nCannot access " << simple_fasta << "!" << endl;    
    return -1;
  }

  /*/ Parsing fasta file /*/
  readFASTA(ifs1, ofs1, ofs2, size); 
    // ifs1: INPUT (original fasta file)
    // ofs1: OUTPUT (annotation file)
    // ofs2: OUTPUT (simplified fasta file)
    // size: # of sequence = row size of sequence similarity matrix

  ofs1.close();
  ofs2.close();
  
  /*PRINT*/ if(!silence) cerr << "Number of sequences:\n  " << size << " sequences" << endl;
  
  /*/ Executing MMSeqs /*/
  /*PRINT*/ if(!silence) cerr << "MMseqs:\n  " << size << "x" << size << " pairwise alignment\n" << "  searching...\r" << flush;
  mmseqs(simple_fasta, mmseqs_result, threads);
    // simple_fasta: INPUT (multiple fasta file) 
    // mmseqs_result: OUTOUT (result file of MMseqs)
    // threads: parameter (threads number for MMseqs)

  ifstream ifs2(mmseqs_result);  if(ifs2.fail()){return -1;} // MMseqs result file
  /*PRINT*/ if(!silence) cerr << "  completed!   " << endl;
  
  /*/ Reading Data /*/
  bl2mat(ifs2, W, size);
    // ifs2: INPUT (result file of MMseqs)
    // W: OUTPUT (sequence similarity matrix)

  /*/ GS method (stepwise spectral clustering) /*/
  /*PRINT*/ if(!silence) cerr << "GS method:\n" << "  searching...\r" << flush;
  GS(W, gs, size);
    // W: INPUT (sequence similarity matrix)
    // gs: OUTPUT (result of stepwise spectral clustering)

  /*PRINT*/ if(!silence) cerr << "  completed!   " << endl;

  /*/ Generating Newick file (GS tree WITHOUT EP values) /*/
  sc2nwk(gs, newick, size);
    // gs: INPUT (result of stepwise spectral clustering)
    // newick: OUTPUT (GS tree [newick format])

  if(ep_num>0){
    unordered_map<string, double> ep;
    string newick_EP; // GS+EP tree

    // Random number generator (Uniform distribution->Mersenne Twister)
    function<double()> R;
    uniform_real_distribution<double> urd(0,1);    // uniform distributed random number

    if(seed>0){
      mt19937 mt(static_cast<unsigned int>(seed)); // mersenne twister
      R = bind(urd, ref(mt));                      // random number generator    
    }
    else{
      random_device rd;                            // random seed
      mt19937 mt(rd());                            // mersenne twister
      R = bind(urd, ref(mt));                      // random number generator        
    }    

    /*PRINT*/ if(!silence) cerr << "EP method:" << endl;

    for(int n=1; n<=ep_num; n++){
      /*PRINT*/ if(!silence) cerr << "  " << n << "/" << ep_num << " iterations" << "\r"<< flush;
      EP(W, ep, R, size);
        // W: INPUT (sequence similarity matrix)
        // ep: OUTPUT (result of Edge Perturbation method)
        // R: random number generator
    }
    
    /*PRINT*/ if(!silence) cerr << "\n  completed!" << endl;
    /*PRINT*/ if(!silence) cerr << "\n--------------------------------------------------\n" << endl;

    addEP(newick, newick_EP, ep, ep_num, size);
      // newick: INPUT (GS tree [newick format])
      // newick_EP: OUTPUT (GS+EP tree [newick format])
      // ep: INPUT (result of Edge Perturbation method)
      // ep_num: INPUT (# of Edge Perturbation method)

    /*/ GS tree (WITH EP values) ->stdout /*/
    cout << newick_EP << flush;
  }
  else{
    /*PRINT*/ if(!silence) cerr << "\n--------------------------------------------------\n" << endl;
    
    /*/ GS tree (WITHOUT EP values) ->stdout /*/
    cout << newick << flush;
  }
  /*PRINT*/ cerr << endl;

  delete[] W;
  delete[] gs;

  return 0;
}
