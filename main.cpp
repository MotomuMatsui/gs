/*=========================================*\
| Graph Splitting Method v1.0 (2018/06/11) |
|                                          |
| Copyright (c) 2018 Motomu Matsui         |
| Systematic Biology, xx:xx-xx, 2018       |
|                                          |
| Web:  http://gs.bs.s.u-tokyo.ac.jp/      |
| Mail: qm at bs.s.u-tokyo.ac.jp           |
\*=========================================*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <regex>
#include <unistd.h>

using namespace std;

//format.cpp (File/Data-handling)
extern void readFASTA(ifstream&, ofstream&, ofstream&, int&);
extern void bl2mat(ifstream&, double*&, int const&);
extern void sc2nwk(int* const&, string&, int const&);
extern void addEP(string const&, string&, unordered_map<string, double>&, int const&, int const&);

//mmseqs.cpp (Wrapper function of MMseqs)
extern void mmseqs(string const&, string const&);

//gs.cpp (Core functions of GS method)
extern void GS(double* const&, int*&, int const&);
extern void EP(double* const&, unordered_map<string, double>&, int const&);

//messages.cpp
extern void print_banner();
extern void print_usage(char*&);


int main(int argc, char* argv[]){

  // Getopt
  int silence = 0;
  int ep_num  = 0;
  int opt;
  regex renum(R"(^\d+$)"); // -e option requires an integer number
  opterr = 0; //default error messages -> off
  while ((opt = getopt(argc, argv, "she:")) != -1){
    if(opt == 'e'){ // OK ($ ./gs -e 100 IN.fst > OUT.nwk)
      if(regex_match(optarg, renum)){
	ep_num = atoi(optarg);	
      }
      else{ // NG ($ ./gs -e hundred IN.fst > OUT.nwk)
	print_banner();
	cerr << "Option -e requires an integer argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
    }
    else if(opt == 'h'){ // HELP message ($ ./gs -h)
      print_banner();
      print_usage(argv[0]);
      return -1;
    }
    else if(opt == 's'){ // SILENT mode ($ ./gs -s -e 100 IN.fst > OUT.nwk)
      silence = 1;
    }
    else if (opt == '?'){
      if(optopt == 'e'){ // NG ($ ./gs IN.fst -e > OUT.nwk)
	print_banner();
	cerr << "Option -e requires an integer argument.\n" << endl;
	print_usage(argv[0]);
	return -1;
      }
      else{ // NG ($ ./gs -Z)
	print_banner();
	cerr << argv[0] << ": invalid option\n" <<  endl;
	print_usage(argv[0]);
	return -1;
      }
    }
  }

  // Get input file name
  string input = "";
  if(optind < argc){ // OK ($ ./gs -e 100 IN.fst > OUT.nwk)
    input = argv[optind];

    if(!silence){
      print_banner();
      cerr << "Number of Edge Perturbating:\n  " << ep_num << endl;
      cerr << "Input file:\n  " << input << endl;
    }
  }
  else{ // NG ($ ./gs -e 100 > OUT.nwk)
    cerr << argv[0] << " requires an input file (fasta format).\n" << endl;
    print_banner();
    print_usage(argv[0]);
    return -1;
  }
  
  // Variables
  double* W;     // Sequence similarity matrix
  int size;      // Row size of W (W is a synmetry matrix)
  int* gs;       // Result of GS method
  string newick; // GS tree (without EP values)

  // File I/O
  regex re(R"(\.[^\.]+$)"); // Extention of the input file

  auto original_fasta = string(input);
  auto annotation_txt = regex_replace(original_fasta, re, "_annotation.txt");
  auto simple_fasta   = regex_replace(original_fasta, re, "_simple.fst");
  auto mmseqs_result  = regex_replace(original_fasta, re, "_mmseqs.txt");
  
  ifstream ifs1(original_fasta); if(ifs1.fail()){return -1;} // Fasta file (original)
  ofstream ofs1(annotation_txt); if(ofs1.fail()){return -1;} // Annotation file
  ofstream ofs2(simple_fasta);   if(ofs2.fail()){return -1;} // Fasta file (simple)
  
  // Parsing fasta file
  // original fasta file (ifs1) -> annotation file (ofs1) & simplified fasta file (ofs2)
  // size: # of sequence = row size of sequence similarity matrix
  readFASTA(ifs1, ofs1, ofs2, size); 
  ofs1.close();
  ofs2.close();

  if(!silence){
    cerr << "Number of sequences:\n  " << size << " sequence(s)" << endl;
  }

  // Executing MMSeqs
  // Input: Multiple fasta file (simple_fasta)
  // Output: Result file of MMseqs (mmseqs_result)
  if(!silence){
    cerr << "MMseqs:\n" << "  searching...\r" << flush;
  }

  mmseqs(simple_fasta, mmseqs_result);
  ifstream ifs2(mmseqs_result);  if(ifs2.fail()){return -1;} // MMseqs result file

  if(!silence){
    cerr << "  completed!   " << endl;
  }

  // Reading Data
  // Input: Result file of MMseqs (mmseqs_result)
  // Output: Sequence similarity matrix (W)
  bl2mat(ifs2, W, size);

  // GS method (stepwise spectral clustering)
  // Input: Sequence similarity matrix (W)
  // Output: Result of stepwise spectral clustering (gs)
  if(!silence){
    cerr << "GS method:\n" << "  searching...\r" << flush;
  }

  GS(W, gs, size);

  if(!silence){
    cerr << "  completed!   " << endl;
  }

  // Generating Newick file (GS tree WITHOUT EP values)
  // Input: Result of stepwise spectral clustering (gs)
  sc2nwk(gs, newick, size);

  if(ep_num>0){
    unordered_map<string, double> ep;
    string newick_EP; // GS+EP tree

    if(!silence){
      cerr << "EP method:" << endl;
    }

    for(int n=1; n<=ep_num; n++){
      if(!silence){
	cerr << "  " << n << "/" << ep_num<< "\r"<< flush;
      }
      EP(W, ep, size);
    }
    
    if(!silence){
      cerr << "\n  completed!" << endl;
      cerr << "\n--------------------------------------------------\n" << endl;
    }

    addEP(newick, newick_EP, ep, ep_num, size);

    //GS tree (WITH EP values) -> stdout
    cout << newick_EP << flush;
  }
  else{
    if(!silence){
      cerr << "\n--------------------------------------------------\n" << endl;
    }
    
    //GS tree (WITHOUT EP values) -> stdout
    cout << newick << flush;
  }
  cerr << endl;

  // Free
  delete[] W;
  delete[] gs;

  // Return
  return 0;
}
