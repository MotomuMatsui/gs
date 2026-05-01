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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>

#include "format.h" // also defines clade_key::hash_sorted_ids

using namespace std;

// Parsing multiple fasta file
int readFASTA(ifstream& ifs, ofstream& ofs1, ofstream& ofs2, int& row, bool& is_nucleotide){

  // 256-byte LUT for nucleotide/protein classification:
  //   class[c] == 1 → A/C/G/T/U/N (case-insensitive) → counts as both nucleotide
  //                   and letter
  //   class[c] == 2 → other alphabetic letter → counts only as letter
  //   class[c] == 0 → ignored (gap, digit, whitespace, punctuation)
  // One LUT lookup replaces an isalpha() + toupper() + 6-way char compare —
  // ~3-5× faster on long sequences and locale-independent.
  static const auto LUT = []{
    std::array<unsigned char, 256> t{};
    for(int c = 'A'; c <= 'Z'; c++) t[c] = 2;
    for(int c = 'a'; c <= 'z'; c++) t[c] = 2;
    for(char c : std::string("ACGTUN")) { t[static_cast<unsigned char>(c)] = 1; t[static_cast<unsigned char>(c) | 0x20] = 1; }
    return t;
  }();

  int info = 0;

  //Input & Output
  int id = 1;
  int f1 = 0;
  int f2 = 0;

  // Counters for nucleotide/amino-acid auto-detection.
  // Only A/C/G/T/U/N (case-insensitive) count as nucleotide characters;
  // any other alphabetic letter counts as non-nucleotide. Gaps, digits and
  // whitespace are ignored so aligned or numbered sequences don't skew the
  // ratio.
  long long nuc_count   = 0;
  long long letter_count = 0;

  string line;
  while(getline(ifs, line)){
    if(line[0] == '>'){
      if(f1==1){
	info = 1; //empty entry
	break;
      }
      else{
	auto name = line.substr(1);
	ofs1 << id << "\t" << name << "\n";
	ofs2 << '>' << id << "\n";

	id ++;
	f1 = 1;
	f2 = 1;
      }
    }
    else{
      if(f2==0){
	info = 2; //empty entry (1st entry lacked ">")
	break;
      }
      else{
	ofs2 << line << "\n";
	f1 = 0;

	for(unsigned char c : line){
	  unsigned char k = LUT[c];
	  if(k){
	    letter_count++;
	    if(k == 1) nuc_count++;
	  }
	}
      }
    }
  }

  // # of sequence (= row size of sequence similarity matrix)
  row = id - 1;

  // Classify as nucleotide when ≥ 90% of letters are A/C/G/T/U/N.
  // The threshold is intentionally strict: proteins are expected to have a
  // small but non-zero share of these letters (≈30% by chance), so requiring
  // 0.9 avoids mis-classifying protein input as DNA.
  if(letter_count > 0 && static_cast<double>(nuc_count) / letter_count >= 0.9){
    is_nucleotide = true;
  }
  else{
    is_nucleotide = false;
  }

  ofs1.close();
  ofs2.close();

  if(f1==1){
    info = 3; // empty entry (last entry lacked a sequence data)
  }
  else if(info==0 && row < 2){
    info = 4; // more than two sequences are required
  }

  return info;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// Parsing mmseqs result file
int bl2mat(ifstream& ifs, double* (&W), int const& size){

  //File I/O
  auto S_buf = std::unique_ptr<double[]>(new double[size*size]()); // Bit score matrix (RAII)
  double* S  = S_buf.get();
  W = new double[size*size](); // Sequence similarity matrix (ownership transfers to caller)

  //Reading lines
  int x = 1; // Sequence IDs (query)
  int y = 1; // Sequence IDs (target)
  double b = 0; // Bit scores calculated by MMseqs
  string line, chr;
  while(getline(ifs, line)){
    //Split lines
    istringstream stream(line);

    int pos = 0;
    while(getline(stream, chr, '\t')){
      if(pos == 0){
	x = stoi(chr) -1;
      } 
      else if(pos == 1){
	y = stoi(chr) -1;
      }
      else if(pos == 11){
	b = stod(chr);
      }
      pos ++;
    }

    auto p = x*size+y;
    auto s = S[p];
    S[p] = s>b? s:b;
  }  

  int same_sequence = 0;

  for(int x = 0; x < size; x ++){
    for(int y = x+1; y < size; y ++){
      auto xy = S[x*size+y];
      auto yx = S[y*size+x];      
      auto xx = S[x*size+x];
      auto yy = S[y*size+y];
      
      auto comp = xy+yx;
      auto self = xx+yy;
      auto sss  = self>0? comp/self:0;

      W[x*size+y] = sss;
      W[y*size+x] = sss;

      // Treat as duplicate when similarity is effectively 1 (avoid FP equality).
      if(sss >= 1.0 - 1e-9){
	same_sequence++;
      }
    }
    W[x*size+x] = 1;
  }

  return same_sequence;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// Parsing distance matrix file and converting to similarity matrix via exp(-d)
// Format: first line = N, then N lines of "label d0 d1 ... dN-1"
int distmat2mat(ifstream& ifs, ofstream& ofs_ann, double* (&W), int& size){

  string line;

  if(!getline(ifs, line)) return -1;
  {
    istringstream ss(line);
    if(!(ss >> size) || size < 2) return -1;
  }

  double* D = new double[size * size]();
  vector<string> labels(size);

  for(int i = 0; i < size; i++){
    if(!getline(ifs, line)){
      delete[] D;
      return -1;
    }
    istringstream row(line);
    if(!(row >> labels[i])){
      delete[] D;
      return -1;
    }
    for(int j = 0; j < size; j++){
      if(!(row >> D[i * size + j])){
        delete[] D;
        return -1;
      }
    }
  }

  for(int i = 0; i < size; i++){
    ofs_ann << (i + 1) << "\t" << labels[i] << "\n";
  }
  ofs_ann.close();

  W = new double[size * size]();
  for(int i = 0; i < size; i++){
    W[i * size + i] = 1.0;
    for(int j = i + 1; j < size; j++){
      double d   = (D[i * size + j] + D[j * size + i]) / 2.0;
      double sim = exp(-d);
      W[i * size + j] = sim;
      W[j * size + i] = sim;
    }
  }

  delete[] D;
  return 0;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// Parsing similarity matrix file (values used directly)
// Format: first line = N, then N lines of "label s0 s1 ... sN-1"
int simmat2mat(ifstream& ifs, ofstream& ofs_ann, double* (&W), int& size){

  string line;

  if(!getline(ifs, line)) return -1;
  {
    istringstream ss(line);
    if(!(ss >> size) || size < 2) return -1;
  }

  double* S = new double[size * size]();
  vector<string> labels(size);

  for(int i = 0; i < size; i++){
    if(!getline(ifs, line)){
      delete[] S;
      return -1;
    }
    istringstream row(line);
    if(!(row >> labels[i])){
      delete[] S;
      return -1;
    }
    for(int j = 0; j < size; j++){
      if(!(row >> S[i * size + j])){
        delete[] S;
        return -1;
      }
    }
  }

  for(int i = 0; i < size; i++){
    ofs_ann << (i + 1) << "\t" << labels[i] << "\n";
  }
  ofs_ann.close();

  W = new double[size * size]();
  for(int i = 0; i < size; i++){
    W[i * size + i] = 1.0;
    for(int j = i + 1; j < size; j++){
      double sim = (S[i * size + j] + S[j * size + i]) / 2.0;
      W[i * size + j] = sim;
      W[j * size + i] = sim;
    }
  }

  delete[] S;
  return 0;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// Parsing result file generated by GS method
void sc2nwk(int* const& W, string& newick, int const& size){

  //Input & Output
  int* tr     = new int[size+1](); 
  int* change = new int[size*3]();

  //Reading lines
  for(int r=0; r<size; r++){
    int cur = 1;
    int old = 1;
    for(int c=0; c<size; c++){
      cur = W[r*size+c];
      if(cur != old){
	change[c*3]   = old; // previous membership
	change[c*3+1] = cur; // current membership
	change[c*3+2] ++;    // # of cluster member
      }
      old = cur;
    }
    tr[W[(r+1)*size-1]] = r+1;
  }

  int* nwk = new int[size+1](); for(int p=0; p<size; p++){nwk[p] = 1;}
  int* bra = new int[2*(size+1)](); // # of brachets
  int* pos = new int[2*(size+1)](); // left and right boundary of each cluster
  int* com = new int[size+1]();     // position of commas

  bra[0]        = 1;
  bra[size*2+1] = 1;
  pos[2]        = 0;
  pos[3]        = size;

  for(int p=1; p<size; p++){
    auto old = change[p*3];
    auto cur = change[p*3+1];
    auto num = change[p*3+2];

    auto ps  = pos[old*2];
    auto ls  = pos[old*2+1];
    
    pos[old*2]   += num; // 1 1 1 1 1 -> (2 2),(1 1 1) ... left boundary of cluster 1 = 2
    pos[old*2+1] -= num; // 1 1 1 1 1 -> (2 2),(1 1 1) ... # of cluster 1 = 3
    pos[cur*2]   = ps;   // 1 1 1 1 1 -> (2 2),(1 1 1) ... left boundary of cluster 2 = 0
    pos[cur*2+1] = num;  // 1 1 1 1 1 -> (2 2),(1 1 1) ... # of cluster 2 = 2
    com[ps+num]  = 1;    // 1 1 1 1 1 -> (2 2),(1 1 1) ... position of comma = 2
    
    for(int n=0; n<num; n++){
      nwk[ps+n] = cur;   // 1 1 1 1 1 -> (2 2),(1 1 1) ... replacing 1 -> 2
    }
    
    bra[ps*2]         ++; // (2 2
    bra[(ps+num)*2+1] ++; //     )
    bra[(ps+num)*2]   ++; //      (1 1 1
    bra[(ps+ ls)*2+1] ++; //            )
  }

  newick = "";
  for(int b=0; b<=size; b++){
    bra[b*2]   --; // (
    bra[b*2+1] --; // )

    string right = (bra[b*2+1]>0)? string(bra[b*2+1], ')') : "";
    string left  = (bra[b*2]  >0)? string(bra[b*2],   '(') : "";
    string comma = (com[b]    >0)? string(com[b],     ',') : "";
    string node  = (tr[nwk[b]]>0)? to_string(tr[nwk[b]])   : "";

    newick += right+comma+left+node;
  }
  newick += ";";

  //Free
  delete[] tr;
  delete[] change;
  delete[] nwk;
  delete[] bra;
  delete[] pos;
  delete[] com;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void addEP(string const& newick, string& newick_EP, unordered_map<std::uint64_t, double>& ep, int const& ep_num, int const& size){

  newick_EP = "";
  int N = newick.size();
  int B = 0; //brachet
  int f = 1;
  string id = "";

  unordered_map<int, vector<int>> clade;
  stringstream ss_EP;

  for(int n=0; n<N; n++){
    auto p = newick[n];
    ss_EP << p;

    if(p == '('){
      B++;
    }
    else if(p == ')'){

      if(id.size() > 0){
	for(int depth=1; depth<=B; depth++){
	  clade[depth].push_back(stoi(id));
	}
	id = "";
      }

      // Hash the sorted clade ID list — same key space as EP_fbs/EP_tbe.
      sort(clade[B].begin(), clade[B].end());
      std::uint64_t const C = clade_key::hash_sorted_ids(clade[B]);
      int const clade_size  = static_cast<int>(clade[B].size());

      if(B>2){
	ss_EP << fixed << setprecision(2) << ep[C]/ep_num;
      }
      else if(B==2 && f==1){
	// Original guarded against the outermost bipartition where one side
	// is the root-clade containing all but one taxon.  Now expressed
	// directly via clade size instead of counting '|' separators.
	if(clade_size < size-1){
	  ss_EP << fixed << setprecision(2) << ep[C]/ep_num;
	}

	f=0;
      }

      clade[B] = {};
      B--;
    }
    else if(p == ','){
      if(id.size() > 0){
	for(int depth=1; depth<=B; depth++){
	  clade[depth].push_back(stoi(id));
	}
	id = "";
      }
    }
    else{
      id += p;
    }
  }

  newick_EP = ss_EP.str();
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void addLABEL(string const& newick, string& newick_ann, string const& annotation_txt, int const& size){

  //File I/O
  ifstream ifs(annotation_txt);   // Annotation file
  std::vector<std::string> L(size); // LABEL: ID->Annotation (RAII)

  string line, chr;
  while(getline(ifs, line)){

    int pos = 0;
    int id  = 0;
    string ann;
    //Split lines
    istringstream stream(line);
    while(getline(stream, chr, '\t')){
      if(pos == 0){
        id = stoi(chr)-1;
      } 
      else if(pos == 1){
        ann = string(chr);
	break;
      }
      pos ++;
    }

    L[id] = ann;
  }

  newick_ann = "";
  int N = newick.size();
  int s = 0;
  string id = "";
  stringstream ss_ann;

  for(int n=0; n<N; n++){
    auto p = newick[n];

    if(p == '('){
      ss_ann << p;

      s = 1;
    }
    else if(p == ')'){
      if(id.size() > 0){
	ss_ann << '"' << L[stoi(id)-1] << '"';
	id = "";
      }
      ss_ann << p;
      s = 0;
    }
    else if(p == ','){
      if(id.size() > 0){
	ss_ann << '"' << L[stoi(id)-1] << '"';
	id = "";
      }
      ss_ann << p;
      s = 1;
    }
    else if(s==1){
      id += p;
    }
    else{
      ss_ann << p;
    }
  }

  newick_ann = ss_ann.str();
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void sc2list(int* const& gs, int* (&list), int const& size){

  list = new int[(size-3)*size]();

  int l = 0;
  for(int i=1; i<size; i++){
    auto m   = i+1;
    int flag = -1;

    vector<int> a;
    vector<int> b;
    
    for(int j=0; j<size; j++){
      if(gs[i+j*size] == m){
        flag = gs[i+j*size-1];
        break;
      }
    }

    for(int j=0; j<size; j++){
      if(gs[i+j*size] == m){
        a.push_back(j);
      }
      else if(gs[i+j*size] == flag){
        b.push_back(j);
      }
    }

    int a_size = a.size();
    int b_size = b.size();

    if(a_size == size-1) continue;

    if(a_size>1){
      for(int n : a){
        list[l*size+n]=1;
      }
      l ++;
    }

    if(i==1) continue;

    if(b_size>1){
      for(int n : b){
        list[l*size+n]=1;
      }
      l ++;
    }
  }
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
