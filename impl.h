#ifndef NUCLEOTIDE_UTILS_H
#define NUCLEOTIDE_UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>  // For std::max
using namespace std;

// Global variable declaration
extern vector<vector<int>> fr;

// Function declarations
void frq(string spacer);
void checklen(string spacer, vector<int> &freq);
void count(string s, vector<int> &v);
void helper(vector<string> spacer, ofstream &output, ofstream &ott, 
            ofstream &sp_data, string A, string B, int o);
void findPromoter(const string &genomeSequence, vector<string> &promoter);
void findSpacers(const string &PromoterSequence, vector<string> &spacers, 
                 int &count, string a, string b, vector<string> &chrom_id, 
                 vector<int> &spacer_start, vector<int> &spacer_end,
                 const string &currentChrom, int offset);
vector<string> findConservedSequences(vector<string> &spacers, double threshold);

#endif // NUCLEOTIDE_UTILS_H
