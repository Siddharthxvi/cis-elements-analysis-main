#include "impl.h"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
vector<vector<int> > fr(4, vector<int>(30, 0));

void frq(string spacer)
{
    for (int i = 0; i < spacer.size(); i++)
    {
        if (spacer[i] == 'A')
            fr[0][i]++;
        if (spacer[i] == 'C')
            fr[1][i]++;
        if (spacer[i] == 'G')
            fr[2][i]++;
        if (spacer[i] == 'T')
            fr[3][i]++;
    }
}
void checklen(string spacer, vector<int> &freq)
{
    freq[spacer.size()]++;
}

void count(string s, vector<int> &v)
{
    for (auto i : s)
    {
        if (i == 'A')
            v[0]++;
        if (i == 'C')
            v[1]++;
        if (i == 'G')
            v[2]++;
        if (i == 'T')
            v[3]++;
    }
}
//void helper(ifstream &input, ofstream &output, ofstream &ott, ofstream &sp_data, string A, string B, int o)
void helper(vector<string> spacer, ofstream &output, ofstream &ott, ofstream &sp_data, string A, string B, int o)
{

    int n = 4;
    vector<int> v(n);
    vector<int> freq(31);
    // ifstream input("spacer.txt");
    // ofstream output("final_output.txt");

   // string line3;
    for(auto it : spacer){
        count(it,v);
        if (isspace(it[0]))
            freq[0]++;
        checklen(it, freq);
        frq(it);
    }
    // while (getline(input, line3))
    // {
    //     count(line3, v);
    //     if (isspace(line3[0]))
    //         freq[0]++;
    //     checklen(line3, freq);
    //     frq(line3);
    // }

    for (int i = 0; i < 4; i++)
    {
        if (i == 0)
            output << "Occurence of A Nucleotide : ";
        if (i == 1)
            output << "Occurence of C Nucleotide : ";
        if (i == 2)
            output << "Occurence of G Nucleotide : ";
        if (i == 3)
            output << "Occurence of T Nucleotide : ";

        output << v[i] << endl;
    }
    output << endl;
    // output << "Total spacer sequences found for orientation "<<A <<"-" <<B<<" : " << total_count << endl<<endl;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0){
             ott << "A "<<endl;
            output << "A "<<endl;
        }
           
        if (i == 1){
            ott << "C "<<endl;
            output << "C "<<endl;
        }
            
        if (i == 2){
            ott << "G "<<endl;
            output << "G "<<endl;
        }
            
        if (i == 3){
               ott << "T "<<endl;
            output << "T "<<endl;
        }
         
        for (int j = 0; j < 30; j++)
        {
            output << " " << j+1 << " : " << fr[i][j] << " ";
            ott << fr[i][j]<<endl;

        }
        output << endl;
    }
    output << endl;

    int a = 0, max_freq = 0;

     sp_data << A+"-"+B << endl;
   
    for (int i = 0; i <= 30; i++)
    {
       output << "Number of sequences with Spacer Length : " << i << " = " << freq[i] << endl;
       sp_data << freq[i] << endl;
        max_freq = max(max_freq, freq[i]);
        a = (freq[i] >= max_freq) ? i : a;
    }
    output << endl;
    output << "Most common spacer length is : " << a << endl
           << endl;
    output << "__________________________________________________________________________________________________________________"<<endl;
}

//Finds conserved spacer sequence for each spacer length
vector<string> findConservedSequences(vector<string> &spacers, double threshold)
{
    vector<string> consensusSequences(21, "");
    map<int, vector<string>> groupedSequences;

    for (const string& seq : spacers) {
        int len = seq.size();
        if (len >= 1 && len <= 21){
            groupedSequences[len].push_back(seq);
        }
    }

    for (auto it = groupedSequences.begin(); it != groupedSequences.end(); ++it) {
        int length = it->first;
        vector<string> sequences = it->second;
        int numSequences = sequences.size();
        vector<unordered_map<char, int>> positionCounts(length);

        // Count nucleotide occurrences at each position
        for (const string& seq : sequences) {
            for (int i = 0; i < length; ++i) {
                positionCounts[i][seq[i]]++;
            }
        }

        string consensus;
        for (int i = 0; i < length; ++i) {
            char conservedNucleotide = '-';
            int maxCount = 0;

            // Most frequent nucleotide at the position
            for (const auto& pair : positionCounts[i]) {
                if (pair.second > maxCount) {
                    maxCount = pair.second;
                    conservedNucleotide = pair.first;
                }
            }

            // Conservation threshold
            if ((double)maxCount / numSequences >= threshold) {
                consensus += conservedNucleotide;
            } else {
                consensus += '-';
            }
        }

        consensusSequences[length - 1] = consensus;
    }

    return consensusSequences;
}
