#include "impl.h"
#include <string>   // For std::string
#include <vector> 
void findPromoter(const string &genomeSequence, vector<string> &promoter)
{
    size_t len = genomeSequence.size();
    size_t searchStart = 0; // Start searching from the beginning of the genome sequence
    while (true)
    {
        // Find the next occurrence of TATA sequence starting from searchStart
        size_t N = genomeSequence.find("TATA", searchStart);
        if (N == string::npos)
            break;

        // Check if TATA sequence is beyond the first 1000 bases
        if (N > 1000)
        {
            // Search for ATG sequence 28-35 bases downstream
            for (int i = N + 28; i <= N + 35; ++i)
            {
                if (i >= 0 && genomeSequence.substr(i, 3) == "ATG")
                {
                     //ot << "Start : " << N-1000 << " End : "<< N <<endl;
                    //outputFile << genomeSequence.substr(N - 1000, 1000) << endl;
                    promoter.push_back(genomeSequence.substr(N - 1000, 1000));
                    break; 
                }
            }
        }

        // Update searchStart to continue searching from the next position after the current TATA sequence
        searchStart = N + 1001;
        if (searchStart >= len)
            break;
    }
}