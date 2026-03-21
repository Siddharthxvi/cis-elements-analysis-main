#include "impl.h"
using namespace std;

void findPromoter(const string &genomeSequence, vector<string> &promoter)
{
    size_t len = genomeSequence.size();
    size_t searchStart = 0;

    while (true)
    {
        size_t N = genomeSequence.find("TATA", searchStart);
        if (N == string::npos)
            break;

        if (N > 1000)
        {
            for (size_t i = N + 28; i <= N + 35 && i + 3 <= len; ++i)
            {
                if (genomeSequence.substr(i, 3) == "ATG")
                {
                    promoter.push_back(genomeSequence.substr(N - 1000, 1000));
                    break;
                }
            }
        }

        searchStart = N + 1001;
        if (searchStart >= len)
            break;
    }
}