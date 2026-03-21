#include "impl.h"
using namespace std;

void findSpacers(const string &PromoterSequence,
                 vector<string> &spacers,
                 int &count,
                 string a, string b,
                 vector<string> &chrom_id,
                 vector<int> &spacer_start,
                 vector<int> &spacer_end,
                 const string &currentChrom,
                 int offset)
{
    size_t n = PromoterSequence.size();
    size_t pos = 0;

    while (pos < n)
    {
        size_t l = PromoterSequence.find(a, pos);
        if (l == string::npos)
            break;

        size_t searchStart = l + a.size();
        if (searchStart >= n)
            break;

        size_t h = PromoterSequence.find(b, searchStart);
        if (h == string::npos)
            break;

        size_t SpaceLen = h - searchStart;

        if (SpaceLen < 21)
        {
            count++;

            spacers.push_back(PromoterSequence.substr(searchStart, SpaceLen));
            chrom_id.push_back(currentChrom);
            spacer_start.push_back(offset + searchStart);
            spacer_end.push_back(offset + h - 1);
        }

        pos = h + b.size();
        if (pos <= l) pos = l + 1;
    }
}