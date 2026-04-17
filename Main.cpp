#include "impl.h"
using namespace std;

int main()
{
    string A, B;
    cout << "Enter first Cis-element: ";
    cin >> A;
    cout << "Enter second Cis-element: ";
    cin >> B;

    for (char &c : A) c = toupper(c);
    for (char &c : B) c = toupper(c);

    cout << "Enter Number of chromosomes in the genome: ";
    int chromosome_num;
    cin >> chromosome_num;

    cout << "Enter name of Folder containing chromosome data: ";
    string foldername;
    cin >> foldername;

    int Start = 1;
    cout << "Enter the starting Chromosome:\t";
    cin >> Start;

    ofstream CombinedSpacersInPromoterFile("SpacersInPromoter_all_chromosomes.fa");

    ofstream gcGenomeAll("ChromosomeGCContentGenome_all.txt");
    ofstream gcPromoterAll("ChromosomeGCContentPromoter_all.txt");

    vector<int> TotalSpacerFreqInGenome(31, 0);
    vector<int> TotalSpacerFreqInPromoter(31, 0);

    long long TotalNumPromoters = 0;

    vector<string> chrom_id_genome;
    vector<int> spacer_start_genome;
    vector<int> spacer_end_genome;

    for (int i = Start; i <= chromosome_num; ++i)
    {
        string filename = "/Users/siddharth/Desktop/sid/3-2/sop/Cis-element-analysis-main/Data/" + foldername + "/chr" + to_string(i) + ".fasta";

        ifstream inputFile(filename);
        if (!inputFile)
        {
            cerr << "Failed to open file: " << filename << endl;
            continue;
        }

        string ChromosomeSeq, line, partialSequence;
        string chrom_id;

        int CountForPartialSeq = 0;
        int CountGenomeSpacers = 0;
        long long genomeOffset = 0;

        vector<string> SpacersPerChromosome;

        while (getline(inputFile, line))
        {
            if (!line.empty() && line[0] == '>')
            {
                string header = line.substr(1);
                size_t spacePos = header.find(' ');
                chrom_id = (spacePos != string::npos) ? header.substr(0, spacePos) : header;
            }
            else if (!line.empty())
            {
                for (char &c : line) c = toupper(c);

                ChromosomeSeq += line;
                partialSequence += line;
                genomeOffset += line.length();
                CountForPartialSeq++;
            }

            if (CountForPartialSeq == 10000)
            {
                findSpacers(partialSequence, SpacersPerChromosome, CountGenomeSpacers,
                            A, B, chrom_id_genome, spacer_start_genome, spacer_end_genome,
                            chrom_id, genomeOffset - partialSequence.length());

                partialSequence.clear();
                partialSequence.shrink_to_fit();
                CountForPartialSeq = 0;
            }
        }

        if (!partialSequence.empty())
        {
            findSpacers(partialSequence, SpacersPerChromosome, CountGenomeSpacers,
                        A, B, chrom_id_genome, spacer_start_genome, spacer_end_genome,
                        chrom_id, genomeOffset - partialSequence.length());
        }

        // ================= GENOME SPACER GC (FIXED) =================
        long long gcGenome = 0, totalGenome = 0;

        vector<int> ChromosomeSpacerFreq(31, 0);
        for (const auto &spacer : SpacersPerChromosome)
        {
            if (spacer.size() <= 30)
                ChromosomeSpacerFreq[spacer.size()]++;

            for (char c : spacer)
            {
                c = toupper(c);

                if (c == 'G' || c == 'C')
                {
                    gcGenome++;
                    totalGenome++;
                }
                else if (c == 'A' || c == 'T')
                {
                    totalGenome++;
                }
            }
        }

        for (int k = 0; k <= 30; ++k)
            TotalSpacerFreqInGenome[k] += ChromosomeSpacerFreq[k];

        float genomeGC = (totalGenome > 0)
            ? (100.0f * gcGenome / totalGenome)
            : 0.0f;

        gcGenomeAll << genomeGC << endl;
        // ===========================================================

        // ================= PROMOTERS =================
        vector<string> PromoterPerChr;
        findPromoter(ChromosomeSeq, PromoterPerChr);
        TotalNumPromoters += PromoterPerChr.size();

        ChromosomeSeq.clear();
        ChromosomeSeq.shrink_to_fit();

        vector<string> chrom_id_promoter;
        vector<int> spacer_start_promoter;
        vector<int> spacer_end_promoter;

        vector<string> SpacersInPromoter;
        int Count = 0;

        for (const auto &p : PromoterPerChr)
        {
            findSpacers(p, SpacersInPromoter, Count,
                        A, B,
                        chrom_id_promoter, spacer_start_promoter, spacer_end_promoter,
                        chrom_id, 0);
        }

        cout << "Chr " << i
             << " Promoters: " << PromoterPerChr.size()
             << " Spacers: " << SpacersInPromoter.size() << endl;

        // ================= CONSERVATION =================
        double conservationThreshold = 0.5;
        vector<string> ConservedSequencesInPromoter =
            findConservedSequences(SpacersInPromoter, conservationThreshold);

        string consPromFile = "ChromosomeConsensusPromoter_" + to_string(i) + ".txt";
        ofstream foutConsProm(consPromFile);

        for (int k = 0; k <= 30; ++k)
        {
            string seq = (k < ConservedSequencesInPromoter.size())
                         ? ConservedSequencesInPromoter[k]
                         : "";
            foutConsProm << seq << endl;
        }
        // =================================================

        // ================= PROMOTER SPACER GC (FIXED) =================
        long long gcProm = 0, totalProm = 0;

        vector<int> ChromosomeSpacerFreqInPromoter(31, 0);
        for (const auto &spacer : SpacersInPromoter)
        {
            if (spacer.size() <= 30)
                ChromosomeSpacerFreqInPromoter[spacer.size()]++;

            for (char c : spacer)
            {
                c = toupper(c);

                if (c == 'G' || c == 'C')
                {
                    gcProm++;
                    totalProm++;
                }
                else if (c == 'A' || c == 'T')
                {
                    totalProm++;
                }
            }
        }

        for (int k = 0; k <= 30; ++k)
            TotalSpacerFreqInPromoter[k] += ChromosomeSpacerFreqInPromoter[k];

        float promoterGC = (totalProm > 0)
            ? (100.0f * gcProm / totalProm)
            : 0.0f;

        gcPromoterAll << promoterGC << endl;
        // =============================================================

        int spacer_id = 1;
        for (const auto &spacer : SpacersInPromoter)
        {
            if (spacer.length() < 6)
                continue;

            int idx = spacer_id - 1;

            if (idx < chrom_id_promoter.size())
            {
                CombinedSpacersInPromoterFile << ">"
                                             << chrom_id_promoter[idx] << ":"
                                             << spacer_start_promoter[idx] << "-"
                                             << spacer_end_promoter[idx]
                                             << "|(" << A << "-" << B << ")|seq" << spacer_id << "\n";

                CombinedSpacersInPromoterFile << spacer << "\n";
            }

            spacer_id++;
        }

        SpacersPerChromosome.clear(); SpacersPerChromosome.shrink_to_fit();
        PromoterPerChr.clear(); PromoterPerChr.shrink_to_fit();
        SpacersInPromoter.clear(); SpacersInPromoter.shrink_to_fit();

        chrom_id_genome.clear();
        spacer_start_genome.clear();
        spacer_end_genome.clear();

        inputFile.close();
    }

    ofstream genomeTotalFile("TotalSpacerFreqGenome.txt");
    for (int k = 0; k <= 30; ++k)
        genomeTotalFile << TotalSpacerFreqInGenome[k] << endl;

    ofstream promoterTotalFile("TotalSpacerFreqPromoter.txt");
    for (int k = 0; k <= 30; ++k)
        promoterTotalFile << TotalSpacerFreqInPromoter[k] << endl;

    CombinedSpacersInPromoterFile.close();

    cout << "Total Promoters: " << TotalNumPromoters << endl;

    return 0;
}