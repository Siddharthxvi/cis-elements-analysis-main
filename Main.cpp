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

    vector<int> TotalSpacerFreqInGenome(21, 0);
    vector<int> TotalSpacerFreqInPromoter(21, 0);

    long int TotalSpacersInGnme = 0, TotalSpacersInPrmtr = 0;
    long long TotalNumPromoters = 0;

    long long totalGenomeGC = 0, totalGenomeBases = 0;
    long long totalPromoterGC = 0, totalPromoterBases = 0;

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
                CountForPartialSeq = 0;
            }
        }

        if (!partialSequence.empty())
        {
            findSpacers(partialSequence, SpacersPerChromosome, CountGenomeSpacers,
                        A, B, chrom_id_genome, spacer_start_genome, spacer_end_genome,
                        chrom_id, genomeOffset - partialSequence.length());
        }

        // GENOME FREQUENCY + GC
        vector<int> ChromosomeSpacerFreq(21, 0);
        int gcCountGenome = 0, baseCountGenome = 0;

        for (const auto &spacer : SpacersPerChromosome)
        {
            if (spacer.size() <= 20)
                ChromosomeSpacerFreq[spacer.size()]++;

            for (char c : spacer)
            {
                baseCountGenome++;
                if (c == 'G' || c == 'C') gcCountGenome++;
            }
        }

        for (int k = 0; k <= 20; ++k)
            TotalSpacerFreqInGenome[k] += ChromosomeSpacerFreq[k];

        string genomeFile = "ChromosomeSpacerFreqGenome_" + to_string(i) + ".txt";
        ofstream fout1(genomeFile);
        for (int k = 0; k <= 20; ++k)
            fout1 << ChromosomeSpacerFreq[k] << endl;

        float gcGenomePercent = (baseCountGenome > 0) ? (100.0 * gcCountGenome / baseCountGenome) : 0;
        gcGenomeAll << gcGenomePercent << endl;

        totalGenomeGC += gcCountGenome;
        totalGenomeBases += baseCountGenome;

        // PROMOTERS
        vector<string> PromoterPerChr;
        findPromoter(ChromosomeSeq, PromoterPerChr);

        TotalNumPromoters += PromoterPerChr.size();

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

        // ALTER THRESHOLD HERE
        double conservationThreshold = 0.5;
        vector<string> ConservedSequencesInPromoter =
            findConservedSequences(SpacersInPromoter, conservationThreshold);

        string consPromFile = "ChromosomeConsensusPromoter_" + to_string(i) + ".txt";
        ofstream foutConsProm(consPromFile);
        for (int k = 0; k <= 20; ++k)
            foutConsProm << ConservedSequencesInPromoter[k] << endl;

        // PROMOTER FREQUENCY + GC
        vector<int> ChromosomeSpacerFreqInPromoter(21, 0);
        int gcCountProm = 0, baseCountProm = 0;

        for (const auto &spacer : SpacersInPromoter)
        {
            if (spacer.size() <= 20)
                ChromosomeSpacerFreqInPromoter[spacer.size()]++;

            for (char c : spacer)
            {
                baseCountProm++;
                if (c == 'G' || c == 'C') gcCountProm++;
            }
        }

        for (int k = 0; k <= 20; ++k)
            TotalSpacerFreqInPromoter[k] += ChromosomeSpacerFreqInPromoter[k];

        string promoterFile = "ChromosomeSpacerFreqPromoter_" + to_string(i) + ".txt";
        ofstream fout2(promoterFile);
        for (int k = 0; k <= 20; ++k)
            fout2 << ChromosomeSpacerFreqInPromoter[k] << endl;

        float gcPromPercent = (baseCountProm > 0) ? (100.0 * gcCountProm / baseCountProm) : 0;
        gcPromoterAll << gcPromPercent << endl;

        totalPromoterGC += gcCountProm;
        totalPromoterBases += baseCountProm;

        // WRITE SPACERS
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

        ChromosomeSeq.clear();
        SpacersPerChromosome.clear();
        PromoterPerChr.clear();
        SpacersInPromoter.clear();

        inputFile.close();
    }

    for (auto &count : TotalSpacerFreqInGenome)
        TotalSpacersInGnme += count;

    for (auto &count : TotalSpacerFreqInPromoter)
        TotalSpacersInPrmtr += count;

    CombinedSpacersInPromoterFile.close();

    cout << "Total Genome Spacers: " << TotalSpacersInGnme << endl;
    cout << "Total Promoter Spacers: " << TotalSpacersInPrmtr << endl;
    cout << "Total Promoters: " << TotalNumPromoters << endl;

    float finalGenomeGC = (totalGenomeBases > 0) ? (100.0 * totalGenomeGC / totalGenomeBases) : 0;
    float finalPromoterGC = (totalPromoterBases > 0) ? (100.0 * totalPromoterGC / totalPromoterBases) : 0;

    cout << "Final Genome GC%: " << finalGenomeGC << endl;
    cout << "Final Promoter GC%: " << finalPromoterGC << endl;

    ofstream genomeTotalFile("TotalSpacerFreqGenome.txt");
    for (int k = 0; k <= 20; ++k)
        genomeTotalFile << TotalSpacerFreqInGenome[k] << endl;

    ofstream promoterTotalFile("TotalSpacerFreqPromoter.txt");
    for (int k = 0; k <= 20; ++k)
        promoterTotalFile << TotalSpacerFreqInPromoter[k] << endl;

    return 0;
}