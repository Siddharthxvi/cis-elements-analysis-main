#include "impl.h"
using namespace std;

int main()
{
    string A, B;
    cout << "Enter first Cis-element: ";
    cin >> A;
    cout << "Enter second Cis-element: ";
    cin >> B;
    // A = "ACGT";
    // B= "AAAG";
    cout << "Enter Number of chromosomes in the genome: ";
    int chromosome_num;
    cin >> chromosome_num;

    cout << "Enter name of Folder containing chromosome data: ";
    string foldername;
    cin >> foldername;
    //foldername = "SweetPotato";

    ofstream outputFile("promoter.txt");
    vector<string> promoter;                                      // TO STORE ALL THE PROMOTERS FOUND
    vector<string> genomeSpacers, genomeSpacersReverse;           // Normal and Reverse Orientation
    vector<string> SpacersPerChromosome, SpacersPerChromosomeRev; // To store Spacers per chromosome;

    // -----------------------------
    vector<int> TotalSpacerFreqInGenome(21, 0);
    vector<int> TotalSpacerFreqInPromoter(21, 0);
    // ---------------------------------

    string genomeSequence;
    long int TotalSpacersInGnme = 0, TotalSpacersInPrmtr = 0;
    long long TotalNumPromoters = 0;

    int Start = 1;
    printf("Enter the starting Chromosome:\t");
    cin >> Start;

    // Combined file containing spacer sequences of all chromosomes
    ofstream CombinedSpacersInPromoterFile("SpacersInPromoter_all_chromosomes.fa", ios::out);

    vector<string> chrom_id_genome;
    vector<int> spacer_start_genome;
    vector<int> spacer_end_genome;

    vector<string> chrom_id_promoter;
    vector<int> spacer_start_promoter;
    vector<int> spacer_end_promoter;

    for (int i = Start; i <= chromosome_num; ++i)
    {

        string filename = "/Users/parikshithmoleyar/Documents/BIO/SOP/Cis-element-analysis-main/Data/" + foldername + "/chr" + to_string(i) + ".fasta";
        ifstream inputFile(filename);
        if (!inputFile)
        {
            cerr << "Failed to open file: " << filename << endl;
            continue;
        }

        string ChromosomeSeq, line, partialSequence;

        int CountForPartialSeq = 0;
        int CountGenomeSpacers = 0;

        string chrom_id; // Ex. NC_021160.1
        long long genomeOffset = 0;  // tracks position in the whole chromosome

        while (getline(inputFile, line))
        {
            if (!line.empty() && line[0] != '>')
            {
                ChromosomeSeq += line;
                partialSequence += line;
                CountForPartialSeq++;

                genomeOffset += line.length();   // increment running position
            }
            else if (!line.empty() && line[0] == '>') {
                size_t first = line.find("|");
                size_t second = line.find("|", first + 1);
                chrom_id = line.substr(first + 1, second - (first + 1));
            }
            if (CountForPartialSeq == 10000) // Partial Sequence is Now 10000 lines
            {
                findSpacers(partialSequence, SpacersPerChromosome, CountGenomeSpacers, A, B, chrom_id_genome, spacer_start_genome, spacer_end_genome, chrom_id, genomeOffset - partialSequence.length());

                CountForPartialSeq = 0;
                partialSequence.clear();
            }
        }
        // TO TAKE CARE OF REMAINING LINES
        if (CountForPartialSeq > 0)
        {
            findSpacers(partialSequence, SpacersPerChromosome, CountGenomeSpacers, A, B, chrom_id_genome, spacer_start_genome, spacer_end_genome, chrom_id, genomeOffset - partialSequence.length());
            CountForPartialSeq = 0;
            partialSequence.clear();
        }
        // ---- I HAVE SPACERS FOR iTH CHROMOSOME ------
        
        //Find conserved spacer sequences of genome
        double conservationThreshold = 0.7; // threshold for nucleotide to be considered conserved 
        vector<string> ConservedSequencesInGenome = findConservedSequences(SpacersPerChromosome, conservationThreshold);

        // --- MAINTAINING A FREQUENCY OF SPACERS AND GC COUNT FOR N .. 0-20
        vector<int> ChromosomeSpacerFreq(21, 0);
        int ChromosomeGCCount;
        int SpacerNucleotideCount; //count of total nucleotides in all spacers
        float ChromosomeGCContent;
        for (const auto &spacer : SpacersPerChromosome)
        {
            ChromosomeSpacerFreq[spacer.size()]++;

            for (const auto &sp : spacer){
                SpacerNucleotideCount++;
                if (sp == 'G' || sp == 'C') ChromosomeGCCount++;
            }
        }

        for (int i = 0; i <= 20; ++i)
        {
            TotalSpacerFreqInGenome[i] += ChromosomeSpacerFreq[i];
        }

        // --------------------------
        // To input fre per chromosome in the file
        string S1 = "ChromosomeSpacerFreq" + to_string(i) + ".txt";
        ofstream ChromosomeSpacerFrqFile(S1);
        for (int k = 0; k <= 20; ++k)
        {
            ChromosomeSpacerFrqFile << ChromosomeSpacerFreq[k] << endl;
        }

        // ---------------------------

        //To input GC Content per chromosome in the file
        string S2 = "ChromosomeGCContent(%)" + to_string(i) + ".txt";
        ofstream ChromosomeGCContentFile(S2);
        ChromosomeGCContent = ((ChromosomeGCCount)/(SpacerNucleotideCount)) * 100;
        ChromosomeGCContentFile << ChromosomeGCContent << endl; // % GC content

        // ---------------------------

        //To input conserved sequence per chromosome in the file
        string S3 = "ChromosomeConsensus" + to_string(i) + ".txt";
        ofstream ChromosomeConsensusFile(S3);
        for (int k = 0; k <= 20; ++k)
        {
            ChromosomeConsensusFile << ConservedSequencesInGenome[k] << endl;
        }

        // ---------------------------

        // ---------------

        vector<string> PromoterPerChr; // TO STORE PROMTOTERS FOUND

        // --------- FIND PROMOTERS PER CHROMOSOME -----
        findPromoter(ChromosomeSeq, PromoterPerChr);
        TotalNumPromoters += PromoterPerChr.size();
        printf("Size of one Promoter Seq : %zu\n", PromoterPerChr[1].size());

        // -------- FIND SPACERS FOR ALL PROMOTERS------
        int Count = 0;
        vector<string> SpacersInPromoter;
        for (const auto &it : PromoterPerChr)
        {
            // cout << it << "\n";
            findSpacers(it, SpacersInPromoter, Count, A, B, chrom_id_promoter, spacer_start_promoter, spacer_end_promoter, chrom_id, 0);
        }

        //Find conserved spacer sequences of promoter
        vector<string> ConservedSequencesInPromoter = findConservedSequences(SpacersInPromoter, conservationThreshold);

        // -----------  GET THE FREQ OF THE SPACERS AND GC COUNT---
        vector<int> ChromosomeSpacerFreqInPromoter(21, 0);
        int ChromosomeGCCountInPromoter;
        int SpacerNucleotideCountInPromoter;
        float ChromosomeGCContentInPromoter; // % GC content
        for (const auto &spacer : SpacersInPromoter)
        {
            // cout << spacer << "\n";
            ChromosomeSpacerFreqInPromoter[spacer.size()]++;

            for (const auto &sp : spacer){
                SpacerNucleotideCountInPromoter++;
                if (sp == 'G' || sp == 'C') ChromosomeGCCountInPromoter++;
            }
        }

        for (int i = 0; i <= 20; ++i)
        {
            TotalSpacerFreqInPromoter[i] += ChromosomeSpacerFreqInPromoter[i];
        }

        // ----------- INPUT THE FREQ IN TO FILE ----
        string St1 = "ChromosomeSpacerFreqPromoter" + to_string(i) + ".txt";
        ofstream ChromosomeSpacerFrqInPromoterFile(St1);
        for (int k = 0; k <= 20; ++k)
        {
            ChromosomeSpacerFrqInPromoterFile << ChromosomeSpacerFreqInPromoter[k] << endl;
        }
        // ---------------------------

        // ----------- INPUT THE GC CONTENT % IN TO FILE ----
        string St2 = "ChromosomeGCContent(%)Promoter" + to_string(i) + ".txt";
        ofstream ChromosomeGCContentInPromoterFile(St2);
        ChromosomeGCContentInPromoter = (ChromosomeGCCountInPromoter * 100)/(SpacerNucleotideCountInPromoter);
        ChromosomeGCContentInPromoterFile << ChromosomeGCContentInPromoter << endl;
        // ---------------------------

        // ----------- INPUT THE CONSERVED SEQUENCE IN TO FILE ----
        string St3 = "ChromosomeConsensusPromoter" + to_string(i) + ".txt";
        ofstream ChromosomeConsensusInPromoterFile(St3);
        for (int k = 0; k <= 20; ++k)
        {
            ChromosomeConsensusInPromoterFile << ConservedSequencesInPromoter[k] << endl;
        }
        // ---------------------------

        // ----------- INPUT THE SPACER SEQUENCES IN TO FILE (spacers of length > 6) ----
        // string St4 = "SpacersInPromoter_chr" + to_string(i) + ".fa";
        // ofstream SpacersInPromoterFile(St4);
        int spacer_id = 1; // To number the spacers
        for (const auto &spacer : SpacersInPromoter){
            if (spacer.length() < 6) continue;
            // SpacersInPromoterFile << ">" << foldername << "_chr" << i << "_seq" << seq_id << "\n";
            // SpacersInPromoterFile << spacer << "\n";
            // CombinedSpacersFile << ">" << foldername << "(" << A << "-" << B << ")_chr" << i << "_seq" << seq_id << "\n";
            CombinedSpacersInPromoterFile << ">" << chrom_id_promoter[spacer_id - 1] << ":" << spacer_start_promoter[spacer_id - 1] << "-" << spacer_end_promoter[spacer_id - 1] << "|(" << A << "-" << B << ")|seq" << spacer_id << "\n";            CombinedSpacersInPromoterFile << spacer << "\n";
            spacer_id++;
        }
        // ---------------------------

        // FREE MEM FOR SPACERS AND PROMOTERS
        PromoterPerChr.clear();
        SpacersInPromoter.clear();
        SpacersPerChromosome.clear();
        ChromosomeSeq.clear();

        // CLOSE THE DATA FILE
        inputFile.close();

        // -------------------------------------------------------------------

        // printf("Number of Genome Spacers in chromosome %d: %ld\n", i, genomeSpacers.size());
        // printf("Number of Genome Spacers in chromosome %d in Reverse-Orientation : %ld\n", i, genomeSpacersReverse.size());
        // cout << "Number of bases in chromosome " << i << ": " << genomeSequence.size() << " Nucleotides" << endl;
    }
    for (auto &count : TotalSpacerFreqInGenome)
        TotalSpacersInGnme += count;
    for (auto &count : TotalSpacerFreqInPromoter)
        TotalSpacersInPrmtr += count;

    CombinedSpacersInPromoterFile.close();

    printf("Total Number of Spacers in Entire Genome: %ld\n", TotalSpacersInGnme);
    printf("Total Number of Spacers in  Promoters: %ld\n", TotalSpacersInPrmtr);

    string Str = "TotalSpacerFreqinGenome.txt";
    ofstream TotalSpacerFreqinGenomeFile(Str);
    for (int i = 0; i <= 20; ++i)
    {
        TotalSpacerFreqinGenomeFile<< TotalSpacerFreqInGenome[i] << endl;
    }

    string PromSpacers = "TotalSpacersInPromoters.txt";
    ofstream SpacerfrqInPromters(PromSpacers);
    for (int i = 0; i <= 20; ++i)
    {
        SpacerfrqInPromters<< TotalSpacerFreqInPromoter[i] << endl;
    }

    printf("TOTAL NUMBER OF PROMOTERS FOUND: %lld\n", TotalNumPromoters);

    return 0;
}
