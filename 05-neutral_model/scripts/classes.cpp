#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>

// Define the codon to aminoacid map. For convinience Z will be the Stop signal
std::map<std::string, char> codon_2_seq = {
    {"UUU", 'F'},
    {"UUC", 'F'},
    {"UUA", 'L'},
    {"UUG", 'L'},
    {"UCU", 'S'},
    {"UCC", 'S'},
    {"UCA", 'S'},
    {"UCG", 'S'},
    {"UAU", 'Y'},
    {"UAC", 'Y'},
    {"UAA", 'Z'},
    {"UAG", 'Z'},
    {"UGU", 'C'},
    {"UGC", 'C'},
    {"UGA", 'Z'},
    {"UGG", 'W'},
    {"CUU", 'L'},
    {"CUC", 'L'},
    {"CUA", 'L'},
    {"CUG", 'L'},
    {"CCU", 'P'},
    {"CCC", 'P'},
    {"CCA", 'P'},
    {"CCG", 'P'},
    {"CAU", 'H'},
    {"CAC", 'H'},
    {"CAA", 'Q'},
    {"CAG", 'Q'},
    {"CGU", 'R'},
    {"CGC", 'R'},
    {"CGA", 'R'},
    {"CGG", 'R'},
    {"AUU", 'I'},
    {"AUC", 'I'},
    {"AUA", 'I'},
    {"AUG", 'M'},
    {"ACU", 'T'},
    {"ACC", 'T'},
    {"ACA", 'T'},
    {"ACG", 'T'},
    {"AAU", 'N'},
    {"AAC", 'N'},
    {"AAA", 'K'},
    {"AAG", 'K'},
    {"AGU", 'S'},
    {"AGC", 'S'},
    {"AGA", 'R'},
    {"AGG", 'R'},
    {"GUU", 'V'},
    {"GUC", 'V'},
    {"GUA", 'V'},
    {"GUG", 'V'},
    {"GCU", 'A'},
    {"GCC", 'A'},
    {"GCA", 'A'},
    {"GCG", 'A'},
    {"GAU", 'D'},
    {"GAC", 'D'},
    {"GAA", 'E'},
    {"GAG", 'E'},
    {"GGU", 'G'},
    {"GGC", 'G'},
    {"GGA", 'G'},
    {"GGG", 'G'}};

std::vector<std::string> seq_2_protein(std::string seq, int initial = 0, bool export_full = false)
{
    // Function that given a sequence, returns the protein(s) sequence
    // If export full is false, then we just keep the initial part of the protein (untill a STOP codon is found)

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper); // Upper case transform
    std::replace(seq.begin(), seq.end(), 'T', 'U');                 // Change the timine by uracil
    std::vector<char> aminos;                                       // Create the container for the proteins

    for (int i = initial; i < seq.length() - 2; i += 3)
    {
        std::string codon = seq.substr(i, 3); // Select the codon
        auto iterator = codon_2_seq.find(codon);
        if (iterator != codon_2_seq.end())
        {
            char aminoacid = iterator->second;
            aminos.push_back(aminoacid);
        }
        else
        {
            std::cout << "Codon " << codon << " not in dictionary. We are in the position i=" << i << std::endl;
        }
    }

    std::string protein;
    for (char amin : aminos)
    {
        protein += amin;
    }

    std::string delimiter = "Z";
    std::vector<std::string> splitParts;

    if (export_full == true)
    {
        size_t pos = 0;
        size_t found;

        while ((found = protein.find(delimiter, pos)) != std::string::npos)
        {
            std::string part = protein.substr(pos, found - pos);
            splitParts.push_back(part);
            pos = found + delimiter.length();
        }

        std::string lastPart = protein.substr(pos);
        splitParts.push_back(lastPart);
    }
    else
    {
        size_t found = protein.find(delimiter);
        if (found != std::string::npos)
        {
            splitParts.push_back(protein.substr(0, found));
        }
    }

    return splitParts;
}

std::map<int, std::vector<std::string>> compare_proteins(std::vector<std::string> amin_1, std::vector<std::string> amin_2)
{
    std::map<int, std::vector<std::string>> results;

    if (amin_1.size() != amin_2.size())
    {
        std::cout << "Amino acids of different lengths!" << std::endl;
    }
    else
    {
        for (int i = 0; i < amin_1.size(); i++)
        {
            if (amin_1[i] != amin_2[i])
            {
                std::vector<std::string> diff = {amin_1[i], amin_2[i]};
                results[i] = diff;
            }
        }
    }

    return results;
}

int main()
{
    std::string our_reference_sequence_qbeta = "CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCGTAAAGCGTTGAAACTTTG";

    std::vector<std::string> protein_wt = seq_2_protein(our_reference_sequence_qbeta);

    // Print the protein
    std::cout << "Our wild type produces the protein ";
    for (const std::string part : protein_wt)
    {
        std::cout << part << std::endl;
    }

    // Now we need to check the protein comparator
    //******************************************

    // std::vector<std::string> aminos_1 = seq_2_amino("CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATACGAGCAATCTTCCGTTCGCTACACTAGATCTTGATACTACCTTTAGTTCGTTTAAACACGTTCTTGATAGTATCTTTTTATTAACCCAACGCGTAAAGCGTTGAAACTTTG");
    // std::vector<std::string> aminos_2 = seq_2_amino("CAACAAGGTCAGCTATATCATAATATCGATATTGTAGACGGCTTTGACAGACGTGACATCCGGCTCAAATCTTTCACCATAAAAGGTGAACGAAATGGGCGGCCTGTTAACGTTTCTGCTAGCCTGTCTGCTGTCGATTTATTTTACAGCCGACTCCATAC
    return 0;
}