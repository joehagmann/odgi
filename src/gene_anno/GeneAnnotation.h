#include <algorithm>
#include <atomic>
#include <cstring>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <exception>
#include "Gene.h"
#include "Interval.h"

#ifndef GENEANNOTATION_H
#define GENEANNOTATION_H

class GeneBin {
using ull_int = unsigned long long;
public:
    // class data
    ull_int start;
    ull_int end = 0;
    std::vector<Gene> genes;

    // add gene to bin
    void add_gene(Gene gene)
    {
        if (genes.size() == 0) {
            start = gene.st;
        }
        genes.push_back(gene);
        // extend boundaries of bin to fully include bin
        if (gene.st < start)
        {
            start = gene.st;
        }
        if (gene.en > end)
        {
            end = gene.en;
        }
    }

    // check if interval overlaps bin
    const bool overlaps(const Interval &it) {
        return !(start > it.en) && !(end < it.st);
    }
};

class GeneBins {
public:
    // class data
    std::vector<GeneBin> gene_bins;

    // get vector of bins overlapping query interval
    std::vector<GeneBin*> get_bins(Interval it)
    {
        std::vector<GeneBin*> overlapped_bins;
        for (auto &bin : gene_bins) {
            if (bin.overlaps(it)) {
                overlapped_bins.push_back(&bin);
            }
        }
        return overlapped_bins;
    }

    // create bins from vector of genes
    void make_bins(std::vector<Gene> &genes) {
        unsigned int bin_index = 0;
        unsigned int count = 0;

        for (auto gene : genes) {
            if (bin_index + 1 > gene_bins.size()) {
                // add bin if required
                gene_bins.resize(gene_bins.size() + 1);
            }

            gene_bins[bin_index].add_gene(gene);
            count++;

            if (count == bin_size) {
                // bin is full, move to next bin
                count = 0;
                bin_index++;
            }
        }
    }
private:
    const unsigned int bin_size = 64;
};

// parse gff3 genome annotation
class GeneAnnotation
{
public:
    std::string anno_source = "";
    std::unordered_set<std::string> recorded_genes;

    std::unordered_map<std::string, std::vector<Gene>> gene_dict;
    std::unordered_map<std::string, GeneBins> bins_dict;

    //get number of genes
    int ngenes();

    //return all gene id as a std::vector
    std::vector<std::string> get_genelist();

    void parse_gff3_annotation(std::string gff3_fn, bool fix_chrname);

    friend std::ostream& operator<< (std::ostream& out, const GeneAnnotation& obj);

private:
    // index variables for gff3 fields
    const int SEQID      = 0;
    const int SOURCE     = 1;
    const int TYPE       = 2;
    const int START      = 3;
    const int END        = 4;
    const int SCORE      = 5;
    const int STRAND     = 6;
    const int PHASE      = 7;
    const int ATTRIBUTES = 8;

    // get attribute from gff3 standard columns
    std::string get_attribute(const std::vector<std::string> &all_attributes, const std::string &target_attribute);
    // convert strand from +- symbols to -1 or 1
    int get_strand(char st);

    const std::string get_parent(const std::vector<std::string> &attributes);
    std::string get_ID(const std::vector<std::string> &attributes);

    // add chr to molecule names if requested
    std::string fix_name(std::string chr_name);

    // parse entry of gff3 annotation and add its information to this object
    void parse_anno_entry(
        const bool &fix_chrname,
        const std::string &line,
        std::unordered_map<std::string, std::unordered_map<std::string, Gene>> &chr_to_genes_dict,
        std::unordered_map<std::string, std::string> &transcript_to_gene_dict
    );

    // generic gene_id getter for gff3 entries
    std::string get_gene_id(const std::vector<std::string> &attributes);

    // specific gene_id getter for gff3 entries
    std::string get_gencode_gene_id(const std::vector<std::string> &attributes);
    std::string get_refseq_gene_id(const std::vector<std::string> &attributes);
    //std::string get_plain_gene_id(const std::vector<std::string> &attributes);

    // guess the source of annotation
    std::string guess_anno_source(std::string gff3_fn);

    const bool parent_is_gene(const std::string &parent);
    const bool parent_is_known_transcript(const std::unordered_map<std::string, std::string> &transcript_to_gene_dict, const std::string &parent);
    const bool is_gene(const std::vector<std::string> &fields, const std::vector<std::string> &attributes);
    const bool is_exon(const std::vector<std::string> &fields, const std::vector<std::string> &attributes);
    const bool is_transcript(const std::vector<std::string> &fields, const std::vector<std::string> &attributes);
};

#endif