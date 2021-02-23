#include "GeneAnnotation.hpp"

using std::atoi;
using std::atomic;
using std::cout;
using std::endl;
using std::fixed;
using std::getline;
using std::ifstream;
using std::ostream;
using std::sort;
using std::string;
using std::stringstream;
using std::thread;
using std::unordered_map;
using std::map;
using std::vector;

namespace odgi {
    namespace gene_anno {

// helper functions
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


//////////////////////////////////


string GeneAnnotation::get_attribute(
    const vector<string> &all_attributes,
    const string &target_attribute
)
{
    for (const string &attr : all_attributes) {
        auto sep_loc = attr.find("=");
        // get key
        const string key = attr.substr(0, sep_loc);
        // get value
        const string val = attr.substr(sep_loc + 1);
        if (key == target_attribute) {
            return val;
        }
    }
    return "";
}

int GeneAnnotation::get_strand(char st)
{
    int strand = 0;
    if (st == '+')
    {
        strand = 1;
    }
    else if (st == '-')
    {
        strand = -1;
    }
    return strand;
}

string GeneAnnotation::get_ID(const vector<string> &attributes)
{
    for (const auto &attr : attributes)
    {
        if (attr.substr(0, 2) == "ID")
        {
            // check for ENSEMBL notation
            if (anno_source == "ensembl")
            {
                return attr.substr(attr.rfind(':') + 1);
            }
            else
            {
                return attr.substr(attr.find('=') + 1);
            }
        }
        else if (attr.substr(0, 6) == "Parent")
        {
            return get_parent(attributes);
        }
    }
    return "";
}

const string GeneAnnotation::get_parent(const vector<string> &attributes)
{
    for (const auto &attr : attributes)
    {
        if (attr.substr(0, 6) == "Parent")
        {
            // clip away isoform
            //string parent = attr.substr(0, attr.find('.'));
            // check for ENSEMBL notation
            if (anno_source == "ensembl")
            {
                return attr.substr(attr.rfind(':') + 1);
                //return parent.substr(parent.rfind(':') + 1);
            }
            else
            {
                return attr.substr(attr.find('=') + 1);
                //return parent.substr(parent.find('=') + 1);
            }
        }
    }
    return "";
}

const string GeneAnnotation::get_parent_of_CDS(
    const vector<string> &attributes,
    const unordered_map<string, string> &transcript_to_gene_dict)
{
    for (const auto &attr : attributes)
    {
        if (anno_source == "ensembl" && attr.substr(0, 7) == "Parent:")
        {
            return attr.substr(attr.rfind(':') + 1);
        } 
        else if (anno_source != "ensembl" && attr.substr(0, 7) == "Parent=")
        {
            std::string s = attr.substr(attr.find('=') + 1);

            // there can be multiple IDs separated by comma. Check them
            // until you find a transcript:
            std::string parent_transcript = s;
            auto start = 0;
            auto end = s.find(',');
            while (end != std::string::npos)
            {
                parent_transcript = s.substr(start, end - start);
                if (transcript_to_gene_dict.find(parent_transcript) != transcript_to_gene_dict.end())
                {
                    return parent_transcript;
                }
                start = end + 1;
                end = s.find(',', start);
            }
            if (parent_transcript.empty())
            {
                std::cerr << "[gene annotation] error: parent '" << attr << "' not found!\n";
            }
            return s;
        }
    }
    return "";
}

string GeneAnnotation::fix_name(string chr_name)
{
    string new_chr_name;
    if (chr_name.compare(0, 3, "chr") == 0)
    {
        return chr_name;
    }
    else if (chr_name.length() > 4) // just fix 1-22, X, Y, MT. ignore contig and ERCC
    {
        return chr_name;
    }
    else
    {
        if (chr_name == "MT")
        {
            new_chr_name = "chrM";
        }
        else
        {
            new_chr_name = "chr" + chr_name;
        }
        return new_chr_name;
    }
}

string GeneAnnotation::get_gene_id(const vector<string> &attributes)
{
    if (anno_source == "gencode")
    {
        return get_gencode_gene_id(attributes);
    }
    else if (anno_source == "refseq")
    {
        return get_refseq_gene_id(attributes);
    }
    else if (anno_source == "plain")
    {
        return get_parent(attributes);
    }
    return "";
}

// string GeneAnnotation::get_plain_gene_id(const vector<string> &attributes)
// {
//     return get_attribute(attributes, "ID");
// }

string GeneAnnotation::get_gencode_gene_id(const vector<string> &attributes)
{
    return get_attribute(attributes, "gene_id");
}

string GeneAnnotation::get_refseq_gene_id(const vector<string> &attributes)
{
    string dbxref = get_attribute(attributes, "Dbxref");

    // GeneID may be missing
    if (dbxref.find("GeneID") == string::npos)
    {
        return "";
    }

    auto start = dbxref.find("GeneID") + 7; // start after "GeneID:"
    auto end = dbxref.find(",", start);
    auto id_length = end - start;

    return dbxref.substr(start, id_length);
}

void GeneAnnotation::parse_anno_entry(
    const bool &fix_chrname,
    const string &line,
    unordered_map<string, unordered_map<string, Gene>> &chr_to_genes_dict,
    unordered_map<string, string> &transcript_to_gene_dict
)
{
    const vector<string> fields = split(line, '\t');
    const vector<string> attributes = split(fields[ATTRIBUTES], ';');

    string chr_name = fields[SEQID];
    const string parent = get_parent(attributes);
    const string type = fields[TYPE];
    const string ID = get_ID(attributes);
    const int strand = get_strand(fields[STRAND][0]);
    const int interval_start = atoi(fields[START].c_str());
    const int interval_end = atoi(fields[END].c_str());

    if (fix_chrname)
    {
        chr_name = fix_name(chr_name);
    }

// DEBUG USE
// cout << "Parsing: " << line << "\n";
// cout << "Type: " << type << " "
//         << "ID: " << ID << " "
//         << "Parent: " << parent << "\n\n";
// DEBUG USE

    string target_gene;
    if (anno_source == "ensembl")
    {
        if (is_exon(fields, attributes))
        {
            if (parent_is_known_transcript(transcript_to_gene_dict, parent))
            {
                target_gene = transcript_to_gene_dict[parent];
            }
            else
            {
                //stringstream err_msg;
                std::cerr << "cannot find grandparent for exon:" << "\n";
                std::cerr << line << "\n";
                throw std::exception();
            }
        }
        else if (is_transcript(fields, attributes))
        {
            if (!ID.empty() && !parent.empty())
            {
                transcript_to_gene_dict[ID] = parent;
            }
            return;
        }
        else if (is_gene(fields, attributes)) {
            recorded_genes.insert(ID);
            return;
        }
    }
    else // anno_source NOT ensembl:
    {
        if (is_transcript(fields, attributes))
        {
            if (!ID.empty() && !parent.empty())
            {
                transcript_to_gene_dict[ID] = parent;
            }
            return;
        }
        else if (is_gene(fields, attributes)) {
            recorded_genes.insert(ID);
            return;
        }

        if (anno_source == "gencode" || anno_source == "refseq")
        {                
            if (type == "exon")
            {
                target_gene = get_gene_id(attributes);
            }
        }
        else // anno_source == plain
        {
            if (type == "CDS")
            {
                target_gene = get_parent_of_CDS(attributes, transcript_to_gene_dict);
                if (recorded_genes.find(target_gene) == recorded_genes.end()) {
                    // it's not a gene, is it a transcript then?
                    if (transcript_to_gene_dict.find(target_gene) != transcript_to_gene_dict.end())
                    {
                        target_gene = transcript_to_gene_dict[target_gene];
                    }
                    else std::cerr << "[gene annotation] warning: target gene for CDS not found\n"
                        << line << std::endl;
                }
            }
        }
    }

    if (!target_gene.empty())
    {
        auto &current_chr = chr_to_genes_dict[chr_name];
        current_chr[target_gene].add_exon(Interval(interval_start, interval_end, strand));
        current_chr[target_gene].set_ID(target_gene);
    }

    return;
}

string GeneAnnotation::guess_anno_source(string gff3_fn)
{
    ifstream infile(gff3_fn);
    string line;

    while (getline(infile, line))
    {
        if (line.find("GENCODE") != string::npos) {
            std::cout << "guessing annotation source: GENCODE" << "\n";
            return "gencode";
        }
        else if (line.find("1\tEnsembl") != string::npos)
        {
            std::cout << "guessing annotation source: ENSEMBL" << "\n";
            return "ensembl";
        }
        else if (line.find("RefSeq\tregion") != string::npos)
        {
            std::cout << "guessing annotation source: RefSeq" << "\n";
            return "refseq";
        }
    }

    std::cout << "Annotation source not recognised, defaulting to 'plain'. Current supported sources: plain, ENSEMBL, GENCODE and RefSeq\n";
    return "plain";
}

const bool GeneAnnotation::parent_is_gene(const string &parent)
{
    return recorded_genes.find(parent) != recorded_genes.end();
}

const bool GeneAnnotation::parent_is_known_transcript(const unordered_map<string, string> &transcript_to_gene_dict, const string &parent)
{
    return transcript_to_gene_dict.find(parent) != transcript_to_gene_dict.end();
}

const bool GeneAnnotation::is_gene(const vector<string> &fields, const vector<string> &attributes)
{
    string type = fields[TYPE];
    if (type.find("gene") != string::npos)
    {
        return true;
    }

    string id = get_attribute(attributes, "ID");
    if (id.find("gene:") != string::npos)
    {
        return true;
    }

    return false;
}

const bool GeneAnnotation::is_exon(const vector<string> &fields, const vector<string> &attributes)
{
    return fields[TYPE] == "exon";
}

const bool GeneAnnotation::is_transcript(const vector<string> &fields, const vector<string> &attributes)
{
    // assume feature is transcript is it has a gene as parent
    return parent_is_gene(get_parent(attributes));
}

void GeneAnnotation::parse_gff3_annotation(string gff3_fn, bool fix_chrname)
{
    ifstream infile(gff3_fn);

    string line;
    unordered_map<string, unordered_map<string, Gene>> chr_to_genes_dict;
    unordered_map<string, string> transcript_to_gene_dict; // store transcript - gene mapping

    // assigned to class member
    anno_source = guess_anno_source(gff3_fn);

    // create transcript-gene mapping
    while (getline(infile, line))
    {
        // skip header lines
        if (line[0] == '#')
        {
            continue;
        }

        parse_anno_entry(fix_chrname, line, chr_to_genes_dict, transcript_to_gene_dict);
    }

    // push genes into annotation class member
    for (auto &chr : chr_to_genes_dict)
    {
        const auto &chr_name = chr.first;
std::cout << "\n\nCHR " << chr_name << std::endl;

        // have one container per path
        vector<Gene> genes_per_chr;

        // merge overlapping exons in each gene
        for (auto &gene : chr.second)
        {
std::cout << "  " << gene.second.gene_id << std::endl;
            gene.second.sort_exon();
            gene.second.flatten_exon();
            genes_per_chr.push_back(gene.second);
            //gene_dict[chr_name].push_back(gene.second);
        }

        //auto &current_genes = gene_dict[chr_name];
        // sort genes based on starting position
        sort(genes_per_chr.begin(), genes_per_chr.end(),
            [] (const Gene &g1, const Gene &g2) { return g1.st < g2.st; }
        );

        // store genes in a path-specific container
        gene_dict[chr_name] = genes_per_chr;

        // create bins of genes
        //bins_dict[chr_name].make_bins(current_genes);

////// output all genes:
std::cout << "  sorted:" << std::endl;
for (Gene &g : gene_dict[chr_name]) {
    std::cout << "  " << g.gene_id << std::endl;
}
//////
    }
}


int GeneAnnotation::ngenes()
{
    int gene_number = 0;
    for (auto iter : gene_dict)
    {
        for (auto sub_iter : iter.second)
        {
            gene_number ++;
        }
    }

    return gene_number;
}


vector<string> GeneAnnotation::get_genelist()
{
    vector<string> gene_list;
    for (auto iter : gene_dict)
    {
        for (auto sub_iter : iter.second)
        {
            gene_list.push_back(sub_iter.gene_id);
        }
    }

    return gene_list;
}


ostream& operator<< (ostream& out, const GeneAnnotation& obj)
{
    out << "annotation statistics:" << "\n";
    for ( const auto& n : obj.gene_dict )
    {
        out << "\t" << "chromosome:[" << n.first << "] number of genes:[" << n.second.size() << "]\n";
    }
    for ( const auto& n : obj.gene_dict )
    {
        out << "first gene in chromosome " << n.first << " :" << "\n";
        out << n.second[0] << "\n";
        //break;
    }
    return out;
}

        }}