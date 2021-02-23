#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/bin_path_info.hpp"
#include "gene_anno/GeneAnnotation.hpp"

#include <regex>
#include <filesystem>

namespace odgi
{

    using namespace odgi::subcommand;

    std::string currentDateTime()
    {
        time_t now = time(0);
        struct tm tstruct;
        char buf[80];
        tstruct = *localtime(&now);
        // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
        // for more information about date/time format
        strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

        return buf;
    }

    int main_bin(int argc, char **argv)
    {

        for (uint64_t i = 1; i < argc - 1; ++i)
        {
            argv[i] = argv[i + 1];
        }
        std::string prog_name = "odgi bin";
        argv[0] = (char *)prog_name.c_str();
        --argc;

        args::ArgumentParser parser("binning of path information in the graph");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
        args::ValueFlag<std::string> fa_out_file(parser, "FILE", "store the pangenome sequence in FASTA format in this file", {'f', "fasta"});
        args::Flag output_json(parser, "write-json", "write JSON format output including additional path positional information, used by Schematize visualization", {'j', "json"});
        args::ValueFlag<std::string> pantograph_file(parser, "FILE", "write JSON format used by PantoGraph in this file; this option activates --no-gap-links", {'p', "pantograph-json"}); // TODO: change into folder
        args::ValueFlagList<std::string> annotation_files(parser, "FILEs", "gene annotation file(s) in gff3 format of the reference genome specified with --ref-genome. Can be repeated", {'A', "anno-file"});
        args::ValueFlagList<std::string> ref_genome_IDs(parser, "IDs", "reference genome ID(s) as specified in the odgi graph, on which the gene annotation is based on. Can be repeated", {'r', "ref-genome"});
        args::ValueFlag<std::string> path_delim(parser, "path-delim", "annotate rows by prefix and suffix of this delimiter", {'D', "path-delim"});
        args::Flag aggregate_delim(parser, "aggregate-delim", "aggregate on path prefix delimiter", {'a', "aggregate-delim"});
        args::ValueFlag<uint64_t> num_bins(parser, "N", "number of bins", {'n', "num-bins"});
        args::ValueFlag<uint64_t> bin_width(parser, "bp", "width of each bin in basepairs along the graph vector", {'w', "bin-width"});
        args::ValueFlag<uint64_t> bins_per_file(parser, "N", "number of bins per path and chunk file, default: 200", {'c', "bins-per-file"});
        args::Flag write_seqs_not(parser, "write-seqs-not", "don't write out the sequences for each bin", {'s', "no-seqs"});
        args::Flag drop_gap_links(parser, "drop-gap-links", "don't include gap links in the output", {'g', "no-gap-links"});
        try
        {
            parser.ParseCLI(argc, argv);
        }
        catch (args::Help)
        {
            std::cout << parser;
            return 0;
        }
        catch (args::ParseError e)
        {
            std::cerr << e.what() << std::endl;
            std::cerr << parser;
            return 1;
        }
        if (argc == 1)
        {
            std::cout << parser;
            return 1;
        }

        if (!dg_in_file)
        {
            std::cerr << "[odgi bin] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
            return 1;
        }

        //////////////////////////////
        std::map<std::string, std::string> ref_paths;
        if (annotation_files)
        {
            if (!ref_genome_IDs)
            {
                std::cerr << "[odgi bin] error: Reference genome IDs needed for annotation files (-r/--ref-genome option missing)" << std::endl;
                return 1;
            }
            std::vector<std::string> anno_files = args::get(annotation_files);
            std::vector<std::string> refIDs = args::get(ref_genome_IDs);
            if (anno_files.size() != refIDs.size())
            {
                std::cerr << "[odgi bin] error: Unequal number of annotation files and reference genome IDs!" << std::endl;
                return 1;
            }
            size_t i = 0;
            for (auto a : anno_files)
            {
                //std::string annofile = (args::get(annotation_file))[0];
                if (std::filesystem::exists(a))
                {
                    ref_paths[refIDs[i]] = a;
                    i++;
                    // odgi::gene_anno::GeneAnnotation anno;
                    // anno.parse_gff3_annotation(a, false);
                    // std::cout << anno << std::endl;
                    // return 0;
                }
                else
                {
                    std::cerr << "[odgi bin] error: The specified gene annotation file (-A, --anno-file) does not exist." << std::endl;
                    return 1;
                }
            }
            // for (auto p: ref_paths) {
            //     std::cout << p.first << ": " << p.second << std::endl;
            // }
            // return 0;
        }
        //////////////////////////////

        graph_t graph;
        assert(argc > 0);
        std::string infile = args::get(dg_in_file);
        if (infile.size())
        {
            if (infile == "-")
            {
                graph.deserialize(std::cin);
            }
            else
            {
                ifstream f(infile.c_str());
                graph.deserialize(f);
                f.close();
            }
        }

        std::string delim = args::get(path_delim);
        bool agg_delim = args::get(aggregate_delim);
        auto get_path_prefix = [&](const std::string &path_name) -> std::string {
            if (agg_delim || delim.empty())
            {
                return "NA";
            }
            else
            {
                return path_name.substr(0, path_name.find(delim));
            }
        };
        auto get_path_suffix = [&](const std::string &path_name) -> std::string {
            if (agg_delim || delim.empty())
            {
                return "NA";
            }
            else
            {
                return path_name.substr(path_name.find(delim) + 1);
            }
        };

        // our aggregation matrix
        std::vector<std::pair<std::string, std::vector<algorithms::path_info_t>>> table;
        if (args::get(num_bins) + args::get(bin_width) == 0)
        {
            std::cerr << "[odgi bin] error: a bin width or a bin count is required" << std::endl;
            return 1;
        }

        // default value for bins_per_file
        uint64_t bins_per_chunk = args::get(bins_per_file);
        if (!args::get(bins_per_file))
        {
            bins_per_chunk = 200;
        }

        bool no_gap_links = args::get(drop_gap_links);
        std::ofstream pantograph_out; // deprecated, only used for pantograph1
        std::ofstream bin2file_out;
        std::vector<std::ofstream> file_handles;
        std::string pantograph_file_name;
        // std::ofstream xoffset_out;
        if (pantograph_file)
        {
            pantograph_file_name = args::get(pantograph_file).c_str();
            pantograph_out = std::ofstream(pantograph_file_name);
            std::string bin2file_file_name = pantograph_file_name + ".bin" + to_string(bin_width) + ".bin2file.json";
            bin2file_out = std::ofstream(bin2file_file_name);
            // std::string xoffset_file_name = pantograph_file_name + ".bin_xoffset.json";
            // xoffset_out = std::ofstream(xoffset_file_name);
            no_gap_links = true;
        }

        // JSON VERSIONS
        const uint64_t ODGI_JSON_VERSION = 12; // v12 brings the exact nucleotide positions for each bin for each path referred to as ranges - used by Schematize
        const uint64_t PANTOGRAPH_JSON_VERSION = 1;

        //////////// FUNCTIONS

        std::function<void(const uint64_t &, const uint64_t &)> write_header_tsv = [&](const uint64_t pangenome_length, const uint64_t bin_width) {
            // no header necessary for tsv so far
        };

        // Schematize json header
        std::function<void(const uint64_t &,
                           const uint64_t &)>
            write_header_json = [&](const uint64_t pangenome_length, const uint64_t bin_width) {
                std::cout << "{\"odgi_version\": " << ODGI_JSON_VERSION << ",";
                std::cout << "\"bin_width\": " << bin_width << ",";
                std::cout << "\"pangenome_length\": " << pangenome_length << "}" << std::endl;
            };

        // PantoGraph json header
        std::function<void(const uint64_t &, const uint64_t &)> write_header_pantograph_json = [&](const uint64_t bin_width, const uint64_t num_bins) {
            bin2file_out << "{" << std::endl
                         << "\t\"version\":" << PANTOGRAPH_JSON_VERSION << "," << std::endl;
            bin2file_out << "\t\"pangenome_len\":" << num_bins << "," << std::endl;
            bin2file_out << "\t\"zoom_level\":" << bin_width << "," << std::endl;
            bin2file_out << "\t\"files\": [" << std::endl;

            for (uint64_t bin = 0; bin < num_bins; bin += bins_per_chunk)
            {
                uint64_t chunk = bin / bins_per_chunk;
                //std::ofstream file = std::ofstream(pantograph_file_name + "." + to_string(chunk));
                std::string filename = pantograph_file_name + "." + to_string(chunk);
                file_handles.emplace_back(std::ofstream{filename}); // https://stackoverflow.com/questions/29004665/create-a-vector-of-ofstreams
                file_handles.back() << "{" << std::endl;
                file_handles.back() << "\"version\":" << PANTOGRAPH_JSON_VERSION << "," << std::endl;
                file_handles.back() << "\"bin_width\":" << bin_width << "," << std::endl;
                file_handles.back() << "\"first_bin\":" << bin + 1 << "," << std::endl;
                file_handles.back() << "\"last_bin\":" << std::min(bin + bins_per_chunk, num_bins) << "," << std::endl;
                file_handles.back() << "\"pangenome_len\":" << num_bins << "," << std::endl;
                file_handles.back() << "\"graph_paths\":[" << std::endl;

                if (bin > 0)
                    bin2file_out << "," << std::endl;
                bin2file_out << "\t\t{" << std::endl;
                bin2file_out << "\t\t\t\"file\": \"" << filename << "\"," << std::endl;
                bin2file_out << "\t\t\t\"first_bin\": \"" << bin + 1 << "\"," << std::endl;
                bin2file_out << "\t\t\t\"last_bin\": \"" << std::min(bin + bins_per_chunk, num_bins) << "\"" << std::endl
                             << "\t\t}";
            }

            bin2file_out << std::endl
                         << "\t]" << std::endl
                         << "}" << std::endl;
        };

        std::function<void(const uint64_t &,
                           const std::string &)>
            write_seq_json = [&](const uint64_t &bin_id, const std::string &seq) {
                if (args::get(write_seqs_not) || fa_out_file)
                {
                    std::cout << "{\"bin_id\":" << bin_id << "}" << std::endl;
                }
                else
                {
                    std::cout << "{\"bin_id\":" << bin_id << ","
                              << "\"sequence\":\"" << seq << "\"}" << std::endl;
                }
            };

        // for Schematize:
        std::function<void(const std::string &)> write_fasta = [&](const std::string &nuc_seq) {
            if (fa_out_file)
            {
                std::ofstream out(args::get(fa_out_file));
                std::string fa_out_name = args::get(fa_out_file).c_str();
                std::regex regex("/");
                std::vector<std::string> splitted(
                    std::sregex_token_iterator(fa_out_name.begin(), fa_out_name.end(), regex, -1),
                    std::sregex_token_iterator());
                fa_out_name = splitted[splitted.size() - 1];
                // Write header
                out << ">" << fa_out_name << std::endl;
                // Write the actual sequences, 80 nucleotides per line
                for (unsigned i = 0; i < nuc_seq.length(); i += 80)
                {
                    std::string sub_nuc_seq = nuc_seq.substr(i, 80);
                    out << sub_nuc_seq << std::endl;
                }
            }
        };

        // for Pantograph:
        std::function<void(const std::string &, const uint64_t &, const uint64_t &)> write_seq_pantograph = [&](const std::string &nuc_seq, const uint64_t bin_width, const uint64_t &pangenome_len) {
            for (uint64_t nuc = 0; nuc < pangenome_len; nuc += bins_per_chunk * bin_width)
            {
                uint64_t chunk = nuc / bin_width / bins_per_chunk;
                //std::cout << " chunk " << chunk << std::endl;
                uint64_t nr_nucs = std::min(bin_width * bins_per_chunk, pangenome_len - nuc);
                file_handles[chunk] << "," << std::endl
                                    << "\"sequence\": \""
                                    << nuc_seq.substr(nuc, nr_nucs) << "\""
                                    << std::endl
                                    << "}" << std::endl;
            }
        };

        // for Schematize:
        std::function<void(const vector<std::pair<uint64_t, uint64_t>> &)> write_ranges_json = [&](const vector<std::pair<uint64_t, uint64_t>> &ranges) {
            std::cout << "[";
            for (int i = 0; i < ranges.size(); i++)
            {
                std::pair<uint64_t, uint64_t> range = ranges[i];
                if (i == 0)
                {
                    std::cout << "[" << range.first << "," << range.second << "]";
                }
                else
                {
                    std::cout << ","
                              << "[" << range.first << "," << range.second << "]";
                }
            }
            std::cout << "]";
        };

        // for PantoGraph:
        std::function<void(const vector<std::pair<uint64_t, uint64_t>> &)> write_ranges_pantograph_json = [&](const vector<std::pair<uint64_t, uint64_t>> &ranges) {
            pantograph_out << "\"ranges\":[";
            for (int i = 0; i < ranges.size(); i++)
            {
                std::pair<uint64_t, uint64_t> range = ranges[i];
                if (i == 0)
                {
                    pantograph_out << "{\"start\":" << range.first << ",\"end\":" << range.second << "}";
                }
                else
                {
                    pantograph_out << ","
                                   << "{\"start\":" << range.first << ",\"end\":" << range.second << "}";
                }
            }
            pantograph_out << "]";
        };

        // for PantoGraph2:
        std::function<std::string(const vector<std::pair<uint64_t, uint64_t>> &)> get_ranges_pantograph_json = [&](const vector<std::pair<uint64_t, uint64_t>> &ranges) {
            std::string out = "";
            out += "\"ranges\":[";
            for (int i = 0; i < ranges.size(); i++)
            {
                std::pair<uint64_t, uint64_t> range = ranges[i];
                if (i == 0)
                {
                    //out += "{\"start\":" + std::to_string(range.first) + ",\"end\":" + std::to_string(range.second) + "}";
                    out += "{" + std::to_string(range.first) + "," + std::to_string(range.second) + "}";
                }
                else
                {
                    //out += ",{\"start\":" + std::to_string(range.first) + ",\"end\":" + std::to_string(range.second) + "}";
                    out += ",{" + std::to_string(range.first) + "," + std::to_string(range.second) + "}";
                }
            }
            out += "]";
            return (out);
        };

        // for Schematize:
        std::function<void(const std::string &,
                           const std::vector<std::pair<uint64_t, uint64_t>> &,
                           const std::map<uint64_t, algorithms::path_info_t> &)>
            write_json = [&](const std::string &path_name,
                             const std::vector<std::pair<uint64_t, uint64_t>> &links,
                             const std::map<uint64_t, algorithms::path_info_t> &bins) {
                std::string name_prefix = get_path_prefix(path_name);
                std::string name_suffix = get_path_suffix(path_name);
                std::cout << "{\"path_name\":\"" << path_name << "\",";
                if (!delim.empty())
                {
                    std::cout << "\"path_name_prefix\":\"" << name_prefix << "\","
                              << "\"path_name_suffix\":\"" << name_suffix << "\",";
                }
                std::cout << "\"bins\":[";
                auto entry_it = bins.begin();
                for (uint64_t i = 0; i < bins.size(); ++i)
                {
                    auto &bin_id = entry_it->first;
                    auto &info = entry_it->second;
                    std::cout << "[" << bin_id << ","
                              << info.mean_cov << ","
                              << info.mean_inv << ","
                              << info.mean_pos << ",";
                    write_ranges_json(info.ranges);
                    std::cout << "]";
                    if (i + 1 != bins.size())
                    {
                        std::cout << ",";
                    }
                    ++entry_it;
                }
                std::cout << "],";
                std::cout << "\"links\":[";
                for (uint64_t i = 0; i < links.size(); ++i)
                {
                    auto &link = links[i];
                    std::cout << "[" << link.first << "," << link.second << "]";
                    if (i + 1 < links.size())
                        std::cout << ",";
                }
                std::cout << "]}" << std::endl;
            };

        // for PantoGraph:
        std::function<void(const std::string &,
                           const std::vector<std::pair<uint64_t, uint64_t>> &,
                           const std::map<uint64_t, algorithms::bin_info_t> &)>
            write_pantograph_json = [&](const std::string &path_name,
                                        const std::vector<std::pair<uint64_t, uint64_t>> &links,
                                        const std::map<uint64_t, algorithms::bin_info_t> &bins) {
                std::string name_prefix = get_path_prefix(path_name);
                std::string name_suffix = get_path_suffix(path_name);
                pantograph_out << "{\"path_name\":\"" << path_name << "\",";
                if (!delim.empty())
                {
                    pantograph_out << "\"path_name_prefix\":\"" << name_prefix << "\","
                                   << "\"path_name_suffix\":\"" << name_suffix << "\",";
                }
                pantograph_out << std::endl
                               << "\"bins\":[" << std::endl;
                auto entry_it = bins.begin();
                for (uint64_t i = 0; i < bins.size(); ++i)
                {
                    auto &bin_id = entry_it->first;
                    auto &info = entry_it->second;
                    pantograph_out << "{\"bin_id\":" << bin_id << ","
                                   << "\"cov\":" << info.mean_cov << ","
                                   << "\"inv\":" << info.mean_inv << ","
                                   << "\"pos\":" << info.mean_pos << ",";
                    write_ranges_pantograph_json(info.ranges);
                    pantograph_out << "}";
                    if (i + 1 != bins.size())
                    {
                        pantograph_out << "," << std::endl;
                    }
                    ++entry_it;
                }
                pantograph_out << std::endl
                               << "]," << std::endl;
                pantograph_out << "\"links\":[" << std::endl;
                for (uint64_t i = 0; i < links.size(); ++i)
                {
                    auto &link = links[i];
                    pantograph_out << "{"
                                   << "\"from\":" << link.first << ","
                                   << "\"to\":" << link.second << "}";
                    if (i + 1 < links.size())
                        pantograph_out << ",";
                }
                pantograph_out << std::endl
                               << "]" << std::endl;
                pantograph_out << "}" << std::endl;
            };

        // for PantoGraph2:
        std::function<void(const std::string &,
                           const std::map<uint64_t, algorithms::bin_info_t> &,
                           const bool &, const uint64_t &)>
            write_pantograph_json2 = [&](const std::string &path_name,
                                         const std::map<uint64_t, algorithms::bin_info_t> &bins,
                                         const bool &first_path, const uint64_t &total_bins) {
                std::cout << "write path " << path_name << " @ " << currentDateTime() << std::endl;
                std::string name_prefix = get_path_prefix(path_name);
                std::string name_suffix = get_path_suffix(path_name);

                for (uint64_t bin = 0; bin < total_bins; bin += bins_per_chunk)
                {
                    uint64_t chunk = bin / bins_per_chunk;
                    //std::cout << "  chunk " << chunk << " @bin " << bin << std::endl;
                    if (!first_path)
                        file_handles[chunk] << "," << std::endl;
                    file_handles[chunk] << "{\"path_name\":\"" << path_name << "\",";
                    if (!delim.empty())
                    {
                        file_handles[chunk] << "\"path_name_prefix\":\"" << name_prefix << "\","
                                            << "\"path_name_suffix\":\"" << name_suffix << "\",";
                    }
                    //std::cout << "  path " << path_name << std::endl;
                    file_handles[chunk] << std::endl
                                        << "\"bins\":[" << std::endl;
                    bool first_bin_of_path = true;
                    for (uint64_t i = bin + 1; i <= std::min(bin + bins_per_chunk, total_bins); ++i)
                    {
                        //uint64_t bin_id = bins[bin].first;
                        if (bins.find(i) != bins.end())
                        {
                            //std::cout << "    bin " << i << std::endl;
                            if (!first_bin_of_path)
                            {
                                file_handles[chunk] << "," << std::endl;
                                first_bin_of_path = false;
                            }
                            algorithms::bin_info_t info = bins.at(i);
                            file_handles[chunk] << "{\"bin\":" << i << ","
                                                << "\"cov\":" << info.mean_cov << ","
                                                << "\"inv\":" << info.mean_inv << ","
                                                << "\"pos\":" << info.mean_pos << ",";
                            file_handles[chunk] << get_ranges_pantograph_json(info.ranges);
                            file_handles[chunk] << ",\"links_to\":[";
                            //std::cout << "    links_to " << info.links_to.size() << std::endl;
                            for (uint64_t li = 0; li < info.links_to.size(); ++li)
                            {
                                auto &link = info.links_to[li];
                                file_handles[chunk] << "{" << link.from << "," << link.to
                                                    << "," << link.times << "}";
                                if (li + 1 < info.links_to.size())
                                    file_handles[chunk] << ",";
                            }
                            //std::cout << "    links_from " << info.links_from.size() << std::endl;
                            file_handles[chunk] << "],\"links_from\":[";
                            for (uint64_t li = 0; li < info.links_from.size(); ++li)
                            {
                                auto &link = info.links_from[li];
                                file_handles[chunk] << "{" << link.from << "," << link.to
                                                    << "," << link.times << "}";
                                if (li + 1 < info.links_from.size())
                                    file_handles[chunk] << ",";
                            }
                            //std::cout << "    links ok " << std::endl;
                            file_handles[chunk] << "]}";
                            // if (i != num_bins && i % bins_per_chunk != 0) {
                            //     file_handles[chunk] << "," << std::endl;
                            // }
                            if (info.genes.size() > 0)
                            {
                                file_handles[chunk] << ",\"genes\":[";
                                for (auto &g : info.genes)
                                {
                                    // g.first = gene ID
                                    // g.second = gene_info_t
                                    if (g.first != info.genes.begin()->first)
                                        file_handles[chunk] << ",";
                                    file_handles[chunk] << "{\"id\":" << g.first << ",\"str\":" << g.second.strand << ",\"ex\":[";
                                    for (auto &ex : g.second.exons)
                                    {
                                        if (ex != *g.second.exons.begin())
                                            file_handles[chunk] << ",";
                                        if (ex > 0)
                                            file_handles[chunk] << ex;
                                    }
                                    file_handles[chunk] << "]}";
                                }
                                file_handles[chunk] << "]";
                            }
                        }
                    }
                    if (!first_bin_of_path)
                        file_handles[chunk] << std::endl;
                    file_handles[chunk] << "]}"; // bins + path
                }
                std::cout << "  finished path " << path_name << " @ " << currentDateTime() << std::endl;
            };

        std::function<void(const uint64_t &,
                           const std::string &)>
            write_seq_noop = [&](const uint64_t &bin_id, const std::string &seq) {
            };

        std::function<void(const std::string &,
                           const std::vector<std::pair<uint64_t, uint64_t>> &,
                           const std::map<uint64_t, algorithms::path_info_t> &)>
            write_tsv = [&](const std::string &path_name,
                            const std::vector<std::pair<uint64_t, uint64_t>> &links,
                            const std::map<uint64_t, algorithms::path_info_t> &bins) {
                std::string name_prefix = get_path_prefix(path_name);
                std::string name_suffix = get_path_suffix(path_name);
                for (auto &entry : bins)
                {
                    auto &bin_id = entry.first;
                    auto &info = entry.second;
                    if (info.mean_cov)
                    {
                        std::cout << path_name << "\t"
                                  << name_prefix << "\t"
                                  << name_suffix << "\t"
                                  << bin_id << "\t"
                                  << info.mean_cov << "\t"
                                  << info.mean_inv << "\t"
                                  << info.mean_pos << "\t"
                                  << info.ranges[0].first << "\t";
                        if (info.ranges[info.ranges.size() - 1].second == 0)
                        {
                            std::cout << info.ranges[info.ranges.size() - 1].first << std::endl;
                        }
                        else
                        {
                            std::cout << info.ranges[info.ranges.size() - 1].second << std::endl;
                        }
                    }
                }
            };

        // for Pantograph:
        std::function<void(const uint64_t &,
                           const uint64_t &,
                           const uint64_t &)>
            write_xoffset = [&](const uint64_t &bin_id, const uint64_t &offset, const uint64_t &max_bin) {
                if (bin_id == 1)
                {
                    pantograph_out << "]," << std::endl; // end of graph_paths
                    pantograph_out << "\"xoffsets\":[";
                }
                else
                {
                    pantograph_out << ",";
                }
                pantograph_out << offset;
                if (bin_id == max_bin)
                {
                    pantograph_out << "]";
                }
            };

        // for Pantograph2:
        std::function<void(const std::map<uint64_t, uint64_t> &,
                           const uint64_t &)>
            write_xoffset2 = [&](const std::map<uint64_t, uint64_t> &offsets, const uint64_t &num_bins) {
                //std::cout << "write xoffsets\n";
                for (uint64_t bin = 0; bin < num_bins; bin++)
                {
                    uint64_t chunk = bin / bins_per_chunk;
                    //std::cout << "chunk " << chunk << " @bin " << bin << std::endl;
                    if (bin % bins_per_chunk == 0)
                    {
                        file_handles[chunk] << "]," << std::endl; // end of graph_paths
                        file_handles[chunk] << "\"xoffsets\":[";
                    }
                    else
                    {
                        file_handles[chunk] << ",";
                    }
                    file_handles[chunk] << offsets.at(bin + 1);
                    //std::cout << "  offset " << offsets.at(bin+1) << std::endl;
                    if ((bin + 1) % bins_per_chunk == 0 || (bin + 1) == num_bins)
                    {
                        file_handles[chunk] << "]";
                    }
                }
            };

        if (args::get(output_json))
        {
            algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                      write_header_json, write_json, write_seq_json, write_fasta,
                                      args::get(num_bins), args::get(bin_width), no_gap_links);
            std::cout << "schematize format completed" << std::endl;
        }
        if (pantograph_file)
        {
            algorithms::bin_path_info_for_pantograph(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                                     write_header_pantograph_json, write_pantograph_json2, write_seq_noop, write_seq_pantograph, write_xoffset2,
                                                     ref_paths, args::get(num_bins), args::get(bin_width));
        }
        if (!args::get(output_json) && !pantograph_file)
        {
            std::cout << "path.name"
                      << "\t"
                      << "path.prefix"
                      << "\t"
                      << "path.suffix"
                      << "\t"
                      << "bin"
                      << "\t"
                      << "mean.cov"
                      << "\t"
                      << "mean.inv"
                      << "\t"
                      << "mean.pos"
                      << "\t"
                      << "first.nucl"
                      << "\t"
                      << "last.nucl" << std::endl;
            algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                      write_header_tsv, write_tsv, write_seq_noop, write_fasta,
                                      args::get(num_bins), args::get(bin_width), no_gap_links);
            std::cout << "pantograph format completed" << std::endl;
        }
        return 0;
    }

    static Subcommand odgi_bin("bin",
                               "bin path information across the graph",
                               PIPELINE, 3, main_bin);

}
