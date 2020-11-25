#include "subcommand.hpp"
#include "odgi.hpp"
#include "args.hxx"
#include "algorithms/bin_path_info.hpp"

#include <regex>

namespace odgi {

using namespace odgi::subcommand;

int main_bin(int argc, char** argv) {

    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi bin";
    argv[0] = (char*)prog_name.c_str();
    --argc;

    args::ArgumentParser parser("binning of path information in the graph");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> fa_out_file(parser, "FILE", "store the pangenome sequence in FASTA format in this file", {'f', "fasta"});
    args::ValueFlag<std::string> path_delim(parser, "path-delim", "annotate rows by prefix and suffix of this delimiter", {'D', "path-delim"});
    args::Flag output_json(parser, "write-json", "write JSON format output including additional path positional information, used by Schematize visualization", {'j', "json"});
    args::ValueFlag<std::string> pantograph_file(parser, "FILE", "write JSON format used by PantoGraph in this file; this option activates --no-gap-links", {'p', "pantograph-json"});
    args::Flag aggregate_delim(parser, "aggregate-delim", "aggregate on path prefix delimiter", {'a', "aggregate-delim"});
    args::ValueFlag<uint64_t> num_bins(parser, "N", "number of bins", {'n', "num-bins"});
    args::ValueFlag<uint64_t> bin_width(parser, "bp", "width of each bin in basepairs along the graph vector", {'w', "bin-width"});
    args::Flag write_seqs_not(parser, "write-seqs-not", "don't write out the sequences for each bin", {'s', "no-seqs"});
    args::Flag drop_gap_links(parser, "drop-gap-links", "don't include gap links in the output", {'g', "no-gap-links"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    if (!dg_in_file) {
        std::cerr << "[odgi bin] error: Please specify an input file from where to load the graph via -i=[FILE], --idx=[FILE]." << std::endl;
        return 1;
    }

    graph_t graph;
    assert(argc > 0);
    std::string infile = args::get(dg_in_file);
    if (infile.size()) {
        if (infile == "-") {
            graph.deserialize(std::cin);
        } else {
            ifstream f(infile.c_str());
            graph.deserialize(f);
            f.close();
        }
    }

    std::string delim = args::get(path_delim);
    bool agg_delim = args::get(aggregate_delim);
    auto get_path_prefix = [&](const std::string& path_name) -> std::string {
        if (agg_delim || delim.empty()) {
            return "NA";
        } else {
            return path_name.substr(0, path_name.find(delim));
        }
    };
    auto get_path_suffix = [&](const std::string& path_name) -> std::string {
        if (agg_delim || delim.empty()) {
            return "NA";
        } else {
            return path_name.substr(path_name.find(delim)+1);
        }
    };

    // our aggregation matrix
    std::vector<std::pair<std::string, std::vector<algorithms::path_info_t>>> table;
    if (args::get(num_bins) + args::get(bin_width) == 0) {
        std::cerr << "[odgi bin] error: a bin width or a bin count is required" << std::endl;
        return 1;
    }

    bool no_gap_links = args::get(drop_gap_links);
    std::ofstream pantograph_out;
    std::ofstream xoffset_out;
    if (pantograph_file) {
        std::string pantograph_file_name = args::get(pantograph_file).c_str();
        pantograph_out = std::ofstream(pantograph_file_name);
        std::string xoffset_file_name = pantograph_file_name + ".bin_xoffset.json";
        xoffset_out = std::ofstream(xoffset_file_name);
        no_gap_links = true;
    }

    // JSON VERSIONS
    const uint64_t ODGI_JSON_VERSION = 12; // v12 brings the exact nucleotide positions for each bin for each path referred to as ranges - used by Schematize
    const uint64_t PANTOGRAPH_JSON_VERSION = 1;

    std::function<void(const uint64_t&, const uint64_t&)> write_header_tsv
    = [&] (const uint64_t pangenome_length, const uint64_t bin_width) {
        // no header necessary for tsv so far
    };

    // Schematize json header
    std::function<void(const uint64_t&,
            const uint64_t&)> write_header_json
    = [&] (const uint64_t pangenome_length, const uint64_t bin_width) {
        std::cout << "{\"odgi_version\": " << ODGI_JSON_VERSION << ",";
        std::cout << "\"bin_width\": " << bin_width << ",";
        std::cout << "\"pangenome_length\": " << pangenome_length << "}" << std::endl;
    };

    // PantoGraph json header
    std::function<void(const uint64_t&,
            const uint64_t&)> write_header_pantograph_json
    = [&] (const uint64_t pangenome_length, const uint64_t bin_width) {
        pantograph_out << "{" << std::endl;
        pantograph_out << "\"version\":" << PANTOGRAPH_JSON_VERSION << std::endl;
        pantograph_out << "\"bin_width\":" << bin_width << std::endl;
        pantograph_out << "\"pangenome_length\":" << pangenome_length << "," << std::endl;
        pantograph_out << "\"graph_paths\":[" << std::endl;
    };

    std::function<void(const uint64_t&,
                       const std::string&)> write_seq_json
        = [&](const uint64_t& bin_id, const std::string& seq) {
        if (args::get(write_seqs_not) || fa_out_file) {
            std::cout << "{\"bin_id\":" << bin_id << "}" << std::endl;
        } else {
            std::cout << "{\"bin_id\":" << bin_id << ","
                      << "\"sequence\":\"" << seq << "\"}" << std::endl;
        }
    };

    std::function<void(const std::string&)> write_fasta
            = [&](const std::string& nuc_seq) {
                if (fa_out_file) {
                    std::ofstream out(args::get(fa_out_file));
                    std::string fa_out_name = args::get(fa_out_file).c_str();
                    std::regex regex("/");
                    std::vector<std::string> splitted(
                            std::sregex_token_iterator(fa_out_name.begin(), fa_out_name.end(), regex, -1),
                            std::sregex_token_iterator()
                            );
                    fa_out_name = splitted[splitted.size() - 1];
                    // Write header
                    out << ">" << fa_out_name << std::endl;
                    // Write the actual sequences, 80 nucleotides per line
                    for (unsigned i = 0; i < nuc_seq.length(); i += 80) {
                        std:: string sub_nuc_seq = nuc_seq.substr(i, 80);
                        out << sub_nuc_seq << std::endl;
                    }
                }
            };

    // for Schematize:
    std::function<void(const vector<std::pair<uint64_t , uint64_t >>&)> write_ranges_json
        = [&](const vector<std::pair<uint64_t , uint64_t >>& ranges) {
        std::cout << "[";
        for (int i = 0; i < ranges.size(); i++) {
            std::pair<uint64_t, uint64_t > range = ranges[i];
            if (i == 0) {
                std::cout << "[" << range.first << "," << range.second << "]";
            } else {
                std::cout << "," << "[" << range.first << "," << range.second << "]";
            }
        }
        std::cout << "]";
    };

    // for PantoGraph:
    std::function<void(const vector<std::pair<uint64_t , uint64_t >>&)> write_ranges_pantograph_json
        = [&](const vector<std::pair<uint64_t , uint64_t >>& ranges) {
        pantograph_out << "\"ranges\":[";
        for (int i = 0; i < ranges.size(); i++) {
            std::pair<uint64_t, uint64_t > range = ranges[i];
            if (i == 0) {
                pantograph_out << "{\"start\":" << range.first << ",\"end\":" << range.second << "}";
            } else {
                pantograph_out << "," << "{\"start\": " << range.first << ",\"end\":" << range.second << "}";
            }
        }
        pantograph_out << "]";
    };

    // for Schematize:
    std::function<void(const std::string&,
                       const std::vector<std::pair<uint64_t, uint64_t>>&,
                       const std::map<uint64_t, algorithms::path_info_t>&)> write_json
        = [&](const std::string& path_name,
              const std::vector<std::pair<uint64_t, uint64_t>>& links,
              const std::map<uint64_t, algorithms::path_info_t>& bins) {
        std::string name_prefix = get_path_prefix(path_name);
        std::string name_suffix = get_path_suffix(path_name);
        std::cout << "{\"path_name\":\"" << path_name << "\",";
        if (!delim.empty()) {
            std::cout << "\"path_name_prefix\":\"" << name_prefix << "\","
                      << "\"path_name_suffix\":\"" << name_suffix << "\",";
        }
        std::cout << "\"bins\":[";
        auto entry_it = bins.begin();
        for (uint64_t i = 0; i < bins.size(); ++i) {
            auto& bin_id = entry_it->first;
            auto& info = entry_it->second;
            std::cout << "[" << bin_id << ","
                      << info.mean_cov << ","
                      << info.mean_inv << ","
                      << info.mean_pos << ",";
            write_ranges_json(info.ranges);
			std::cout  << "]";
            if (i+1 != bins.size()) {
                std::cout << ",";
            }
            ++entry_it;
        }
        std::cout << "],";
        std::cout << "\"links\":[";
        for (uint64_t i = 0; i < links.size(); ++i) {
            auto& link = links[i];
            std::cout << "[" << link.first << "," << link.second << "]";
            if (i+1 < links.size()) std::cout << ",";
        }
        std::cout << "]}" << std::endl;
    };

    // for PantoGraph:
    std::function<void(const std::string&,
                       const std::vector<std::pair<uint64_t, uint64_t>>&,
                       const std::map<uint64_t, algorithms::path_info_t>&)> write_pantograph_json
        = [&](const std::string& path_name,
              const std::vector<std::pair<uint64_t, uint64_t>>& links,
              const std::map<uint64_t, algorithms::path_info_t>& bins) {
        std::string name_prefix = get_path_prefix(path_name);
        std::string name_suffix = get_path_suffix(path_name);
        pantograph_out << "{\"path_name\":\"" << path_name << "\",";
        if (!delim.empty()) {
            pantograph_out << "\"path_name_prefix\":\"" << name_prefix << "\","
                      << "\"path_name_suffix\":\"" << name_suffix << "\",";
        }
        pantograph_out << std::endl << "\"bins\":[" << std::endl;
        auto entry_it = bins.begin();
        for (uint64_t i = 0; i < bins.size(); ++i) {
            auto& bin_id = entry_it->first;
            auto& info = entry_it->second;
            pantograph_out << "{\"bin_id\":" << bin_id << ","
                      << "\"cov\":" << info.mean_cov << ","
                      << "\"inv\":" << info.mean_inv << ","
                      << "\"pos\":" << info.mean_pos << ",";
            write_ranges_pantograph_json(info.ranges);
			pantograph_out << "}";
            if (i+1 != bins.size()) {
                pantograph_out << "," << std::endl;
            }
            ++entry_it;
        }
        pantograph_out << std::endl << "]," << std::endl;
        pantograph_out << "\"links\":[" << std::endl;
        for (uint64_t i = 0; i < links.size(); ++i) {
            auto& link = links[i];
            pantograph_out << "{" << "\"from\":" << link.first << "," 
                      << "\"to\":" << link.second << "}";
            if (i+1 < links.size()) pantograph_out << ",";
        }
        pantograph_out << std::endl << "]" << std::endl;
        pantograph_out << "}" << std::endl;
    };

    std::function<void(const uint64_t&,
                       const std::string&)> write_seq_noop
        = [&](const uint64_t& bin_id, const std::string& seq) {
    };

    std::function<void(const std::string&,
                       const std::vector<std::pair<uint64_t, uint64_t>>&,
                       const std::map<uint64_t, algorithms::path_info_t>&)> write_tsv
        = [&](const std::string& path_name,
              const std::vector<std::pair<uint64_t, uint64_t>>& links,
              const std::map<uint64_t, algorithms::path_info_t>& bins) {
        std::string name_prefix = get_path_prefix(path_name);
        std::string name_suffix = get_path_suffix(path_name);
        for (auto& entry : bins) {
            auto& bin_id = entry.first;
            auto& info = entry.second;
            if (info.mean_cov) {
                std::cout << path_name << "\t"
                          << name_prefix << "\t"
                          << name_suffix << "\t"
                          << bin_id << "\t"
                          << info.mean_cov << "\t"
                          << info.mean_inv << "\t"
                          << info.mean_pos << "\t"
						  << info.ranges[0].first << "\t";
                if (info.ranges[info.ranges.size() - 1].second == 0) {
                    std::cout << info.ranges[info.ranges.size() - 1].first << std::endl;
                } else {
                    std::cout << info.ranges[info.ranges.size() - 1].second << std::endl;
                }
            }
        }
    };

    // for Schematize:
    std::function<void(const uint64_t&,
                       const uint64_t&,
                       const uint64_t&)> write_xoffset_noop
        = [&](const uint64_t& bin_id, const uint64_t& offset, const uint64_t& max_bin) {
    };

    // for Pantograph:
    std::function<void(const uint64_t&,
                       const uint64_t&,
                       const uint64_t&)> write_xoffset
        = [&](const uint64_t& bin_id, const uint64_t& offset, const uint64_t& max_bin) {
        if (bin_id == 0) {
            xoffset_out << "{" << std::endl << "[";
        } else {
            xoffset_out << ",";
        }
        xoffset_out << offset;
        if (bin_id == max_bin) {
            xoffset_out << "]" << std::endl << "}" << std::endl;
        }
    };

    if (args::get(output_json)) {
        algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                  write_header_json, write_json, write_seq_json, write_fasta, write_xoffset_noop,
                                  args::get(num_bins), args::get(bin_width), no_gap_links);
    }
    if (pantograph_file) {
        algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                  write_header_pantograph_json, write_pantograph_json, write_seq_noop, write_fasta, write_xoffset,
                                  args::get(num_bins), args::get(bin_width), no_gap_links);
        pantograph_out << "]}" << std::endl;
    }
    if (!args::get(output_json) && !pantograph_file) {
        std::cout << "path.name" << "\t"
                  << "path.prefix" << "\t"
                  << "path.suffix" << "\t"
                  << "bin" << "\t"
                  << "mean.cov" << "\t"
                  << "mean.inv" << "\t"
                  << "mean.pos" << "\t"
                  << "first.nucl" << "\t"
                  << "last.nucl" << std::endl;
        algorithms::bin_path_info(graph, (args::get(aggregate_delim) ? args::get(path_delim) : ""),
                                  write_header_tsv,write_tsv, write_seq_noop, write_fasta, write_xoffset_noop,
                                  args::get(num_bins), args::get(bin_width), no_gap_links);
    }
    return 0;
}

static Subcommand odgi_bin("bin", "bin path information across the graph",
                              PIPELINE, 3, main_bin);


}
