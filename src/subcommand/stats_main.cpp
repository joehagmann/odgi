#include "subcommand.hpp"
#include "odgi.hpp"
#include "IntervalTree.h"
//#include "gfakluge.hpp"
#include "args.hxx"
#include "split.hpp"
//#include "io_helper.hpp"
#include <omp.h>
#include "algorithms/layout.hpp"

//#define debug_odgi_stats

namespace odgi {

using namespace odgi::subcommand;

int main_stats(int argc, char** argv) {

    // trick argumentparser to do the right thing with the subcommand
    for (uint64_t i = 1; i < argc-1; ++i) {
        argv[i] = argv[i+1];
    }
    std::string prog_name = "odgi stats";
    argv[0] = (char*)prog_name.c_str();
    --argc;
    
    args::ArgumentParser parser("metrics describing variation graphs and their path relationships");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> dg_in_file(parser, "FILE", "load the variation graph from this file", {'i', "idx"});
    args::ValueFlag<std::string> layout_in_file(parser, "FILE", "read the layout coordinates from this file", {'c', "coords-in"});
    args::Flag summarize(parser, "summarize", "summarize the graph properties and dimensions", {'S', "summarize"});
    args::Flag base_content(parser, "base-content", "describe the base content of the graph", {'b', "base-content"});
    args::Flag path_coverage(parser, "coverage", "provide a histogram of path coverage over bases in the graph", {'C', "coverage"});
    args::Flag path_setcov(parser, "setcov", "provide a histogram of coverage over unique sets of paths", {'V', "set-coverage"});
    args::Flag path_setcov_count(parser, "setcountcov", "provide a histogram of coverage over counts of unique paths", {'Q', "set-count-coverage"});
    args::Flag path_multicov(parser, "multicov", "provide a histogram of coverage over unique multisets of paths", {'M', "multi-coverage"});
    args::Flag path_multicov_count(parser, "multicountcov", "provide a histogram of coverage over counts of paths", {'L', "multi-count-coverage"});
    args::ValueFlag<std::string> path_bedmulticov(parser, "BED", "for each BED entry, provide a table of path coverage over unique multisets of paths in the graph. Each unique multiset of paths overlapping a given BED interval is described in terms of its length relative to the total interval, the number of path traversals, and unique paths involved in these traversals.", {'B', "bed-multicov"});
    args::ValueFlag<std::string> path_delim(parser, "CHAR", "the part of each path name before this delimiter is a group identifier, which when specified will cause stats to be collected in a group-wise rather than path-wise fashion", {'D', "delim"});

    args::Flag mean_links_length(parser, "mean_links_length", "calculate the mean links length", {'l', "mean-links-length"});
    args::Flag dont_penalize_gap_links(parser, "dont-penalize-gap-links", "don't penalize gap links in the mean links length", {'g', "no-gap-links"});
    args::Flag sum_of_path_node_distances(parser, "sum_of_path_node_distances", "calculate the sum of path nodes distances", {'s', "sum-path-nodes-distances"});
    args::Flag penalize_diff_orientation(parser, "penalize_diff_orientation", "penalize links which connect nodes with different orientation", {'d', "penalize-different-orientation"});
    args::Flag path_statistics(parser, "path_statistics", "display the statistics for each path", {'P', "path-statistics"});

    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});

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

    if (!args::get(mean_links_length) && !args::get(sum_of_path_node_distances)){
        if (args::get(path_statistics)){
            std::cerr
                    << "[odgi stats] error: Please specify the -l/--mean-links-length and/or the -s/--sum-path-nodes-distances options to use the -P/--path-statistics option."
                    << std::endl;
            return 1;
        }

        if (layout_in_file){
            std::cerr
                    << "[odgi stats] error: Please specify the -l/--mean-links-length and/or the -s/--sum-path-nodes-distances options to specify the layout coordinates (-c/--coords-in option)."
                    << std::endl;
            return 1;
        }
    }

    if (args::get(penalize_diff_orientation) && !args::get(sum_of_path_node_distances)){
        std::cerr
                << "[odgi stats] error: Please specify the -s/--sum-path-nodes-distances option to use the -d/--penalize-different-orientation option."
                << std::endl;
        return 1;
    }
    if (args::get(dont_penalize_gap_links)){
        if (!args::get(mean_links_length)) {
            std::cerr
                    << "[odgi stats] error: Please specify the -l/--mean-links-length option to use the -g/--no-gap-links option."
                    << std::endl;
            return 1;
        } else if (layout_in_file){
            std::cerr
                    << "[odgi stats] error: The -g/--no-gap-links option can be specified together with the layout coordinates (-c/--coords-in option)."
                    << std::endl;
            return 1;
        }
    }

    if (args::get(threads)) {
        omp_set_num_threads(args::get(threads));
    } else {
        omp_set_num_threads(1);
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
    ///graph.display();

    if (args::get(summarize)) {
        uint64_t length_in_bp = 0, node_count = 0, edge_count = 0, path_count = 0;
        graph.for_each_handle([&](const handle_t& h) {
                length_in_bp += graph.get_length(h);
                ++node_count;
            });
        graph.for_each_edge([&](const edge_t& e) {
                ++edge_count;
                return true;
            });
        graph.for_each_path_handle([&](const path_handle_t& p) {
                ++path_count;
            });
        std::cout << "#length\tnodes\tedges\tpaths" << std::endl;
        std::cout << length_in_bp << "\t" << node_count << "\t" << edge_count << "\t" << path_count << std::endl;
    }
    if (args::get(base_content)) {
        std::vector<uint64_t> chars(256);
        graph.for_each_handle([&](const handle_t& h) {
                std::string seq = graph.get_sequence(h);
                for (auto c : seq) {
                    ++chars[c];
                }
            });
        for (uint64_t i = 0; i < 256; ++i) {
            if (chars[i]) {
                std::cout << (char)i << "\t" << chars[i] << std::endl;
            }
        }
    }
    if (args::get(path_coverage)) {
        std::map<uint64_t, uint64_t> full_histogram;
        std::map<uint64_t, uint64_t> unique_histogram;
        graph.for_each_handle([&](const handle_t& h) {
                std::vector<uint64_t> paths_here;
                graph.for_each_step_on_handle(h, [&](const step_handle_t& occ) {
                        paths_here.push_back(as_integer(graph.get_path(occ)));
                    });
                std::sort(paths_here.begin(), paths_here.end());
                std::vector<uint64_t> unique_paths = paths_here;
                unique_paths.erase(std::unique(unique_paths.begin(), unique_paths.end()), unique_paths.end());
                full_histogram[paths_here.size()] += graph.get_length(h);
                unique_histogram[unique_paths.size()] += graph.get_length(h);
            });
        std::cout << "type\tcov\tN" << std::endl;
        for (auto& p : full_histogram) {
            std::cout << "full\t" << p.first << "\t" << p.second << std::endl;
        }
        for (auto& p : unique_histogram) {
            std::cout << "uniq\t" << p.first << "\t" << p.second << std::endl;
        }
    }

    if (args::get(mean_links_length) || args::get(sum_of_path_node_distances)) {
        // This vector is needed for computing the metrics in 1D and for detecting gap-links
        std::vector<uint64_t> position_map(graph.get_node_count() + 1);

        // These vectors are needed for computing the metrics in 2D
        std::vector<double> X, Y;

        if (layout_in_file) {
            auto& infile = args::get(layout_in_file);
            if (infile.size()) {
                algorithms::layout::Layout layout;

                if (infile == "-") {
                    layout.load(std::cin);
                } else {
                    ifstream f(infile.c_str());
                    layout.load(f);
                    f.close();
                }

                X = layout.get_X();
                Y = layout.get_Y();
            }
        }

        uint64_t len = 0;
        nid_t last_node_id = graph.min_node_id();
        graph.for_each_handle([&](const handle_t &h) {
            nid_t node_id = graph.get_id(h);
            if (node_id - last_node_id > 1) {
                std::cerr << "[odgi stats] error: The graph is not optimized. Please run 'odgi sort' using -O, --optimize" << std::endl;
                exit(1);
            }
            last_node_id = node_id;

            position_map[number_bool_packing::unpack_number(h)] = len;

#ifdef debug_odgi_stats
            std::cerr << "SEGMENT ID: " << graph.get_id(h) << " - " << as_integer(h) << " - index_in_position_map (" << number_bool_packing::unpack_number(h) << ") = " << len << std::endl;
#endif

            uint64_t hl = graph.get_length(h);
            len += hl;
        });
        position_map[position_map.size() - 1] = len;

        if (args::get(mean_links_length)){
            bool _dont_penalize_gap_links = args::get(dont_penalize_gap_links);

            uint64_t sum_all_node_space = 0;
            uint64_t sum_all_nt_space = 0;
            double sum_all_2D_space = 0.0;
            uint64_t num_all_links = 0;
            uint64_t num_all_gap_links = 0;

            std::cout << "#mean_links_length" << std::endl;
            if (layout_in_file) {
                std::cout << "path\tin_2D_space\tnum_links_considered" << std::endl;
            }else{
                std::cout << "path\tin_node_space\tin_nucleotide_space\tnum_links_considered";

                if (dont_penalize_gap_links){
                    std::cout << "\tnum_gap_links_not_penalized" << std::endl;
                }else{
                    std::cout << std::endl;
                }
            }

            graph.for_each_path_handle([&](const path_handle_t &path) {
                std::set<uint64_t> ordered_unpacked_numbers_in_path; // ordered set to retain the positions order

                uint64_t sum_node_space = 0;
                uint64_t sum_nt_space = 0;
                double sum_2D_space = 0.0;
                uint64_t num_links = 0;
                uint64_t num_gap_links = 0;

                if (_dont_penalize_gap_links){
                    graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                        ordered_unpacked_numbers_in_path.insert(number_bool_packing::unpack_number(graph.get_handle_of_step(occ)));
                    });
                }

                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);

                    if (graph.has_next_step(occ)){
                        handle_t i = graph.get_handle_of_step(graph.get_next_step(occ));

                        uint64_t unpacked_h = number_bool_packing::unpack_number(h);
                        uint64_t unpacked_i = number_bool_packing::unpack_number(i);

                        // The position map includes the start and end of the node in successive entries.
                        // Edges leave from the end (or start) of one node (depending on whether they are on the
                        // forward or reverse strand) and go to the start (or end) of the other side of the link.
                        uint64_t _info_a = unpacked_h + !number_bool_packing::unpack_bit(h);
                        uint64_t _info_b = unpacked_i + number_bool_packing::unpack_bit(i);

                        if (!dont_penalize_gap_links || next(ordered_unpacked_numbers_in_path.find(unpacked_h)) != ordered_unpacked_numbers_in_path.find(unpacked_i)){
                            if (_info_b < _info_a){
                                _info_a = _info_b;
                                _info_b = unpacked_h + !number_bool_packing::unpack_bit(h);
                            }

                            if (layout_in_file) {
                                // 2D metric
                                double dx = X[2 * unpacked_h + number_bool_packing::unpack_number(h)] - X[2 * unpacked_i + number_bool_packing::unpack_bit(i)];
                                double dy = Y[2 * unpacked_h + number_bool_packing::unpack_number(h)] - Y[2 * unpacked_i + number_bool_packing::unpack_bit(i)];

                                sum_2D_space += sqrt(dx * dx + dy * dy);
                            }else{
                                // 1D metric (in node space and int nucleotide space)
                                sum_node_space += _info_b - _info_a;
                                sum_nt_space += position_map[_info_b] - position_map[_info_a];

#ifdef debug_odgi_stats
                                std::cerr << _info_b << " - " << _info_a << ": " << position_map[_info_b] - position_map[_info_a] << std::endl;
#endif
                            }
                        }else if (dont_penalize_gap_links){
                            num_gap_links++;
                        }

                        num_links++;
                    }
                });

                if (args::get(path_statistics)) {
                    double ratio_node_space = 0;
                    double ratio_nt_space = 0;
                    double ratio_2D_space = 0;
                    if (num_links > 0){
                        if (layout_in_file) {
                            ratio_2D_space = sum_2D_space / (double)num_links;
                        }else{
                            ratio_node_space = (double)sum_node_space / (double)num_links;
                            ratio_nt_space = (double)sum_nt_space / (double)num_links;
                        }
                    }

                    if (layout_in_file) {
                        std::cout << graph.get_path_name(path) << "\t" << ratio_2D_space << "\t" << num_links << std::endl;
                    }else{
                        std::cout << graph.get_path_name(path) << "\t" << ratio_node_space << "\t" << ratio_nt_space << "\t" << num_links;

                        if (dont_penalize_gap_links){
                            std::cout << "\t" << num_gap_links << std::endl;
                        }else{
                            std::cout << std::endl;
                        }
                    }
                }

                sum_all_node_space += sum_node_space;
                sum_all_nt_space += sum_nt_space;
                sum_all_2D_space += sum_2D_space;
                num_all_links += num_links;
                num_all_gap_links += num_gap_links;
            });

            double ratio_node_space = 0;
            double ratio_nt_space = 0;
            double ratio_2D_space = 0;
            if (num_all_links > 0) {
                if (layout_in_file) {
                    ratio_2D_space = sum_all_2D_space / (double)num_all_links;
                }else{
                    ratio_node_space = (double)sum_all_node_space / (double)num_all_links;
                    ratio_nt_space = (double)sum_all_nt_space / (double)num_all_links;
                }
            }

            if (layout_in_file) {
                std::cout << "all_paths\t" << ratio_2D_space << "\t" << num_all_links << std::endl;
            }else{
                std::cout << "all_paths\t" << ratio_node_space << "\t" << ratio_nt_space << "\t" << num_all_links;

                if (_dont_penalize_gap_links){
                    std::cout << "\t" << num_all_gap_links << std::endl;
                }else{
                    std::cout << std::endl;
                }
            }
        }

        if (args::get(sum_of_path_node_distances)){
            bool _penalize_diff_orientation = args::get(penalize_diff_orientation);

            uint64_t sum_all_path_node_dist_node_space = 0;
            uint64_t sum_all_path_node_dist_nt_space = 0;
            double sum_all_path_node_dist_2D_space = 0.0;
            uint64_t len_all_path_node_space = 0;
            uint64_t len_all_path_nt_space = 0;
            uint64_t num_all_penalties = 0;
            uint64_t num_all_penalties_diff_orientation = 0;

            std::cout << "#sum_of_path_node_distances" << std::endl;

            if (layout_in_file) {
                std::cout << "path\tin_2D_space_by_nodes\tin_2D_space_by_nucleotides\tnodes\tnucleotides";
            }else{
                std::cout << "path\tin_node_space\tin_nucleotide_space\tnodes\tnucleotides\tnum_penalties";
            }

            if (_penalize_diff_orientation){
                std::cout << "\tnum_penalties_different_orientation" << std::endl;
            }else{
                std::cout << std::endl;
            }

            graph.for_each_path_handle([&](const path_handle_t &path) {
#ifdef debug_odgi_stats
                std::cerr << "path_name: " << graph.get_path_name(path) << std::endl;
#endif
                uint64_t sum_path_node_dist_node_space = 0;
                uint64_t sum_path_node_dist_nt_space = 0;
                double sum_path_node_dist_2D_space = 0.0;
                uint64_t len_path_node_space = 0;
                uint64_t len_path_nt_space = 0;
                uint64_t num_penalties = 0;
                uint64_t num_penalties_diff_orientation = 0;

                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);

                    if (graph.has_next_step(occ)){
                        handle_t i = graph.get_handle_of_step(graph.get_next_step(occ));

                        uint64_t unpacked_a = number_bool_packing::unpack_number(h);
                        uint64_t unpacked_b = number_bool_packing::unpack_number(i);

                        double euclidean_distance_2D;

                        if (layout_in_file) {
                            // 2D metric
                            double dx = X[2 * unpacked_a + number_bool_packing::unpack_number(h)] - X[2 * unpacked_b + number_bool_packing::unpack_bit(i)];
                            double dy = Y[2 * unpacked_a + number_bool_packing::unpack_number(h)] - Y[2 * unpacked_b + number_bool_packing::unpack_bit(i)];

                            euclidean_distance_2D = sqrt(dx * dx + dy * dy);
                            sum_path_node_dist_2D_space += euclidean_distance_2D;
                        }else{
                            uint8_t weight = 1;
                            if (unpacked_b < unpacked_a){
                                unpacked_a = unpacked_b;
                                unpacked_b = number_bool_packing::unpack_number(h);

                                // When a path goes back in terms of pangenomic order, this is punished
                                weight = 3;
                                num_penalties++;
                            }

#ifdef debug_odgi_stats
                            std::cerr << unpacked_b << " - " << unpacked_a << ": " << position_map[unpacked_b] - position_map[unpacked_a] << " * " << weight << std::endl;
#endif

                            sum_path_node_dist_node_space += weight * (unpacked_b - unpacked_a);
                            sum_path_node_dist_nt_space += weight * (position_map[unpacked_b] - position_map[unpacked_a]);
                        }

                        if (_penalize_diff_orientation && (number_bool_packing::unpack_bit(h) != number_bool_packing::unpack_bit(i))){
                            if (layout_in_file) {
                                sum_path_node_dist_2D_space += 2 * euclidean_distance_2D;
                            }else{
                                sum_path_node_dist_node_space += 2 * (unpacked_b - unpacked_a);
                                sum_path_node_dist_nt_space += 2 * (position_map[unpacked_b] - position_map[unpacked_a]);
                            }

                            num_penalties_diff_orientation++;
                        }
                    }

                    len_path_node_space++;
                    len_path_nt_space += graph.get_length(h);
                });

                if (args::get(path_statistics)) {
                    if (layout_in_file) {
                        std::cout << graph.get_path_name(path) << "\t" << (double)sum_path_node_dist_2D_space / (double)len_path_node_space << "\t" << (double)sum_path_node_dist_2D_space / (double)len_path_nt_space << "\t" << len_path_node_space << "\t" << len_path_nt_space;
                    }else{
                        std::cout << graph.get_path_name(path) << "\t" << (double)sum_path_node_dist_node_space / (double)len_path_node_space << "\t" << (double)sum_path_node_dist_nt_space / (double)len_path_nt_space << "\t" << len_path_node_space << "\t" << len_path_nt_space  << "\t" << num_penalties;
                    }

                    if (_penalize_diff_orientation){
                        std::cout << "\t" << num_penalties_diff_orientation << std::endl;
                    }else{
                        std::cout << std::endl;
                    }
                }

                sum_all_path_node_dist_node_space += sum_path_node_dist_node_space;
                sum_all_path_node_dist_nt_space += sum_path_node_dist_nt_space;
                sum_all_path_node_dist_2D_space += sum_path_node_dist_2D_space;
                len_all_path_node_space += len_path_node_space;
                len_all_path_nt_space += len_path_nt_space;
                num_all_penalties += num_penalties;
                num_all_penalties_diff_orientation += num_penalties_diff_orientation;
            });

            if (layout_in_file) {
                std::cout << "all_paths\t" << (double)sum_all_path_node_dist_2D_space / (double)len_all_path_node_space << "\t" << (double)sum_all_path_node_dist_2D_space / (double)len_all_path_nt_space << "\t" << len_all_path_node_space << "\t" << len_all_path_nt_space;
            }else{
                std::cout << "all_paths\t" << (double)sum_all_path_node_dist_node_space / (double)len_all_path_node_space << "\t" << (double)sum_all_path_node_dist_nt_space / (double)len_all_path_nt_space << "\t" << len_all_path_node_space << "\t" << len_all_path_nt_space << "\t" << num_all_penalties;
            }

            if (_penalize_diff_orientation){
                std::cout << "\t" << num_all_penalties_diff_orientation << std::endl;
            }else{
                std::cout << std::endl;
            }
        }
    }

    bool using_delim = !args::get(path_delim).empty();
    char delim = '\0';
    ska::flat_hash_map<std::string, uint64_t> path_group_ids;
    std::vector<std::string> path_groups;
    if (using_delim) {
        delim = args::get(path_delim).at(0);
        uint64_t i = 0;
        graph.for_each_path_handle(
            [&](const path_handle_t& p) {
                std::string group_name = split(graph.get_path_name(p), delim)[0];
                auto f = path_group_ids.find(group_name);
                if (f == path_group_ids.end()) {
                    path_group_ids[group_name] = i++;
                    path_groups.push_back(group_name);
                }
            });
    }

    auto get_path_name
        = (using_delim ?
           (std::function<std::string(const uint64_t&)>)
           [&](const uint64_t& id) { return path_groups[id]; }
           :
           (std::function<std::string(const uint64_t&)>)
           [&](const uint64_t& id) { return graph.get_path_name(as_path_handle(id)); });

    auto get_path_id
        = (using_delim ?
           (std::function<uint64_t(const path_handle_t&)>)
           [&](const path_handle_t& p) {
               std::string group_name = split(graph.get_path_name(p), delim)[0];
               return path_group_ids[group_name];
           }
           :
           (std::function<uint64_t(const path_handle_t&)>)
           [&](const path_handle_t& p) {
               return as_integer(p);
           });

    if (args::get(path_setcov)) {
        uint64_t total_length = 0;
        std::map<std::set<uint64_t>, uint64_t> setcov;
        graph.for_each_handle(
            [&](const handle_t& h) {
                std::set<uint64_t> paths_here;
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t& occ) {
                        paths_here.insert(get_path_id(graph.get_path(occ)));
                    });
                size_t l = graph.get_length(h);
#pragma omp critical (setcov)
                {
                    setcov[paths_here] += l;
                    total_length += l;
                }
            }, true);
        std::cout << "length\tgraph.frac\tn.paths\tpath.set" << std::endl;
        for (auto& p : setcov) {
            std::cout << p.second << "\t"
                      << (double)p.second/(double)total_length << "\t"
                      << p.first.size() << "\t";
            for (auto& i : p.first) {
                std::cout << get_path_name(i) << ",";
            }
            std::cout << std::endl;
        }
    }

    if (args::get(path_setcov_count)) {
        uint64_t total_length = 0;
        std::map<uint64_t, uint64_t> setcov_count;
        graph.for_each_handle(
            [&](const handle_t& h) {
                std::set<uint64_t> paths_here;
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t& occ) {
                        paths_here.insert(get_path_id(graph.get_path(occ)));
                    });
                size_t l = graph.get_length(h);
#pragma omp critical (setcov)
                {
                    setcov_count[paths_here.size()] += l;
                    total_length += l;
                }
            }, true);
        std::cout << "length\tgraph.frac\tunique.path.count" << std::endl;
        for (auto& p : setcov_count) {
            std::cout << p.second << "\t"
                      << (double)p.second/(double)total_length
                      << "\t" << p.first << std::endl;
        }
    }

    if (args::get(path_multicov)) {
        uint64_t total_length = 0;
        std::map<std::vector<uint64_t>, uint64_t> multisetcov;
        graph.for_each_handle(
            [&](const handle_t& h) {
                std::vector<uint64_t> paths_here;
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t& occ) {
                        paths_here.push_back(get_path_id(graph.get_path(occ)));
                    });
                std::sort(paths_here.begin(), paths_here.end());
                size_t l = graph.get_length(h);
#pragma omp critical (multisetcov)
                {
                    multisetcov[paths_here] += l;
                    total_length += l;
                }
            }, true);
        std::cout << "length\tgraph.frac\nn.paths\tpath.multiset" << std::endl;
        for (auto& p : multisetcov) {
            std::cout << p.second << "\t"
                      << (double)p.second/(double)total_length << "\t"
                      << p.first.size() << "\t";
            bool first = true;
            for (auto& i : p.first) {
                std::cout << (first ? (first=false, "") :",") << get_path_name(i);
            }
            std::cout << std::endl;
        }
    }

    if (args::get(path_multicov_count)) {
        uint64_t total_length = 0;
        std::map<uint64_t, uint64_t> multisetcov_count;
        graph.for_each_handle(
            [&](const handle_t& h) {
                uint64_t paths_here = 0;
                graph.for_each_step_on_handle(
                    h,
                    [&](const step_handle_t& occ) {
                        ++paths_here;
                    });
                size_t l = graph.get_length(h);
#pragma omp critical (multisetcov)
                {
                    multisetcov_count[paths_here] += l;
                    total_length += l;
                }
            }, true);
        std::cout << "length\tgraph.frac\tpath.step.count" << std::endl;
        for (auto& p : multisetcov_count) {
            std::cout << p.second << "\t"
                      << (double)p.second/(double)total_length << "\t"
                      << p.first << std::endl;
        }
    }

    if (!args::get(path_bedmulticov).empty()) {
        std::string line;
        typedef std::map<std::vector<uint64_t>, uint64_t> setcov_t;
        typedef IntervalTree<uint64_t, std::pair<std::string, setcov_t*> > itree_t;
        map<std::string, itree_t::interval_vector> intervals;
        auto& x = args::get(path_bedmulticov);
        std::ifstream bed_in(x);
        while (std::getline(bed_in, line)) {
            // BED is base-numbered, 0-origin, half-open.  This parse turns that
            // into base-numbered, 0-origin, fully-closed for internal use.  All
            // coordinates used internally should be in the latter, and coordinates
            // from the user in the former should be converted immediately to the
            // internal format.
            std::vector<string> fields = split(line, '\t');
            intervals[fields[0]].push_back(
                itree_t::interval(
                    std::stoul(fields[1]),
                    std::stoul(fields[2]),
                    make_pair(fields[3], new setcov_t())));
        }

        std::vector<std::string> path_names;
        path_names.reserve(intervals.size());
        for(auto const& i : intervals) {
            path_names.push_back(i.first);
        }

        // the header
        std::cout << "path.name" << "\t"
                  << "bed.name" << "\t"
                  << "bed.start" << "\t"
                  << "bed.stop" << "\t"
                  << "bed.len" << "\t"
                  << "path.set.state.len" << "\t"
                  << "path.set.frac" << "\t"
                  << "path.traversals" << "\t"
                  << "uniq.paths.in.state" << "\t"
                  << "path.multiset" << std::endl;

#pragma omp parallel for
        for (uint64_t k = 0; k < path_names.size(); ++k) {
            auto& path_name = path_names.at(k);
            auto& path_ivals = intervals[path_name];
            path_handle_t path = graph.get_path_handle(path_name);
            // build the intervals for each path we'll query
            itree_t itree(std::move(path_ivals)); //, 16, 1);
            uint64_t pos = 0;
            graph.for_each_step_in_path(graph.get_path_handle(x), [&](const step_handle_t& occ) {
                    std::vector<uint64_t> paths_here;
                    handle_t h = graph.get_handle_of_step(occ);
                    graph.for_each_step_on_handle(
                        h,
                        [&](const step_handle_t& occ) {
                            paths_here.push_back(get_path_id(graph.get_path(occ)));
                        });
                    uint64_t len = graph.get_length(h);
                    // check each position in the node
                    auto hits = itree.findOverlapping(pos, pos+len);
                    if (hits.size()) {
                        for (auto& h : hits) {
                            auto& q = *h.value.second;
                            // adjust length for overlap length
                            uint64_t ovlp = len - (h.start > pos ? h.start - pos : 0) - (h.stop < pos+len ? pos+len - h.stop : 0);
                            q[paths_here] += ovlp;
                        }
                    }
                    pos += len;
                });
            itree.visit_all([&](const itree_t::interval& ival) {
                    auto& name = ival.value.first;
                    auto& setcov = *ival.value.second;
                    for (auto& p : setcov) {
                        std::set<uint64_t> u(p.first.begin(), p.first.end());
#pragma omp critical (cout)
                        {
                            std::cout << path_name << "\t" << name << "\t" << ival.start << "\t" << ival.stop << "\t"
                                      << ival.stop - ival.start << "\t"
                                      << p.second << "\t"
                                      << (float)p.second/(ival.stop-ival.start) << "\t"
                                      << p.first.size() << "\t"
                                      << u.size() << "\t";
                            bool first = true;
                            for (auto& i : p.first) {
                                std::cout << (first ? (first=false, "") :",") << get_path_name(i);
                            }
                            std::cout << std::endl;
                        }
                    }
            });
            for (auto& ival : path_ivals) {
                delete ival.value.second;
            }
        }
    }

    return 0;
}

static Subcommand odgi_stats("stats", "describe the graph and its path relationships",
                              PIPELINE, 3, main_stats);


}
