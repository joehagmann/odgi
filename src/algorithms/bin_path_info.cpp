#include "bin_path_info.hpp"

// #define  debug_bin_path_info

namespace odgi {
    namespace algorithms {

        template<int first, int second>
        void updatePair(std::pair<uint64_t, uint64_t> &p, uint64_t nucleotide_count) {
#ifdef debug_bin_path_info
            std::cerr << "UPDATE_PAIR: " << "<" << std::get<first>(p) << "," << std::get<second>(p) << ">" << std::endl;
            std::cerr << "NUC: " << nucleotide_count << std::endl;
#endif
            if (std::get<second>(p) == 0) {
                std::get<second>(p) = std::get<first>(p);

            }
            std::get<first>(p) = nucleotide_count;
#ifdef debug_bin_path_info
            std::cerr << "AFTER_UPDATE_PAIR: " << "<" << std::get<first>(p) << "," << std::get<second>(p) << ">"
                      << std::endl;
#endif
        }

        void bin_path_info(const PathHandleGraph &graph,
                           const std::string &prefix_delimiter,
                           const std::function<void(const uint64_t &, const uint64_t &)> &handle_header,
                           const std::function<void(const std::string &,
                                                    const std::vector<std::pair<uint64_t, uint64_t>> &,
                                                    const std::map<uint64_t, algorithms::path_info_t> &)> &handle_path,
                           const std::function<void(const uint64_t &, const std::string &)> &handle_sequence,
                           const std::function<void(const string &)> &handle_fasta,
                           uint64_t num_bins,
                           uint64_t bin_width,
                           bool drop_gap_links) {
            // the graph must be compacted for this to work
            std::vector<uint64_t> position_map(graph.get_node_count() + 1);
            uint64_t len = 0;
            std::string graph_seq;
            graph.for_each_handle([&](const handle_t &h) {
                position_map[number_bool_packing::unpack_number(h)] = len;
                uint64_t hl = graph.get_length(h);
                graph_seq.append(graph.get_sequence(h));
                len += hl;
            });
            if (!num_bins) {
                num_bins = len / bin_width + (len % bin_width ? 1 : 0);
            } else if (!bin_width) {
                bin_width = len / num_bins;
                num_bins = len / bin_width + (len % bin_width ? 1 : 0);
            }
            position_map[position_map.size() - 1] = len;
            // write header
            handle_header(len, bin_width);
            // collect bin sequences
            for (uint64_t i = 0; i < num_bins; ++i) {
                handle_sequence(i + 1, graph_seq.substr(i * bin_width, bin_width));
            }
            // write out pangenome sequence if wished so
            handle_fasta(graph_seq);
            graph_seq.clear(); // clean up
            std::unordered_map<path_handle_t, uint64_t> path_length;
            uint64_t gap_links_removed = 0;
            uint64_t total_links = 0;
            graph.for_each_path_handle([&](const path_handle_t &path) {
                std::vector<std::pair<uint64_t, uint64_t>> links;
                std::map<uint64_t, path_info_t> bins;
                // walk the path and aggregate
                uint64_t path_pos = 0;
                int64_t last_bin = 0; // flag meaning "null bin"
                uint64_t last_pos_in_bin = 0;
                uint64_t nucleotide_count = 0;
                bool last_is_rev = false;
                graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                    handle_t h = graph.get_handle_of_step(occ);
                    bool is_rev = graph.get_is_reverse(h);
                    uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                    uint64_t hl = graph.get_length(h);
                    // detect bin crossings
                    // make contects for the bases in the node
                    for (uint64_t k = 0; k < hl; ++k) {
                        int64_t curr_bin = (p + k) / bin_width + 1;
                        uint64_t curr_pos_in_bin = (p + k) - (curr_bin * bin_width);
                        if (curr_bin != last_bin && std::abs(curr_bin - last_bin) > 1 || last_bin == 0) {
                            // bin cross!
                            links.push_back(std::make_pair(last_bin, curr_bin));
                        }
                        ++bins[curr_bin].mean_cov;
                        if (is_rev) {
                            ++bins[curr_bin].mean_inv;
                        }
                        bins[curr_bin].mean_pos += path_pos++;
                        nucleotide_count += 1;
                        if ((bins[curr_bin].ranges.size() == 0) ||
                            ((nucleotide_count - bins[curr_bin].ranges.back().second) > 1 &&
                             (nucleotide_count - bins[curr_bin].ranges.back().first) > 1) ||
                            (is_rev != last_is_rev)) {
                            std::pair<uint64_t, uint64_t> p = std::make_pair(0, 0);
                            if (is_rev) {
                                std::get<0>(p) = nucleotide_count;
                            } else {
                                std::get<1>(p) = nucleotide_count;
                            }
                            bins[curr_bin].ranges.push_back(p);
#ifdef debug_bin_path_info
                            std::cerr << "PUSHED PAIR: " << "<" << std::get<0>(p) << "," << std::get<1>(p) << ">"
                                      << std::endl;
#endif
                        } else {
                            std::pair<uint64_t, uint64_t> &p = bins[curr_bin].ranges.back();
                            if (is_rev) {
                                updatePair<0, 1>(p, nucleotide_count);
                            }
                            else {
                                updatePair<1, 0>(p, nucleotide_count);
                            }
                        }
                        last_bin = curr_bin;
                        last_is_rev = is_rev;
                        last_pos_in_bin = curr_pos_in_bin;
                    }
                });
                links.push_back(std::make_pair(last_bin, 0));
                uint64_t path_length = path_pos;
                uint64_t end_nucleotide = nucleotide_count;
                for (auto &entry : bins) {
                    auto &v = entry.second;
                    v.mean_inv /= (v.mean_cov ? v.mean_cov : 1);
                    v.mean_cov /= bin_width;
                    v.mean_pos /= bin_width * path_length * v.mean_cov;
                }

                if (drop_gap_links) {
                    std::vector<uint64_t> bin_ids;
                    for (const auto &entry: bins) {
                        bin_ids.push_back(entry.first);
                    }
                    std::sort(bin_ids.begin(), bin_ids.end());
                    total_links += links.size();

                    uint64_t fill_pos = 0;

                    for (uint64_t i = 0; i < links.size(); ++i) {
                        auto link = links[i];

                        if (link.first == 0 || link.second == 0)
                            continue;

                        if (link.first > link.second) {
                            links[fill_pos++] = link;
                            continue;
                        }

                        auto left_it = std::lower_bound(bin_ids.begin(), bin_ids.end(), link.first + 1);
                        auto right_it = std::lower_bound(bin_ids.begin(), bin_ids.end(), link.second);
                        if (right_it > left_it) {
                            links[fill_pos++] = link;
                        }
                    }

                    gap_links_removed += links.size() - fill_pos;
                    links.resize(fill_pos);
                }

                handle_path(graph.get_path_name(path), links, bins);
            });

            if (drop_gap_links) {
                uint64_t path_count = graph.get_path_count();

                std::cerr << std::setprecision(4) << "Gap links removed: " << (100.0 *  ((double)gap_links_removed / (double)total_links))
                          << "%, that is " << gap_links_removed << " gap links (" << path_count << " path start links + "
                          << path_count << " path end links + " << (gap_links_removed - path_count * 2) << " inner gap links) of "
                          << total_links << " total links" << std::endl;
            }
        }
    

        int get_count_of_link_in_bin(const std::vector<link_info_t> &links, const link_info_t &link) {
            int count=0;
            for (const auto &l : links) {
                if (l.from == link.from && l.to == link.to) count++;
            }
            return(count);
        }

        void print_link(link_info_t l) {
            std::cout << "  link from " << l.from << " to " << l.to << ", " <<
              l.times << "x, " << l.nr_passes << " passes" << std::endl;
        }

        bool is_not_adjacent_link(link_info_t link, uint64_t adj_link) {
            if (link.from == adj_link || link.to == adj_link) return(false);
            else return(true);
        }


        void bin_path_info_for_pantograph(const PathHandleGraph &graph,
                            const std::string &prefix_delimiter,
                            const std::function<void(const uint64_t &, const uint64_t &)> &handle_header,
                            const std::function<void(const std::string &,
                                                     const std::map<uint64_t, algorithms::bin_info_t> &,
                                                     const bool)> &handle_path,
                            const std::function<void(const uint64_t &, const std::string &)> &handle_sequence,
                            const std::function<void(const string &)> &handle_fasta,
                            const std::function<void(const uint64_t &, const uint64_t &, const uint64_t &)> &handle_xoffset,
                            uint64_t num_bins,
                            uint64_t bin_width) {
                // the graph must be compacted for this to work

                // Containers:
                // position_map: for each node stores the pangenome position
                // link_columns
                // bin_pfreq: bin presence/path frequency

                std::vector<uint64_t> position_map(graph.get_node_count() + 1);
                uint64_t len = 0;
                std::string graph_seq;
                graph.for_each_handle([&](const handle_t &h) {
                    position_map[number_bool_packing::unpack_number(h)] = len;
                    uint64_t hl = graph.get_length(h);
                    graph_seq.append(graph.get_sequence(h));
                    len += hl;
                });
                if (!num_bins) {
                    num_bins = len / bin_width + (len % bin_width ? 1 : 0);
                } else if (!bin_width) {
                    bin_width = len / num_bins;
                    num_bins = len / bin_width + (len % bin_width ? 1 : 0);
                }
                position_map[position_map.size() - 1] = len;
                // write header
                handle_header(len, bin_width);

                // deprecated? ///////////////////////////
                    std::vector<std::unordered_set<std::string>> bin_pfreq(num_bins);
                    //std::vector<std::vector<std::pair<uint64_t, uint64_t>>> links_per_bin(num_bins);
                    std::vector<std::unordered_map<std::string, std::vector<int64_t>>> links_per_bin(num_bins);

                    for (uint64_t i = 0; i < num_bins; ++i) {
                        std::unordered_set<std::string> s;
                        bin_pfreq.push_back(s);
                        std::unordered_map<std::string, std::vector<int64_t>> m;
                        links_per_bin.push_back(m);
                    }
                    std::unordered_map<uint64_t, uint64_t> link_columns;
                /////////////////////////////////////////

                // !!!! LATEST:
                //std::unordered_map<std::string, std::unordered_map< int64_t, std::vector<link_info_t> >> links_per_path_n_bin; // per path!
                std::unordered_map<int64_t, std::unordered_set<std::string>> uniq_links;
                std::unordered_map<std::string, int> links_nr_passes;
                
                //std::unordered_map<path_handle_t, uint64_t> path_length;
                uint64_t gap_links_removed = 0;
                uint64_t total_links = 0;
                bool first_path = true;
                graph.for_each_path_handle([&](const path_handle_t &path) {
                    std::string path_name = graph.get_path_name(path);
                    std::cout << "path: " << path_name << std::endl;
                    //std::vector<std::pair<uint64_t, uint64_t>> links;
                    //std::vector<link_info_t> links;
                    //std::vector<int> links_nr_passes;
                    std::map<uint64_t, bin_info_t> bins;
                    std::map<uint64_t, bool> bins_revisited; // TODO: I think not needed, remove all occs

                    // walk the path and aggregate
                    uint64_t path_pos = 0;
                    int64_t last_bin = 0; // flag meaning "null bin"
                    uint64_t last_pos_in_bin = 0;
                    uint64_t nucleotide_count = 0;
                    bool last_is_rev = false;
                    graph.for_each_step_in_path(path, [&](const step_handle_t &occ) {
                        handle_t h = graph.get_handle_of_step(occ);
                        bool is_rev = graph.get_is_reverse(h);
                        uint64_t p = position_map[number_bool_packing::unpack_number(h)];
                        uint64_t hl = graph.get_length(h);

                        // fill each bin
                        for (uint64_t k = 0; k < hl; ++k) {
                            int64_t curr_bin = (p + k) / bin_width + 1;
                            uint64_t curr_pos_in_bin = (p + k) - (curr_bin * bin_width);

                            // detect non-consecutive bin ordering for the path
                            std::pair<uint64_t, uint64_t> p = std::make_pair(last_bin, curr_bin);
                            std::vector<int64_t> v = {last_bin, curr_bin, 1, 1}; // from, to, times, nr_passes
                            // check if previously inserted link is the same as this one.
                            //   If so, increase variable times
                            //   If not, insert a new link, even if it was already traversed in an 
                            //     earlier pass through the bin
                            //if (links_per_bin[curr_bin].size() > 0 && curr_bin != last_bin) {
                            //if (links_per_bin[curr_bin].find(path_name) != links_per_bin[curr_bin].end()) {
                                //if (links_per_bin[curr_bin][path_name].back() == v) {
                            if (curr_bin != last_bin && last_bin != 0) {
                            //if (bins[curr_bin].links.size() > 0 && curr_bin != last_bin) {
                                if (bins[curr_bin].links_to.size() > 0 && 
                                    bins[curr_bin].links_to.back().from == last_bin) { //&&
                                    //bins[curr_bin].links.back().to   == curr_bin) {
                                        // increase the number of 'times' the link is traversed consecutively
                                        ++bins[curr_bin].links_to.back().times;
                                } else {
                                    link_info_t l = { last_bin, curr_bin, 1, 1 };
                                    bins[curr_bin].links_to.push_back(l);
                                }
                                
                                if (bins[last_bin].links_from.size() > 0 && 
                                    bins[last_bin].links_from.back().to == curr_bin) {
                                        // increase the number of 'times' the link is traversed consecutively
                                        ++bins[last_bin].links_from.back().times; // TODO: last link at bin last_bin must not be this link! MAJOR ISSUE, we have to search that link first
                                        // TODO: storing the links once (as pointers) and pointing to them would save space, esp. when I store all neighboring links temporarily, then the search for the link in bin[last_bin] would be obsolete ------ might be solved by the check for bin[last_bin].links_from.back() ?
                                        // TODO #1: alternative: log nr. visits of bin for that path, and if nr. passes is increased, check if sum of nr.passes is nr.visits. if not, insert link to neighbor
                                } else {
                                    link_info_t l = { last_bin, curr_bin, 1, 1 };
                                    bins[last_bin].links_from.push_back(l);
                                }
                        

                                // bins_revisited.insert(std::make_pair(curr_bin, false));
                                // int nr_passes = get_count_of_link_in_bin(links, l);
                                // link_t link = {last_bin, curr_bin};
                                // bins[curr_bin].links_nr_passes[link] = nr_passes;

                                // if (nr_passes > 1) bins_revisited[curr_bin] = true;
                                // bins[curr_bin].links.back().nr_passes = nr_passes;
                                //      bins[curr_bin].links_nr_passes[links.size()-1] = get_count_of_link_in_bin(links, l);

                                // std::string l_str = to_string(last_bin) + "_" + to_string(curr_bin);
                                // if (links_nr_passes.find(l_str) != links_nr_passes.end()) links_nr_passes[l_str] = 0;
                                // links_nr_passes[l_str]++;

                                //bins[curr_bin].last_link = l;

                                // TODO insert links at bin struct of other end of link as well!!!!
                                
                            }
                            // else if (curr_bin != last_bin && std::abs(curr_bin - last_bin) > 1 || last_bin == 0) {
                            //     // make a link
                            //     links.push_back(std::make_pair(last_bin, curr_bin));
                            //     std::pair<uint64_t, uint64_t> p = std::make_pair(last_bin, curr_bin);
                            //     if (links_per_bin[curr_bin].find(p) != links_per_bin[curr_bin].end()) {
                                    
                            //     }

                            //     links_per_bin[curr_bin].push_back(std::make_pair(last_bin, curr_bin));

                            //     // link has to be generated independent of curr_bin==last_bin+1, and further below it might get removed if there is no real link!
                            //     // for cases where node is traversed multiple times to tell the different passes apart
                                
                            // }
                            ++bins[curr_bin].mean_cov;
                            if (is_rev) {
                                ++bins[curr_bin].mean_inv;
                            }
                            bins[curr_bin].mean_pos += path_pos++;
                            // fill bin_maf
                            bin_pfreq[curr_bin-1].insert(graph.get_path_name(path));
//std::cout << curr_bin << ": " << bin_pfreq[curr_bin-1].size() << std::endl;

                            // fill bin ranges
                            nucleotide_count += 1;
                            if ((bins[curr_bin].ranges.size() == 0) ||
                                ((nucleotide_count - bins[curr_bin].ranges.back().second) > 1 &&
                                (nucleotide_count - bins[curr_bin].ranges.back().first) > 1) ||
                                (is_rev != last_is_rev)) {
                                std::pair<uint64_t, uint64_t> p = std::make_pair(0, 0);
                                if (is_rev) {
                                    std::get<0>(p) = nucleotide_count;
                                } else {
                                    std::get<1>(p) = nucleotide_count;
                                }
                                bins[curr_bin].ranges.push_back(p);
    #ifdef debug_bin_path_info
                                std::cerr << "PUSHED PAIR: " << "<" << std::get<0>(p) << "," << std::get<1>(p) << ">"
                                        << std::endl;
    #endif
                            } else {
                                std::pair<uint64_t, uint64_t> &p = bins[curr_bin].ranges.back();
                                if (is_rev) {
                                    updatePair<0, 1>(p, nucleotide_count);
                                }
                                else {
                                    updatePair<1, 0>(p, nucleotide_count);
                                }
                            }
                            last_bin = curr_bin;
                            last_is_rev = is_rev;
                            last_pos_in_bin = curr_pos_in_bin;
                        }
                    });
                    //links.push_back(std::make_pair(last_bin, 0));
                    
                    // calculate mean values for each bin
                    // uint64_t path_length = path_pos;
                    // uint64_t end_nucleotide = nucleotide_count;
                    // for (auto &entry : bins) {
                    //     auto &v = entry.second;
                    //     v.mean_inv /= (v.mean_cov ? v.mean_cov : 1);
                    //     v.mean_cov /= bin_width;
                    //     v.mean_pos /= bin_width * path_length * v.mean_cov;
                    // }

                    // drop gap links:
                    // std::vector<uint64_t> bin_ids;
                    // for (const auto &entry: bins) {
                    //     bin_ids.push_back(entry.first);
                    // }
                    // std::sort(bin_ids.begin(), bin_ids.end());
                    //total_links += links.size();
                    // TODO merge previous for loop with the one before

                    uint64_t path_length = path_pos;
                    uint64_t fill_pos = 0;
                    std::map<uint64_t, bin_info_t>::iterator it = bins.begin();
                    for (it; it!=bins.end(); ++it) {
                        uint64_t bin_id = it->first;
                        auto &v = it->second;
                        v.mean_inv /= (v.mean_cov ? v.mean_cov : 1);
                        v.mean_cov /= bin_width;
                        v.mean_pos /= bin_width * path_length * v.mean_cov;
std::cout << "bin " << bin_id << std::endl;

    // TODO do the following for links_to and links_from...
                        std::vector<link_info_t> links;
                        int which_links = 1;
                        //for (int which_links=0; which_links<=1; ++which_links) {
                            if (which_links==0) links = bins[bin_id].links_to;
                            else                links = bins[bin_id].links_from;

std::cout << "  nr links: " << links.size() << std::endl;
                            if (links.size() == 0) continue;
                            for (std::vector<link_info_t>::iterator it1 = links.begin(); it1 != links.end() - 1; ++it1) {
                                for (std::vector<link_info_t>::iterator it2 = std::next(it1, 1); it2 != links.end(); ++it2) {
                                    if (it1->from == it2->from && it1->to == it2->to) {
                                        it1->nr_passes++;
                                        it2->nr_passes++;
std::cout << "nr passes incr" << std::endl;
                                    }
                                }
                            }
                        
                            // if (it == bins.begin()) continue; // first bin 
                            // std::map<uint64_t, bin_info_t>::iterator prev_it = std::prev(it, 1);

                            // std::map<uint64_t, bin_info_t>::iterator next_it = std::next(it, 1);
    if (it != bins.begin()) std::cout << "  previous bin: " << std::prev(it, 1)->first << std::endl;
    if (it != std::prev(bins.end(), 1)) std::cout << "  next bin: " << std::next(it, 1)->first << std::endl;
                            uint64_t fill_pos = 0;
                            for (uint64_t i=0; i<links.size(); ++i) {
    print_link(links[i]);
                                //should the link be retained? remove all links connecting adjacent bins
                                bool retain_link = false;
                                if (it != bins.begin() && it != std::prev(bins.end()) &&
                                        is_not_adjacent_link(links[i], std::prev(it, 1)->first) &&
                                        is_not_adjacent_link(links[i], std::next(it, 1)->first))
                                    retain_link = true;
                                if (it == bins.begin() &&
                                        is_not_adjacent_link(links[i], std::next(it, 1)->first))
                                    retain_link = true;
                                if (it == std::prev(bins.end()) && 
                                        is_not_adjacent_link(links[i], std::prev(it, 1)->first))
                                    retain_link = true;
                                if (links[i].nr_passes > 1)
                                    retain_link = true;

                                // TODO: don't check links_to ???

                                if (retain_link) {

                                // bool not_right_adj = (it != bins.begin() && is_not_adjacent_link(links[i], std::prev(it, 1)->first));
                                // if ((it != bins.begin() && is_not_adjacent_link(links[i], std::prev(it, 1)->first)) ||
                                //     (it != std::prev(bins.end()) && is_not_adjacent_link(links[i], std::next(it, 1)->first)) ||
                                //     //     links[i].to   != std::prev(it, 1)->first && 
                                //     //     links[i].from != std::prev(it, 1)->first) ||
                                //     //  (it != std::prev(bins.end()) &&
                                //     //     links[i].to   != std::next(it, 1)->first && 
                                //     //     links[i].from != std::next(it, 1)->first) ||
                                //      links[i].nr_passes > 1) {
    std::cout << "    condition true" << std::endl;
                                    //bins[bin_id].links_nr_passes[linkpair] > 1) {
                                        if (which_links == 0) {
                                            bins[bin_id].links_to[fill_pos] = bins[bin_id].links_to[i];
                                        } else {
                                            bins[bin_id].links_from[fill_pos] = bins[bin_id].links_from[i];
                                        }
                                        fill_pos++;

                                        link_columns[links[i].from]++;
                                        link_columns[links[i].to]++;
                                } else {
                                    // TODO: this is complicated, we here remove the link from the other bin, if it's still in there. Could likely be deleted if TODO #1 above is implemented
                                    if (which_links == 1) {
                                        uint64_t other_bin = links[i].to;
                                        uint64_t fill = 0;
                                        for (size_t j=0; j<bins[other_bin].links_to.size(); ++j) {
                                            link_info_t other_link = bins[other_bin].links_to[j];
                                            if (other_link.from != links[i].from) {
                                                bins[other_bin].links_to[fill] = other_link;
                                                fill++;
                                            }
                                        }
                                        if (other_bin > bin_id) {
                                            total_links += bins[other_bin].links_to.size();
                                        }
                                        gap_links_removed +=  bins[other_bin].links_to.size() - fill;
                                        bins[other_bin].links_to.resize(fill);
                                    } // TODO: if this stays, and links_to should be also checked, copy and paste previous 'if paragraph' from above
                                }
                            }
                            total_links += links.size();
                            gap_links_removed += links.size() - fill_pos;

                            std::cout << "AFTER DROP GAP LINKS (fill_pos: " << fill_pos << ")" << std::endl;
                            if (which_links == 0) {
                                bins[bin_id].links_to.resize(fill_pos);
                                for (uint64_t i=0; i!=bins[bin_id].links_to.size(); ++i) {
                                    print_link(bins[bin_id].links_to[i]);
                                }
                            } else {
                                bins[bin_id].links_from.resize(fill_pos);
                                for (uint64_t i=0; i!=bins[bin_id].links_from.size(); ++i) {
                                    print_link(bins[bin_id].links_from[i]);
                                }
                            }
                        }
                    //}


                    // for (uint64_t i = 0; i < links.size(); ++i) {
                    //     auto link = links[i];

                    //     if (link.first == 0 || link.second == 0)
                    //         continue;

                    //     if (link.first > link.second) {
                    //         links[fill_pos++] = link;
                    //         link_columns[link.first]++;
                    //         link_columns[link.second]++;
                    //         continue;
                    //     }

                    //     auto left_it = std::lower_bound(bin_ids.begin(), bin_ids.end(), link.first + 1);
                    //     auto right_it = std::lower_bound(bin_ids.begin(), bin_ids.end(), link.second);
                    //     if (right_it > left_it) {
                    //         links[fill_pos++] = link;
                    //         link_columns[link.first]++;
                    //         link_columns[link.second]++;
                    //     }
                    // }

                    // gap_links_removed += links.size() - fill_pos;
                    // links.resize(fill_pos);

                    handle_path(graph.get_path_name(path), bins, first_path);
                    first_path = false;
                });

                // iterate sorted bins and cumulatively sum up values in link_columns:
                std::unordered_map<uint64_t, uint64_t> cumsum_links;
                cumsum_links[0] = 0;
                for (uint64_t bin = 1; bin <= num_bins; ++bin) { // bin IDs start with 1
                    cumsum_links[bin] = cumsum_links[bin-1];
                    if (link_columns.find(bin-1) != link_columns.end()) {
        std::cout << "link_columns[" << bin << "] = " << link_columns[bin] << std::endl;
                        cumsum_links[bin] += link_columns[bin-1]/2; // divided by 2 because link is stored in both connecting bins, TODO check this is not true anymore?
                    }
                    handle_xoffset(bin, cumsum_links[bin], num_bins);
                }

                // collect bin sequences // TODO can be deleted
                for (uint64_t i = 0; i < num_bins; ++i) {
                    handle_sequence(i + 1, graph_seq.substr(i * bin_width, bin_width));
                    // 
                    std::cout << "bin " << i+1 << ": " << bin_pfreq[i].size() << std::endl;
                }
                // write out pangenome sequence
                handle_fasta(graph_seq);
                graph_seq.clear(); // clean up

                // entries in link_columns indicate 'component breakpoints'! can be used to split info at these points into chunk files

                uint64_t path_count = graph.get_path_count();

                std::cerr << std::setprecision(4) << "Gap links removed: " << (100.0 *  ((double)gap_links_removed / (double)total_links))
                            << "%, that is " << gap_links_removed << " gap links (" << path_count << " path start links + "
                            << path_count << " path end links + " << (gap_links_removed - path_count * 2) << " inner gap links) of "
                            << total_links << " total links" << std::endl;

        }
    }
}
