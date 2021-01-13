#pragma once

#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <iomanip> // std::setprecision
#include <handlegraph/handle_graph.hpp>
#include <handlegraph/util.hpp>
#include <handlegraph/path_handle_graph.hpp>

namespace odgi {
    namespace algorithms {

        using namespace std;
        using namespace handlegraph;

        struct link_info_t {
            int64_t from;
            int64_t to;
            int times;
            int nr_passes;
        };

        //int get_count_of_link_in_bin(const std::vector<link_info_t> &links, const link_info_t &link);

        struct bin_info_t {
            double mean_cov;
            double mean_inv;
            double mean_pos;
            std::vector<std::pair<uint64_t,uint64_t>> ranges;
            std::vector<link_info_t> links_to;
            std::vector<link_info_t> links_from;
            //std::unordered_set<std::map<std::pair<uint64_t, uint64_t>, int>> links_nr_passes;
            //std::unordered_map<link_t, int> links_nr_passes;
            //link_info_t last_link;
            // long int first_nucleotide;
            // long int last_nucleotide;
        };

        struct path_info_t {
            double mean_cov;
            double mean_inv;
            double mean_pos;
            std::vector<std::pair<uint64_t,uint64_t>> ranges;
        };

        void bin_path_info(const PathHandleGraph &graph,
                           const std::string &prefix_delimiter,
                           const std::function<void(const uint64_t &, const uint64_t &)> &handle_header,
                           const std::function<void(const std::string &,
                                                    const std::vector<std::pair<uint64_t, uint64_t>> &,
                                                    const std::map<uint64_t, algorithms::path_info_t> &)> &handle_path,
                           const std::function<void(const uint64_t &, const std::string &)> &handle_sequence,
                           const std::function<void(const std::string&)> &handle_fasta,
                           uint64_t num_bins = 0,
                           uint64_t bin_width = 0,
                           bool drop_gap_links = false);

        void bin_path_info_for_pantograph(const PathHandleGraph &graph,
                           const std::string &prefix_delimiter,
                           const std::function<void(const uint64_t &, const uint64_t &, const uint64_t &)> &handle_header,
                           const std::function<void(const std::string &,
                                                    const std::map<uint64_t, algorithms::bin_info_t> &,
                                                    const bool &, const uint64_t &)> &handle_path,
                           const std::function<void(const uint64_t &, const std::string &)> &handle_sequence,
                           const std::function<void(const std::string&, const uint64_t&)> &handle_fasta,
                           const std::function<void(const uint64_t &, const uint64_t &, const uint64_t &)> &handle_xoffset,
                           uint64_t num_bins = 0,
                           uint64_t bin_width = 0);
    }
}
