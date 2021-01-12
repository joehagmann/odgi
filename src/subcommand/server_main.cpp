#include "subcommand.hpp"
#include "args.hxx"
#include "algorithms/xp.hpp"
#include <httplib.h>
//#include <filesystem.h>
//#include <experimental>
//namespace fs = std::experimental::filesystem;
#include <sys/types.h>
#include <dirent.h>

namespace odgi {

    using namespace odgi::subcommand;
    using namespace xp;
    using namespace httplib;

    int main_server(int argc, char** argv) {

        for (uint64_t i = 1; i < argc-1; ++i) {
            argv[i] = argv[i+1];
        }
        const std::string prog_name = "odgi server";
        argv[0] = (char*)prog_name.c_str();
        --argc;

        args::ArgumentParser parser("start a HTTP server with a given index file to query a pangenome position");
        args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
        args::ValueFlag<std::string> folder(parser, "FOLDER", "load the index from this folder. Must contain .px files.", {'f', "folder"});
        args::ValueFlag<std::string> port(parser, "N", "run the server under this port", {'p', "port"});
        args::ValueFlag<std::string> ip_address(parser, "IP", "run the server under this IP address", {'a', "ip"});

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

        if (!folder) {
            std::cerr << "[odgi server]: Please enter an existing folder containing path index files via -f=[FOLDER], --folder=[FOLDER]." << std::endl;
            exit(1);
        }

        if (!port) {
            std::cerr << "[odgi server]: Please enter a port for the server via -p=[N], --port=[N]." << std::endl;
            exit(1);
        }

        //std::unordered_map<std::string, XP*> data_sets;
        std::unordered_map<std::string, XP> data_sets;

        DIR* dirp = opendir(args::get(folder).c_str());
        struct dirent * dp;
        while ((dp = readdir(dirp)) != NULL) {
            std::string file = dp->d_name;
            if (file.find(".px") == std::string::npos) continue;
            
            std::string full_filename = args::get(folder) + "/" + file;
            std::cout << "File: " << full_filename << std::endl;
            
            XP path_index;
            std::ifstream in;
            in.open(full_filename);
                //std::cerr << "[odgi server]: Cannot load path index from " << file.path() << std::endl;
                //exit(1);
            //}
            // path_index.load(in);
            // in.close();

            size_t pos = file.find_last_of("."); 
            std::string dataset_name = file.substr(0, pos);
            std::cout << "Data set " << dataset_name << " stored." << std::endl;
            //data_sets[dataset_name] = &path_index;
            //data_sets[dataset_name] = new XP();
            data_sets[dataset_name].load(in);
            in.close();

            std::cout << "Loaded. " << std::endl;
            std::cout << "Pangenome pos of path 1, pos 16: " << data_sets[dataset_name].get_pangenome_pos("1", 16) << std::endl;

            //data_sets.insert(std::make_pair<std::string, XP>(filename.c_str(), path_index));
        }
        closedir(dirp);
        

        /*
        const char* pattern = R"(/(\d+)/(\w+))";
        std::regex regexi = std::regex(pattern);
        std::cmatch cm;    // same as std::match_results<const char*> cm;
        std::regex_match ("/3/test",cm,regexi);
        std::cout << "the matches were: " << std::endl;
        for (unsigned i=0; i<cm.size(); ++i) {
            std::cout << "[" << cm[i] << "] " << std::endl;
        }
        std::cout << std::endl;

        // const char* pattern1 = R"(/(\w+)/(\d+))";
        // const char* pattern1 = R"(/([a-zA-Z]*[0-9]*)/(\d+))";
        const char* pattern1 = R"(/(\w*.*)/(\d+))";
        std::regex regexi1 = std::regex(pattern1);
        std::cmatch cm1;    // same as std::match_results<const char*> cm;
        std::regex_match ("/5-/3",cm1,regexi1); // /1741.hr2/3
        std::cout << "the matches were: " << std::endl;
        for (unsigned i=0; i<cm1.size(); ++i) {
            std::cout << "[" << cm1[i] << "] " << std::endl;
        }
        */

        Server svr;

        svr.Get("/hi", [](const Request& req, Response& res) {
            res.set_header("Access-Control-Allow-Origin", "*");
            res.set_header("Access-Control-Expose-Headers", "text/plain");
            res.set_header("Access-Control-Allow-Methods", "GET, POST, DELETE, PUT");
            res.set_content("Hello World!", "text/plain");
            std::cout << "GOT REQUEST : HELLO WORLD!" << std::endl;
        });

        svr.Get(R"(/(\w*)/(\w*.*)/(\d+))", [&](const Request& req, Response& res) {
            
            for (size_t i = 0; i < req.matches.size(); i++) {
                std::cout << "req.matches " << req.matches[i] << std::endl;
            }
            
            auto data_set  = req.matches[1];
            auto path_name = req.matches[2];
            auto nuc_pos_1 = req.matches[3];
            std::cout << "GOT REQUEST : data set: " << data_set << " path name: " << path_name << "; 1-based nucleotide position: " << nuc_pos_1 << std::endl;
            size_t nuc_pos_0 = std::stoi(nuc_pos_1) - 1;
            size_t pan_pos = 0;

            // path_index is stored in data_sets[data_set]

            if (data_sets[data_set].has_path(path_name)) {
                if (data_sets[data_set].has_position(path_name, nuc_pos_0)) {
                    pan_pos = data_sets[data_set].get_pangenome_pos(path_name, nuc_pos_0) + 1;
                }
            }
            std::cout << "SEND RESPONSE: pangenome position: " << pan_pos << std::endl;
            res.set_header("Access-Control-Allow-Origin", "*");
            res.set_header("Access-Control-Expose-Headers", "text/plain");
            res.set_header("Access-Control-Allow-Methods", "GET, POST, DELETE, PUT");
            res.set_content(std::to_string(pan_pos), "text/plain");
        });

        svr.Get("/stop", [&](const Request& req, Response& res) {
            svr.stop();
        });

        const int p = std::stoi(args::get(port));
        std::string ip;
        if (!ip_address) {
            ip = "localhost";
        } else {
            ip = args::get(ip_address);
        }

        std::cout << "http server listening on http://" << ip << ":" << args::get(port) << std::endl;
        std::cout << "Pangenome pos of path 1, pos 16: " << data_sets["graph_MMk11_SQk10r10.Ygssorted"].get_pangenome_pos("1", 16) << std::endl;
        svr.listen(ip.c_str(), p);
        std::cout << "Pangenome pos of path 1, pos 50: " << data_sets["graph_MMk11_SQk10r10.Ygssorted"].get_pangenome_pos("1", 50) << std::endl;

        /*
        // we have a 0-based positioning
        uint64_t nucleotide_pos = args::get(nuc_pos) - 1;
        std::string p_name = args::get(path_name);

        if (!path_index.has_path(p_name)) {
            std::cerr << "The given path name " << p_name << " is not in the index." << std::endl;
            exit(1);
        }

        if (!path_index.has_position(p_name, nucleotide_pos)) {
            std::cerr << "The given path " << p_name << " with nucleotide position " << nuc_pos << " is not in the index." << std::endl;
            exit(1);
        }

        size_t pangenome_pos = path_index.get_pangenome_pos(p_name, nucleotide_pos) + 1;
        cout << pangenome_pos << endl;

        */
        return 0;
    }

    static Subcommand odgi_server("server",
                                  "start a HTTP server with a given folder to path index files to query a pangenome position",
                                  PIPELINE, 3, main_server);

}
