#include "subcommand.hpp"
#include "args.hxx"
#include "algorithms/xp.hpp"
#include <httplib.h>
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
            data_sets[dataset_name].load(in);
            in.close();

            std::cout << "Loaded. " << std::endl;
            std::cout << "Pangenome pos of path 1, pos 50: " << data_sets[dataset_name].get_pangenome_pos("1", 50) << std::endl;
        }
        closedir(dirp);
        

        Server svr;

        svr.Get("/hi", [](const Request& req, Response& res) {
            res.set_header("Access-Control-Allow-Origin", "*");
            res.set_header("Access-Control-Expose-Headers", "text/plain");
            res.set_header("Access-Control-Allow-Methods", "GET, POST, DELETE, PUT");
            res.set_content("Hello World!", "text/plain");
            std::cout << "GOT REQUEST : HELLO WORLD!" << std::endl;
        });

        svr.Get(R"(/(.+)/(.+)/(\d+))", [&](const Request& req, Response& res) {
            
            for (size_t i = 0; i < req.matches.size(); i++) {
                std::cout << "req.matches " << req.matches[i] << std::endl;
            }
            
            auto data_set  = req.matches[1];
            auto path_name = req.matches[2];
            auto nuc_pos_1 = req.matches[3];
            std::cout << "GOT REQUEST : data set: " << data_set << " path name: " << path_name << "; 1-based nucleotide position: " << nuc_pos_1 << std::endl;

            std::string output_json = "{\"pan_pos\":";
            size_t nuc_pos_0 = std::stoi(nuc_pos_1) - 1;
            size_t pan_pos = 0;

            // path_index is stored in data_sets[data_set]

            if (data_sets.find(data_set) != data_sets.end()) {
                if (data_sets[data_set].has_path(path_name)) {
                    if (data_sets[data_set].has_position(path_name, nuc_pos_0)) {
                        pan_pos = data_sets[data_set].get_pangenome_pos(path_name, nuc_pos_0) + 1;
                        if (pan_pos != 0) output_json += std::to_string(pan_pos) + ",\"exitcode\":0,\"errmsg\":\"\"}\n";
                    }
                    else {
                        output_json += "0,\"exitcode\":1,\"errmsg\":\"Path " + static_cast<std::string>(path_name) + 
                                       " does not have position " + static_cast<std::string>(nuc_pos_1) + "\"}\n";
                    }
                }
                else {
                    output_json += "0,\"exitcode\":1,\"errmsg\":\"Dataset " + static_cast<std::string>(data_set) + 
                                   " does not have path " + static_cast<std::string>(path_name) + "\"}\n";
                }
            }
            else {
                output_json += "0,\"exitcode\":1,\"errmsg\":\"There is no path index for dataset " + 
                               static_cast<std::string>(data_set) + "\"}\n";
            }

            std::cout << "SEND RESPONSE: pangenome position: " << pan_pos << std::endl;
            res.set_header("Access-Control-Allow-Origin", "*");
            res.set_header("Access-Control-Expose-Headers", "text/plain");
            res.set_header("Access-Control-Allow-Methods", "GET, POST, DELETE, PUT");
            res.set_content(output_json, "text/plain");
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
        //std::cout << "Pangenome pos of path 1, pos 16: " << data_sets["graph_MMk11_SQk10r10.Ygssorted.2"].get_pangenome_pos("1", 16) << std::endl;
        svr.listen(ip.c_str(), p);

        return 0;
    }

    static Subcommand odgi_server("server",
                                  "start a HTTP server with a given folder to path index files to query a pangenome position",
                                  PIPELINE, 3, main_server);

}
