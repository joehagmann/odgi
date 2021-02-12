#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include "GeneAnnotation.h"

using std::string;

int main()
{
	string gff3_fn = "ath.gff3";

    GeneAnnotation anno;
    std::cout << "source:" << anno.anno_source << ":" << std::endl;

    try {
        anno.parse_gff3_annotation(gff3_fn, false);
    }
    catch(const std::exception&) {
        return 1;
    }

    // std::ifstream infile(file);
	// std::string line;
   	// while (getline(infile, line))
	// {
	// 	// skip header lines
    //     if (line[0] == '#')
    //     {
    //         continue;
    //     }

	// 	std::cout << line << std::endl;

	// 	const std::vector<std::string> fields = split(line, '\t');
	// 	const std::vector<std::string> attributes = split(fields[ATTRIBUTES], ';');


	// }

   return 0;
}
