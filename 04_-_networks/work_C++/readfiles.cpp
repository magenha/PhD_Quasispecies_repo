
//This code is intended to compute the distance matrix of a population's GN
//Compute the mean distance using the abundance of nodes

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <sys/stat.h>
#include <list>
namespace fs = std::filesystem;
using namespace std; 

list<string> files_in_directory(string path_dir){
    //Read files
    std::string path = path_dir;
    struct stat sb;
    list<string> files;
 
    for (const auto& entry : fs::directory_iterator(path)) {
        std::filesystem::path outfilename = entry.path();
        std::string outfilename_str = outfilename.string();
        const char* path = outfilename_str.c_str();
 
        // Testing whether the path points to a
        // non-directory or not If it does, displays path
        if (stat(path, &sb) == 0 && !(sb.st_mode & S_IFDIR))
            //std ::cout << path << std::endl;
            files.push_back(path);
    }

    return files;
}

list<string> read_file(string file){
    std::ifstream myfile;
    list<string> lines;
    myfile.open(file);
    std::string row;
    if (myfile.is_open()){
        while (myfile){
            std::getline (myfile, row);
            lines.push_back(row);
        }
    }
    return lines;
}

//
class Graph{
   
    public:
    list<int> nodes_list;
    int adj_matrix;
    

    void load_edges_list(list<string> Input_list){
        for (auto const &i: Input_list){
            std::string s = i;
            std::string delimiter = ", ";
            size_t pos = 0;
            std::string token;
            list<int> edges;
            while ((pos = s.find(delimiter)) != std::string::npos){
                token = s.substr(0,pos);
                edges.push_back(stoi(token));
                s.erase(0,pos+delimiter.length());
            }
            list<int>::iterator it=edges.begin();
            int e1 = *it;
            advance(it,1);
            int e2 = *it;
            adj_matrix[e1][e2] = {1};
        }
    }
};






int main()
{
    //Define path for data
    
    std::string path = "./SubNetworks_nNodes-10";

    
    // 1- Read the edges of the network
    //Read files
    list<string> files = files_in_directory(path);
    std::string file_sample = path + "/edges.csv";
    list<string> lines_in_file = read_file(file_sample);
    //Create graph with certain edges
    Graph G;
    G.load_edges_list(lines_in_file);
    //G.print_edges();

    return 0;
} 




// 2- Select the biggest connected component
// 3- Compute distance matrix
// 4- Compute distance dispersion
// 4.1 - Read abundances for each sequence
// 4.2 - Do computations
// 5- Save results