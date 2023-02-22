#include <iostream>
#include <string>
#include <fstream>

//const std::string input_file = "D:/Geomatics/Q3/GEO1004 3D Modelling/assignments/1/data/station.hw1";
const std::string input_file = "station.hw1";


int main() {
    std::cout << "input file name: " << input_file << std::endl;

    std::ifstream input_stream;
    input_stream.open(input_file);

    bool is_file_open {true};
    is_file_open = input_stream.is_open();
    std::cout << std::boolalpha;
    std::cout << "is file open: " << is_file_open << std::endl;

    std::string line1;
    std::getline(input_stream, line1);
    std::cout << "line1: " << line1 << std::endl;

    return 0;
}