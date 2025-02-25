#include "FDOptionPricer.h"
#include "DateUtils.h"
#include <iostream>
#include <exception>
#include <string>

// We also declare the external I/O functions from DataIO.cpp
BSParams readParameters(const std::string& filename);
void writeOutput(const std::string& filename, const FDSolution& sol);

int main()
{
try
{
// Example usage: read from an "input.txt", run solver, write "output.txt"
const char* inFile = "/Users/jameskouchade/Desktop/M2QF/Cours C++/Projet/Dereck/input.txt";
std::string outFile = "/Users/jameskouchade/Desktop/M2QF/Cours C++/Projet/Dereck/output.txt";

    BSParams p ;
    readParameters(inFile, &p);
    FDSolution* sol = PriceVanillaOptionCN(p);

    writeOutput(outFile, *sol);

    std::cout << "Success! Wrote results to " << outFile << std::endl;
}
catch(std::exception &e)
{
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
}
return 0;
}
