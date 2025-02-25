#include "FDOptionPricer.h"
#include "DateUtils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>

// -----------------------------------------------------------------------------
// Example function to read our BSParams from a text file that includes real date
// fields and rate curve as lines, e.g.:
//
// calcDate=2023-10-01
// maturityDate=2024-10-01
// S0=100
// K=100
// sigma=0.2
// nTimeSteps=200
// nSpaceSteps=200
// Smax=400
// optionType=call
// exerciseType=american
// rCurve.dates=2023-10-01,2024-01-01,2024-10-01
// rCurve.rates=0.03,0.035,0.04
//
// You can adapt as needed for your style of input file.
// -----------------------------------------------------------------------------
extern "C" void readParameters(const char* filename, BSParams* p)
{
if (!p) {
    std::cerr << "Error in readParameters()" << std::endl;
    return;
}

std::ifstream ifs(filename);
if(!ifs.is_open()) {
    throw std::runtime_error("Cannot open input file: "+std::string(filename));
}
while(!ifs.eof())
{
    std::string line;
    if(!std::getline(ifs,line)) break;
    if(line.empty()) continue;

    // parse "key=value"
    std::string key, val;
    {
        std::istringstream iss(line);
        if(!std::getline(iss,key,'=')) continue;
        if(!std::getline(iss,val)) continue;
    }

    if(key=="Date")           p->calcDate        = val;
    else if(key=="Maturity Date")  p->maturityDate    = val;
    else if(key=="S0")           p->S0              = std::stod(val);
    else if(key=="K")            p->K               = std::stod(val);
    else if(key=="Sigma")        p->sigma           = std::stod(val);
    else if(key=="Time Steps")    p->nTimeSteps      = std::stoi(val);
    else if(key=="Space Steps")   p->nSpaceSteps     = std::stoi(val);
    else if(key=="Smax")         p->Smax            = std::stod(val);
    else if(key=="Option Type")    p->optionType      = val;
    else if(key=="Exercise Type")  p->exerciseType    = val;
    else if(key=="Curve Dates")
    {
        // parse comma-separated dates
        std::istringstream issDates(val);
        while(!issDates.eof())
        {
            std::string dStr;
            if(!std::getline(issDates,dStr,',')) break;
            if(!dStr.empty()) p->rCurve.dateStrings.push_back(dStr);
        }
    }
    else if(key=="Curve Rate")
    {
        // parse comma-separated rates
        std::istringstream issRates(val);
        while(!issRates.eof())
        {
            std::string rStr;
            if(!std::getline(issRates,rStr,',')) break;
            if(!rStr.empty()) p->rCurve.rates.push_back(std::stod(rStr));
        }
    }
    // else ignore unrecognized keys
}
}

// -----------------------------------------------------------------------------
// Example function to write results (including all Greeks at t=0) to an output file
// so that the user can plot them in Excel or another tool.
// -----------------------------------------------------------------------------
extern "C" void writeOutput(const std::string& filename, const FDSolution& sol)
{
std::ofstream ofs(filename);
if(!ofs.is_open())
throw std::runtime_error("Cannot open output file: "+filename);

ofs << "1) Prix theorique P(S0) = " << sol.priceAtS0 << "\n";
ofs << "2) Delta at S0 not separately shown here, see arrays below.\n";
ofs << "3) Gamma at S0 likewise.\n";
ofs << "4) Theta at S0 likewise.\n";
ofs << "5) Rho at S0 = " << sol.rhoAtS0 << "\n";
ofs << "6) Vega at S0 = " << sol.vegaAtS0 << "\n\n";

// Now dump the arrays for each S => we can do S, P, Delta, Gamma, Theta
ofs << "Spot, Price(t=0), Delta, Gamma, Theta\n";
int n = (int)sol.spotGrid.size();
for(int i=0; i<n; ++i)
{
    ofs << sol.spotGrid[i] << ","
        << sol.priceGridAtT0[i] << ","
        << sol.greeksAtT0.Delta[i] << ","
        << sol.greeksAtT0.Gamma[i] << ","
        << sol.greeksAtT0.Theta[i] << "\n";
}

// For an American option, also write the early-ex boundary
if(!sol.earlyExerciseBoundary.empty())
{
    ofs << "\nAmerican Early-Exercise Boundary per time step:\n";
    double finalTime = 1.0; // It's not necessarily 1.0 if T != 1. But we can’t see T directly. 
    // If we wanted to be thorough, we’d store T in FDSolution. For demonstration, we just show j indices.
    ofs << "StepIndex, BoundaryInSpot\n";
    for(size_t j=0; j<sol.earlyExerciseBoundary.size(); ++j)
    {
        ofs << j << "," << sol.earlyExerciseBoundary[j] << "\n";
    }
}
}
