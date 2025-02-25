#ifndef FDOPTIONPRICER_H
#define FDOPTIONPRICER_H

#include <vector>
#include <string>
#include <stdexcept>
#include <chrono>

// -----------------------------------------------------------------------------
// A structure for storing date/rate pairs in real calendar format.
// For example, date[i] might be "2024-01-15" and rates[i] = 0.05.
// We'll convert them internally into "year fraction" from a given start date.
// -----------------------------------------------------------------------------
struct PiecewiseRate
{
// Each element is a (date, rate).
// We'll store them in parallel arrays for simplicity.
std::vector<std::string> dateStrings; // e.g. {"2023-10-01", "2024-01-15", "2025-01-01"}
std::vector<double> rates; // e.g. {0.03, 0.035, 0.04}

// After reading from file, we convert dateStrings to double times (year fractions)
// measured from the calculation start date. We'll do that behind the scenes.
std::vector<double>      yearFractions;  // same size as dateStrings & rates
};

// -----------------------------------------------------------------------------
// A structure to hold the input parameters for our PDE solver
// -----------------------------------------------------------------------------
struct BSParams
{
// Real dates as strings, e.g. "2024-01-15"
// We'll convert them into double year-fractions (T).
std::string calcDate; // The calculation start date (default = "today")
std::string maturityDate; // The option maturity date

// Spot and option parameters
double S0;     // Current spot price
double K;      // Strike
double sigma;  // Volatility (constant, though we can bump it)

// The total time to maturity (in years) once we parse the dates => T
// We compute T = (maturityDate - calcDate) in year-fraction.
// We'll store it in T internally, or compute on the fly.
// (No direct field "T" here, to avoid confusion. We'll compute it once we have real dates.)

// Piecewise linear risk-free rate curve, in real calendar format
PiecewiseRate rCurve;

// Grid parameters
int    nTimeSteps;     // Number of time steps
int    nSpaceSteps;    // Number of space (spot) steps
double Smax;           // Maximum spot boundary

// Contract details
std::string optionType;    // "call" or "put"
std::string exerciseType;  // "european" or "american"
};

// -----------------------------------------------------------------------------
// A structure that holds the Greeks at a given time layer. For t=0, we have
// Delta[i] = ∂P/∂S at S[i]
// Gamma[i] = ∂²P/∂S²
// Theta[i] = ∂P/∂t
// If you wish to store these for every time step, you’d keep a 2D array of PDEGreeks.
// -----------------------------------------------------------------------------
struct PDEGreeks
{
// Same length as the spatial grid
std::vector<double> Delta;
std::vector<double> Gamma;
std::vector<double> Theta;
};

// -----------------------------------------------------------------------------
// A structure that holds the final PDE solution plus requested outputs
// -----------------------------------------------------------------------------
struct FDSolution
{
// Price at S0, t=0
double priceAtS0;

// Rho at S0 => ∂P/∂r at S0, t=0
double rhoAtS0;

// Vega at S0 => ∂P/∂σ at S0, t=0
double vegaAtS0;

// Entire spot grid
std::vector<double> spotGrid;

// PDE solution at t=0 for each S => P( S[i], t=0 )
std::vector<double> priceGridAtT0;

// Greeks at t=0 for each S => for plotting ∆(S), Γ(S), Θ(S)
PDEGreeks greeksAtT0;

// For American options: early-exercise boundary as a function of time step j
// earlyExerciseBoundary[j] = that boundary in spot
// If none is exercised, store -1 or 0 as you see fit.
std::vector<double> earlyExerciseBoundary;
};

// -----------------------------------------------------------------------------
// The main PDE solver function to be exposed
// -----------------------------------------------------------------------------
#ifdef _WIN32
extern "C" __declspec(dllexport)
#else

#endif
extern "C" {
    void readParameters(const char* filename, BSParams* p);
    void writeOutput(const std::string& filename, const FDSolution& sol);
    FDSolution* PriceVanillaOptionCN(const BSParams& params);
}

#endif // FDOPTIONPRICER_H
