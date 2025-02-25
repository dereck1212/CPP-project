#include "FDOptionPricer.h"
#include "DateUtils.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

// -----------------------------------------------------------------------------
// 1) Convert real calendar dates to year fractions for:
// • calculation start date => t=0
// • maturity date => T = yearFraction_Act365(calcDate, maturityDate)
// • piecewise rate dates => yearFractions[i] = yearFraction_Act365(calcDate, dateStrings[i])
// -----------------------------------------------------------------------------
static double computeTimeToMaturity(BSParams& p)
{
double T = yearFraction_Act365(p.calcDate, p.maturityDate);
if (T <= 0.0)
throw std::runtime_error("Maturity date is not after the calculation date.");
return T;
}

// -----------------------------------------------------------------------------
// Validate and convert the piecewise date/rates into yearFraction from p.calcDate
// times must be strictly ascending
// -----------------------------------------------------------------------------
static void validateAndComputeYearFractions(BSParams& p)
{
auto n = p.rCurve.dateStrings.size();
if (n == 0)
throw std::runtime_error("Rate curve: no points provided.");
if (n != p.rCurve.rates.size())
throw std::runtime_error("Rate curve: mismatch between dateStrings and rates.");

p.rCurve.yearFractions.resize(n);

for (size_t i=0; i<n; ++i)
{
    double tau = yearFraction_Act365(p.calcDate, p.rCurve.dateStrings[i]);
    p.rCurve.yearFractions[i] = tau;
}
// Check ascending
for (size_t i=1; i<n; ++i)
{
    if (p.rCurve.yearFractions[i] < p.rCurve.yearFractions[i-1])
        throw std::runtime_error("Rate curve times are not strictly ascending by date.");
}
}

// -----------------------------------------------------------------------------
// Interpolate the piecewise rate curve at time t.
// Extrapolate flat if t < 0 or t > last
// -----------------------------------------------------------------------------
static double interpolateRate(const PiecewiseRate& pc, double t)
{
// If only 1 point, return that rate
if (pc.yearFractions.size() == 1)
return pc.rates[0];

if (t <= pc.yearFractions.front())
    return pc.rates.front();
if (t >= pc.yearFractions.back())
    return pc.rates.back();

// Otherwise, find interval
for (size_t i=0; i<pc.yearFractions.size()-1; ++i)
{
    double t0=pc.yearFractions[i];
    double t1=pc.yearFractions[i+1];
    if (t>=t0 && t<=t1)
    {
        double r0=pc.rates[i];
        double r1=pc.rates[i+1];
        double w = (t - t0)/(t1 - t0);
        return r0 + w*(r1 - r0);
    }
}
// fallback
return pc.rates.back();
}

// -----------------------------------------------------------------------------
// Thomas solver for tridiagonal systems
// -----------------------------------------------------------------------------
static void thomasSolve(std::vector<double>& a,
std::vector<double>& b,
std::vector<double>& c,
std::vector<double>& d)
{
int n = (int)b.size();
// Forward elimination
for (int i = 1; i < n; ++i)
{
double w = a[i]/b[i-1];
b[i] -= w * c[i-1];
d[i] -= w * d[i-1];
}
// Back-substitution
d[n-1] = d[n-1]/b[n-1];
for (int i=n-2; i>=0; --i)
{
d[i] = (d[i] - c[i]*d[i+1]) / b[i];
}
}

// -----------------------------------------------------------------------------
// Compute Delta, Gamma, Theta at t=0 in a discrete sense
// P0[i] = P(S[i], t=0),
// P1[i] = P(S[i], t=dt)
// -----------------------------------------------------------------------------
static void computeGreeksAtT0(const std::vector<double>& P0,
const std::vector<double>& P1,
const std::vector<double>& S,
double dt,
PDEGreeks& out)
{
int n = (int)S.size();
out.Delta.resize(n, 0.0);
out.Gamma.resize(n, 0.0);
out.Theta.resize(n, 0.0);

// central differences for i in [1..n-2]
for (int i=1; i<n-1; ++i)
{
    double dS1 = (S[i+1] - S[i-1]);
    double delta = (P0[i+1] - P0[i-1]) / dS1;

    // gamma
    double dSleft  = (S[i] - S[i-1]);
    double dSright = (S[i+1] - S[i]);
    double denom   = dSleft*dSright; // if uniform same as (dS^2)
    double gamma   = 0.0;
    if (denom>1.0e-14)
    {
        gamma = (P0[i+1] - 2.0*P0[i] + P0[i-1]) / denom;
    }

    out.Delta[i] = delta;
    out.Gamma[i] = gamma;
}
// boundary i=0 => one-sided
out.Delta[0] = (P0[1] - P0[0])/(S[1]-S[0]);
out.Gamma[0] = 0.0;  
// boundary i=n-1 => one-sided
out.Delta[n-1] = (P0[n-1] - P0[n-2])/(S[n-1]-S[n-2]);
out.Gamma[n-1] = 0.0;

// Theta via forward diff in time: (P1 - P0)/dt
if (dt>1.0e-14)
{
    for(int i=0; i<n; ++i)
    {
        out.Theta[i] = (P1[i] - P0[i])/dt;
    }
}
}

// -----------------------------------------------------------------------------
// findEarlyExerciseBoundary: largest S[i] for which P[i]==payoff[i] within tol
// -----------------------------------------------------------------------------
static double findEarlyExerciseBoundary(const std::vector<double>& P,
const std::vector<double>& payoff,
const std::vector<double>& S)
{
double boundary = -1.0;
int n = (int)S.size();
for(int i=0; i<n; ++i)
{
if(std::fabs(P[i] - payoff[i])<1.0e-12)
{
boundary = std::max(boundary, S[i]);
}
}
return boundary;
}

// -----------------------------------------------------------------------------
// runCrankNicolsonPDE: single PDE run, returning 2D array P[j][i]
// -----------------------------------------------------------------------------
static std::vector<double> runCrankNicolsonPDE(BSParams& p, double T)
{
if (p.nTimeSteps<1 || p.nSpaceSteps<2)
throw std::runtime_error("Invalid grid steps.");
if (p.Smax<=0.0)
throw std::runtime_error("Smax must be positive.");
if (T<=0.0)
throw std::runtime_error("Time to maturity must be > 0.");
if (p.sigma<=0.0)
throw std::runtime_error("Volatility must be > 0.");

double dt = T/p.nTimeSteps;
double dS = p.Smax/p.nSpaceSteps;

// Build the space grid
std::vector<double> S(p.nSpaceSteps+1);
for(int i=0; i<=p.nSpaceSteps; ++i)
    S[i] = i*dS;

// 2D array for solution
// Pdata[j*(nSpaceSteps+1) + i] => P[j][i]
std::vector<double> Pdata((p.nTimeSteps+1)*(p.nSpaceSteps+1), 0.0);
auto idx = [&](int j,int i){return j*(p.nSpaceSteps+1)+i;};

// payoff
std::vector<double> payoff(p.nSpaceSteps+1);
for(int i=0; i<=p.nSpaceSteps; ++i)
{
    double callVal = std::max(S[i] - p.K, 0.0);
    double putVal  = std::max(p.K - S[i], 0.0);
    payoff[i] = (p.optionType=="call") ? callVal : putVal;
}

// terminal condition => j=nTimeSteps
for(int i=0; i<=p.nSpaceSteps; ++i)
    Pdata[idx(p.nTimeSteps,i)] = payoff[i];

// backward in time
for(int j=p.nTimeSteps-1; j>=0; --j)
{
    // time midpoint for the interest rate => tMid
    double tMid = ( (j+0.5)*dt );
    double rMid = interpolateRate(p.rCurve, tMid);

    double sigma2 = p.sigma*p.sigma;

    // build tri-diag
    std::vector<double> a(p.nSpaceSteps+1, 0.0);
    std::vector<double> b(p.nSpaceSteps+1, 0.0);
    std::vector<double> c(p.nSpaceSteps+1, 0.0);
    std::vector<double> d(p.nSpaceSteps+1, 0.0);

    auto alpha = [&](int i){
        return 0.5*dt*(sigma2*S[i]*S[i]/(dS*dS) - rMid*S[i]/dS);
    };
    auto beta = [&](int i){
        return - dt*(sigma2*S[i]*S[i]/(dS*dS) + rMid);
    };
    auto gamma = [&](int i){
        return 0.5*dt*(sigma2*S[i]*S[i]/(dS*dS) + rMid*S[i]/dS);
    };

    // boundary conditions
    // use approximate
    double discFactor = std::exp(-rMid*(T - j*dt));
    double bcLeft=0.0, bcRight=0.0;
    if(p.optionType=="call")
    {
        // call => P(0)=0, P(Smax)= Smax - K e^{-r(T-t)}
        bcLeft = 0.0;
        bcRight= S[p.nSpaceSteps] - p.K*discFactor;
        if(bcRight<0.0) bcRight=0.0;
    }
    else
    {
        // put => P(0)= K e^{-r(T-t)}, P(Smax)=0
        bcLeft = p.K*discFactor;
        bcRight=0.0;
    }

    // fill a,b,c
    for(int i=0; i<=p.nSpaceSteps; ++i)
    {
        b[i] = 1.0 - 0.5*beta(i);
        a[i] = -0.5*alpha(i);
        c[i] = -0.5*gamma(i);
    }

    // build RHS from old => j+1
    for(int i=1; i<p.nSpaceSteps; ++i)
    {
        double A = 0.5*alpha(i);
        double B = 1.0 + 0.5*beta(i);
        double C = 0.5*gamma(i);
        double P_im1 = Pdata[idx(j+1,i-1)];
        double P_i   = Pdata[idx(j+1,i)];
        double P_ip1 = Pdata[idx(j+1,i+1)];
        d[i] = B*P_i + A*P_im1 + C*P_ip1;
    }

    // boundary i=0 => value known
    d[0] = bcLeft;
    b[0] = 1.0; c[0]=0.0; a[0]=0.0;

    // boundary i=nSpaceSteps => value known
    d[p.nSpaceSteps] = bcRight;
    b[p.nSpaceSteps] = 1.0; a[p.nSpaceSteps]=0.0; c[p.nSpaceSteps]=0.0;

    // solve tri-diag
    thomasSolve(a,b,c,d);

    // write solution => j
    for(int i=0; i<=p.nSpaceSteps; ++i)
        Pdata[idx(j,i)] = d[i];

    // if American => max(P, payoff)
    if(p.exerciseType=="american")
    {
        for(int i=0; i<=p.nSpaceSteps; ++i)
        {
            if(Pdata[idx(j,i)]<payoff[i])
                Pdata[idx(j,i)] = payoff[i];
        }
    }
} // end time loop

return Pdata;
}

// -----------------------------------------------------------------------------
// Bump-and-reprice for Rho or Vega at S0
// -----------------------------------------------------------------------------
static double bumpAndRepriceRho(BSParams& p, double T, double bumpSize)
{
// Copy p
BSParams bumped = p;
// Bump the entire curve by bumpSize
for(size_t i=0; i<bumped.rCurve.rates.size(); ++i)
bumped.rCurve.rates[i] += bumpSize;

std::vector<double> PDE      = runCrankNicolsonPDE(p,       T);
std::vector<double> PDEbump  = runCrankNicolsonPDE(bumped,  T);

// same grid => just get price at S0
auto idx = [&](int j,int i){return j*(p.nSpaceSteps+1)+i;};
double dS = p.Smax/p.nSpaceSteps;
auto spotIndexPrice = [&](const std::vector<double>& arr)->double {
    if(p.S0<=0.0) return arr[idx(0,0)];
    if(p.S0>=p.Smax) return arr[idx(0,p.nSpaceSteps)];
    int iLow = (int)(p.S0/dS);
    int iHigh= iLow+1;
    double w=(p.S0-(iLow*dS))/dS;
    return (1.0-w)*arr[idx(0,iLow)] + w*arr[idx(0,iHigh)];
};

double base = spotIndexPrice(PDE);
double bump= spotIndexPrice(PDEbump);

return (bump-base)/bumpSize;
}

static double bumpAndRepriceVega(BSParams& p, double T, double bumpSize)
{
BSParams bumped = p;
bumped.sigma+=bumpSize;
if(bumped.sigma<=0.0) bumped.sigma=1.0e-8;

std::vector<double> PDE     = runCrankNicolsonPDE(p,      T);
std::vector<double> PDEbump = runCrankNicolsonPDE(bumped, T);

auto idx = [&](int j,int i){return j*(p.nSpaceSteps+1)+i;};
double dS = p.Smax/p.nSpaceSteps;
auto spotIndexPrice = [&](const std::vector<double>& arr)->double {
    if(p.S0<=0.0) return arr[idx(0,0)];
    if(p.S0>=p.Smax) return arr[idx(0,p.nSpaceSteps)];
    int iLow = (int)(p.S0/dS);
    int iHigh= iLow+1;
    double w=(p.S0-(iLow*dS))/dS;
    return (1.0-w)*arr[idx(0,iLow)] + w*arr[idx(0,iHigh)];
};

double base = spotIndexPrice(PDE);
double bump= spotIndexPrice(PDEbump);

return (bump-base)/bumpSize;
}

// -----------------------------------------------------------------------------
// Exposed function: PriceVanillaOptionCN
// Reads real dates, converts them to T, builds PDE, etc.
// -----------------------------------------------------------------------------
FDSolution* PriceVanillaOptionCN(const BSParams& origParams)
{
// We copy the input param struct because we need to mutate some fields
BSParams p = origParams;

// 1) Convert date strings => T, yearFractions
double T= computeTimeToMaturity(p);
validateAndComputeYearFractions(p);

// 2) Solve PDE
std::vector<double> PDE = runCrankNicolsonPDE(p, T);

// 3) Spot grid
double dS = p.Smax/p.nSpaceSteps;
std::vector<double> S(p.nSpaceSteps+1);
for(int i=0; i<=p.nSpaceSteps; ++i)
    S[i] = i*dS;

// Helper to index PDE
auto idx = [&](int j,int i){return j*(p.nSpaceSteps+1)+i;};

// Price at S0 => j=0
double priceAtS0=0.0;
// Linear interpolation
if(p.S0<=0.0)
    priceAtS0 = PDE[idx(0,0)];
else if(p.S0>=p.Smax)
    priceAtS0 = PDE[idx(0,p.nSpaceSteps)];
else
{
    int iLow  = (int)(p.S0/dS);
    int iHigh = iLow+1;
    double w = (p.S0-(iLow*dS))/dS;
    priceAtS0= (1.0-w)*PDE[idx(0,iLow)] + w*PDE[idx(0,iHigh)];
}

// Extract price array at t=0 => j=0 for all i
std::vector<double> price0(p.nSpaceSteps+1);
for(int i=0; i<=p.nSpaceSteps; ++i)
    price0[i] = PDE[idx(0,i)];

// For Theta at t=0, we need j=1 => if nTimeSteps=1, j=1 => T => still possible
std::vector<double> price1(p.nSpaceSteps+1);
if(p.nTimeSteps>=1)
{
    for(int i=0; i<=p.nSpaceSteps; ++i)
        price1[i]= PDE[idx(1,i)];
}
else
{
    throw std::runtime_error("Not enough time steps for Theta.");
}

// compute ∆, Γ, Θ at t=0 (for each S)
PDEGreeks gre;
double dt= T/p.nTimeSteps;
computeGreeksAtT0(price0,price1,S,dt,gre);

// 4) Rho, Vega at S0
double dr=1.0e-4;
double dSigma=1.0e-4;
double rho  = bumpAndRepriceRho(p, T, dr);
double vega = bumpAndRepriceVega(p, T, dSigma);

// 5) If American => early boundary
std::vector<double> eebound(p.nTimeSteps+1, -1.0);
if(p.exerciseType=="american")
{
    // payoff
    std::vector<double> payoff(p.nSpaceSteps+1);
    for(int i=0; i<=p.nSpaceSteps; ++i)
    {
        double callVal = std::max(S[i]-p.K,0.0);
        double putVal  = std::max(p.K-S[i],0.0);
        payoff[i]=(p.optionType=="call")?callVal:putVal;
    }
    for(int j=0; j<=p.nTimeSteps; ++j)
    {
        std::vector<double> slice(p.nSpaceSteps+1);
        for(int i=0; i<=p.nSpaceSteps; ++i)
            slice[i]= PDE[idx(j,i)];
        double b = findEarlyExerciseBoundary(slice,payoff,S);
        eebound[j]= b;
    }
}

// build final solution
FDSolution* sol = new FDSolution();
sol->priceAtS0 = priceAtS0;
sol->rhoAtS0   = rho;
sol->vegaAtS0  = vega;
sol->spotGrid  = S;
sol->priceGridAtT0 = price0;
sol->greeksAtT0    = gre;
sol->earlyExerciseBoundary = eebound;

return sol;
}
