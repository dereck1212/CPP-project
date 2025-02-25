#ifndef DATEUTILS_H
#define DATEUTILS_H

#include <string>
#include <chrono>
#include <stdexcept>

// -----------------------------------------------------------------------------
// Very simple date-parsing function. We assume format "YYYY-MM-DD" strictly.
// Throws if parsing fails, or if out of range. (No leap-year checks here.)
// -----------------------------------------------------------------------------
inline std::tm parseDate(const std::string& dateStr)
{
// Expect dateStr like "2024-01-15"
// We'll parse into a std::tm struct.
// NOTE: This is a simplistic approach that does not handle time zones, etc.
std::tm tm_val = {};
int year, month, day;
if (sscanf(dateStr.c_str(), "%d-%d-%d", &year, &month, &day) != 3)
{
throw std::runtime_error("parseDate: Invalid date format: " + dateStr);
}
tm_val.tm_year = year - 1900; // tm_year=years since 1900
tm_val.tm_mon = month - 1; // 0-based
tm_val.tm_mday = day; // 1-based
// The rest fields are left as 0
return tm_val;
}

// -----------------------------------------------------------------------------
// Convert a std::tm to a time_t, ignoring time zones. We set the hour, minute, second to 0
// -----------------------------------------------------------------------------
inline std::time_t toTimeT(const std::tm& tm_val)
{
std::tm copy = tm_val;
copy.tm_hour = 0; copy.tm_min = 0; copy.tm_sec=0; // midnight
return std::mktime(&copy); // system-dependent, but typically OK
}

// -----------------------------------------------------------------------------
// Compute a day count fraction between two dates using a simple "Actual/365" convention
// or something similar. We'll do Actual/365 for demonstration.
// -----------------------------------------------------------------------------
inline double yearFraction_Act365(const std::string& date1, const std::string& date2)
{
std::tm tm1 = parseDate(date1);
std::tm tm2 = parseDate(date2);
std::time_t t1 = toTimeT(tm1);
std::time_t t2 = toTimeT(tm2);

// Number of days difference => double
double diff = std::difftime(t2, t1)/ (60.0*60.0*24.0);
// Actual/365 fraction
return diff / 365.0;
}

#endif // DATEUTILS_H
