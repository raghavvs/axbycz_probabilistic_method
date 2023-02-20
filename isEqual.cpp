#include <cmath>

bool isequalf(double a, double b, double thresh = std::numeric_limits<double>::epsilon() * 100) {
    double diff = std::abs(a - b);
    return diff <= thresh;
}

/*
Note that the std::numeric_limits<double>::epsilon() method provides 
the smallest possible difference that can be detected between 
two floating-point numbers, which can be used as the default 
value for the thresh parameter. The implementation simply checks 
if the absolute difference between a and b is less than or equal 
to the threshold.
*/