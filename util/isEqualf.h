/*
DESCRIPTION:

The program defines a function named isequalf that takes two double precision 
floating point numbers a and b as input arguments and returns a boolean value 
indicating whether they are equal within a specified threshold thresh. The 
function first checks if either of the input numbers is not a number (NaN) 
and returns false if so. Otherwise, it calculates the absolute difference 
between the two numbers, and checks if it is less than or equal to the specified 
threshold thresh. If the difference is within the threshold, the function 
returns true, otherwise it returns false. 
*/

#ifndef ISEQUAL_H
#define ISEQUAL_H

#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>

bool isEqualf(double a,
              double b,
              double thresh = std::numeric_limits<double>::epsilon()*100)
{
    if (std::isnan(a) || std::isnan(b)) return false;

    double m = std::max(std::abs(a-b), std::abs(b-a));

    if (m > thresh)
        return false;
    else
        return true;
}

#endif