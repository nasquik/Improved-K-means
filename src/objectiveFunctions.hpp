#ifndef OBJECTIVE_FUNCTIONS_HPP
#define OBJECTIVE_FUNCTIONS_HPP

#include <unordered_set>
using std::unordered_set;

double minDistObjective(class Point*, unordered_set<class Point*>);
double minDistObjective(class Curve*, unordered_set<class Curve*>);

#endif