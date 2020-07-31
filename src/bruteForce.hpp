#ifndef _BRUTE_FORCE_
#define _BRUTE_FORCE_

#include <vector>
using std::vector;

#include "dataStructs.hpp"

using namespace std;

pair<class Point *, double> *bruteForce(vector<class Point *> *, class Point *);
int bruteForceAll(vector<class Point *> *, vector<class Point *> *);
pair<class Curve *, double> *bruteForceCurve(vector<class Curve *> *, class Curve *);
int calculateW(vector<class Point *> *points, int percentage);

#endif