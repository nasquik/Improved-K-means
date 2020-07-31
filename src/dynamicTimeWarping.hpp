#ifndef DYNAMICTIMEWARPING_HPP
#define DYNAMICTIMEWARPING_HPP

#include <vector>
using std::vector;
using std::pair;

#include <utility>
using std::pair;

double dtwDist(class Curve *, class Curve *);
vector<pair<int, int> >* dtwBestTraversal(class Curve*, class Curve*);
vector <pair<int, int> > * optimalPath(double**, int, int);

#endif