#ifndef FILEREADING_HPP
#define FILEREADING_HPP

#include <vector>
#include <string>

using std::vector;
using std::string;

class Reading

{
public:
    vector<class Point*>* readPoints(string);
    vector<class Curve*>* readCurves(string);
};

int readConfig(string, int *, int *, int *, int *);

#endif