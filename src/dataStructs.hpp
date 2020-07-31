#ifndef DATASTRUCTS_HPP
#define DATASTRUCTS_HPP

#include <string>
using std::string;
#include <utility>
using std::pair;
#include <vector>
using std::vector;

class Curve
{
    string ID;
    pair<double, double> *coordinates;
    int size;
    char isCentroid;

public:
    Curve(){};
    ~Curve();
    Curve(string, pair<double, double> *, int, int);
    pair<double, double> *getCoord();
    int getSize();
    string getID();
    int Centroid();
};

class Point
{
    string ID;
    double *coordinates;
    int size;
    class Curve* curvePtr;
    char isCentroid;

public:
    Point();
    Point(string ID, double *coord, int siz);
    Point(string ID, double *coord, int siz, class Curve*);

    Point(string ID, double *coord, int siz, int cent);
    ~Point();
    string getID();
    double *getCoord();
    int getSize();
    class Curve *getCurvePtr();
    int Centroid();
};

#endif