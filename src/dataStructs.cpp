#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <random>

#include "dataStructs.hpp"

using namespace std;

/////// POINTS //////

Point::Point(string ID, double *coord, int siz)
{
    this->ID = ID;
    coordinates = coord;
    size = siz;
    isCentroid = 0;
}

Point::Point(string ID, double *coord, int siz,class Curve* mycurve)
{
    this->ID = ID;
    coordinates = coord;
    size = siz;
    isCentroid = 0;
    curvePtr = mycurve;
}


Point::Point(string ID, double *coord, int siz, int cent)
{
    this->ID = ID;
    coordinates = coord;
    size = siz;
    isCentroid = 1;
}

Point::~Point()
{
    delete[] coordinates;
}

string Point::getID()
{
    return ID;
}

int Point::getSize()
{
    return size;
}

double *Point::getCoord()
{
    return coordinates;
}

class Curve *Point::getCurvePtr(){
    return curvePtr;
}

int Point::Centroid()
{
    return isCentroid;
}
////// CURVES //////

Curve::Curve(string givenID, pair<double, double> *coords, int size, int isCentroid)
{
    ID = givenID;
    coordinates = coords;
    this->size = size;
    this->isCentroid = isCentroid;
}

pair<double, double> *Curve::getCoord()
{
    return coordinates;
}

int Curve::getSize()
{
    return size;
}

Curve::~Curve()
{
    delete [] coordinates;
}

string Curve::getID()
{
    return ID;
}

int Curve::Centroid()
{
    return isCentroid;
}
