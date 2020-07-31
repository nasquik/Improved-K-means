#include <cmath>
#include <iostream>

#include "dataStructs.hpp"
#include "metrics.hpp"

#include <unordered_set>

using namespace std;

double manhattanDist(class Point *Point1, class Point *Point2)
{
    double distance = 0;
    unsigned int dimensions = Point1->getSize();
    double *coord_1 = Point1->getCoord();
    double *coord_2 = Point2->getCoord();

    for (int i = 0; i < dimensions; i++)
    {
        distance += abs(coord_1[i] - coord_2[i]);
    }

    return distance;
}

double manhattanDist(pair<double, double> Point1, pair<double, double> Point2)
{
    double distance = 0;
    distance += abs(Point1.first - Point2.first);
    distance += abs(Point1.second - Point2.second);

    return distance;
}

double euclideanDist(pair<double, double> Point1, pair<double, double> Point2)
{
    double distance;
    double x, y;
    x = Point1.first - Point2.first;
    y = Point1.second - Point2.second;

    distance = pow(x, 2) + pow(y, 2);
    distance = sqrt(distance);

    return distance;
}

double euclideanDist(class Point* Point1, class Point* Point2)
{
    double distance = 0;

    unsigned int dimensions1 = Point1->getSize();

    double *coord_1 = Point1->getCoord();
    double *coord_2 = Point2->getCoord();

    for (int i = 0; i < dimensions1; i++)
    {
        distance += pow(coord_1[i] - coord_2[i], 2);
    }

    distance = sqrt(distance);

    return distance;
}