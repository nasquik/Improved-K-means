#include <iostream>
#include <fstream>
#include <random>
#include <string>

#include "dataStructs.hpp"
#include "metrics.hpp"
#include "dynamicTimeWarping.hpp"
#include "bruteForce.hpp"

using namespace std;

pair<class Point *, double> *bruteForce(vector<class Point *> *points, class Point *query)
{
    //find the nearest neighbor from the dataset of a given query
    size_t npoints = points->size();

    pair<class Point *, double> *nearestNeighbor = new pair<class Point *, double>;

    class Point *closestPoint = points->at(0);
    double minDistance = manhattanDist(points->at(0), query);

    class Point *tempPoint;
    double tempDist;

    for (int i = 1; i < npoints; i++)
    {
        tempPoint = points->at(i);
        tempDist = manhattanDist(tempPoint, query);
        if (tempDist < minDistance)
        {
            minDistance = tempDist;
            closestPoint = tempPoint;
        }
    }
    nearestNeighbor->first = closestPoint;
    nearestNeighbor->second = minDistance;

    return nearestNeighbor;
}

int bruteForceAll(vector<class Point *> *points, vector<class Point *> *queries)
{
    //same as brute force but operates with a set of queries
    size_t npoints = points->size();
    size_t nqueries = queries->size();

    ofstream outfile;
    outfile.open("brute.dat");
    double avg = 0;

    //for each nqueries find the one with the smallest manhattan distance
    for (int i = 0; i < nqueries; i++)
    {
        double temp_dist, min_distance = manhattanDist(points->at(0), queries->at(i));
        class Point *closest_point = points->at(0);
        for (int j = 1; j < npoints; j++)
        {
            temp_dist = manhattanDist(points->at(j), queries->at(i));
            if (temp_dist < min_distance)
            {
                min_distance = temp_dist;
                closest_point = points->at(j);
            }
        }
        //write the result to the file
        outfile << "Query Point: " << queries->at(i)->getID() << "\tNearest Neighbor: " << closest_point->getID() << "\tDistance: " << min_distance << endl;
        avg += min_distance;
    }
    avg /= nqueries;
    cout << "Average Distance: " << avg << endl;
    outfile.close();
    return 1;
}

pair<class Point *, double> *bruteForceTweaked(vector<class Point *> *points, class Point *query, string ID)
{
    //same as bruteforce but ignores points with the same id as the given query
    size_t npoints = points->size();

    pair<class Point *, double> *nearestNeighbor = new pair<class Point *, double>;

    class Point *closestPoint = points->at(0);
    double minDistance = manhattanDist(points->at(0), query);

    class Point *tempPoint;
    double tempDist;

    for (int i = 1; i < npoints; i++)
    {
        tempPoint = points->at(i);
        if (ID == tempPoint->getID())
            continue;
        tempDist = manhattanDist(tempPoint, query);
        if (tempDist < minDistance)
        {
            minDistance = tempDist;
            closestPoint = tempPoint;
        }
    }
    nearestNeighbor->first = closestPoint;
    nearestNeighbor->second = minDistance;

    return nearestNeighbor;
}

pair<class Curve *, double> *bruteForceCurve(vector<class Curve *> *curves, class Curve *query)
{
    //finds the nearest neighbor of a curve with dtw distance
    size_t npoints = curves->size();

    pair<class Curve *, double> *nearestNeighbor = new pair<class Curve *, double>;

    class Curve *closestCurve = curves->at(0);
    double minDistance = dtwDist(curves->at(0), query);

    class Curve *tempCurve;
    double tempDist;

    for (int i = 1; i < npoints; i++)
    {
        tempCurve = curves->at(i);
        tempDist = dtwDist(tempCurve, query);
        if (tempDist < minDistance)
        {
            minDistance = tempDist;
            closestCurve = tempCurve;
        }
    }
    nearestNeighbor->first = closestCurve;
    nearestNeighbor->second = minDistance;

    return nearestNeighbor;
}

int calculateW(vector<class Point *> *points, int percentage)
{
    //estimates the value of w used in the algorithms
    //uses a percentage of random points from the datasets as queries
    //and returns their average distance from their nearest neighbor
    
    int toCheck = (points->size() * percentage) / 100;

    random_device rd;
    mt19937 gen(rd());
    pair<class Point *, double> *temp;
    double avgDist = 0;

    uniform_int_distribution<> dis(0, points->size() - 1);

    for (int i = 0; i < toCheck; i++)
    {
        class Point *temp_point = points->at(dis(gen)); //generate a random point
        temp = bruteForceTweaked(points, temp_point, temp_point->getID());
        avgDist += temp->second;
        delete temp;
    }

    avgDist = ceil(avgDist / toCheck);
    cout << "ESTIMATED W IS : " << (int)avgDist << endl;
    return (int)avgDist;
}