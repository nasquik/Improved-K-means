#include "metrics.hpp"
#include "dynamicTimeWarping.hpp"
#include "objectiveFunctions.hpp"

double minDistObjective(class Point* medoid, unordered_set<class Point*> points)
{
    double distanceTotal = 0;
    double distanceOnePoint;

    for(auto iter: points)
    {
        distanceOnePoint = manhattanDist(medoid, iter);
        distanceTotal += distanceOnePoint;
    }

    return distanceTotal;
}

double minDistObjective(class Curve* medoid, unordered_set<class Curve*> curves)
{
    double distanceTotal = 0;
    double distanceOnePoint;

    for(auto iter: curves)
    {

        distanceOnePoint = dtwDist(medoid, iter);
        distanceTotal += distanceOnePoint;
    }

    return distanceTotal;
}