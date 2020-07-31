#include <string>
#include <cmath>
#include <random>
#include <iostream>

#include "dataStructs.hpp"
#include "bruteForce.hpp"
#include "gridCurve.hpp"
#include "dynamicTimeWarping.hpp"
#include "LSH.hpp"
#include "metrics.hpp"

using namespace std;
////// GRID CURVES //////

template <class T>
gridCurve<T>::gridCurve(vector<class Curve *> *curves, int k, int L, int w, int probes, int minCurvePoints, int maxCurvePoints, double maxCoord)
{
    this->curves = curves;
    this->k = k;
    this->L = L;
    this->w = w;
    this->probes = probes;
    this->minCurvePoints = minCurvePoints;
    this->maxCurvePoints = maxCurvePoints;
    this->maxCoord = maxCoord * 1000;

    //use a multiplication of average distance of consecutive points for better results
    this->delta = 6 * calculateDelta();
    this->maxCurvePoints = maxCurvePoints * 2;

    this->displacement = new double *[L];

    //initialize the tau vectors
    for (int i = 0; i < L; i++)
    {
        this->displacement[i] = new double[2];
    }

    //create the tau vectors
    random_device rd;
    mt19937 generator(rd());
    uniform_real_distribution<> dis(0.0, this->delta);

    this->points = new vector<class Point *> *[L];

    //create the points produces from the snapped curves for each lsh structure

    for (int i = 0; i < L; i++)
    {
        points[i] = new vector<class Point *>;
        points[i]->reserve(curves->size());
    }

    for (int i = 0; i < this->L; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            this->displacement[i][j] = dis(generator);
        }

        for (int t = 0; t < curves->size(); t++)
        {
            points[i]->push_back(createVector(curves->at(t), i));
        }
    }
    cout << "Constructor done!" << endl;
}

template <>
int gridCurve<class LSH<class Curve *> >::initializeAlgorithm()
{

    this->algorithm = new class LSH<class Curve *> *[L];

    //get an estimation of the average distance of a set of points to use as w
    for (int i = 0; i < L; i++)
    {
        this->w = 15 * calculateW(points[i], 5);
        if(this->w == 0)
            this->w = 15;
        this->algorithm[i] = new class LSH<class Curve *>(this->k, 1, this->w, points[i]);
    }

    cout << "Initialized algorithm!" << endl;
    return 0;
}



template <class T>
int gridCurve<T>::initializeAlgorithm()
{
    return 0;
}

template <class T>
class Point *gridCurve<T>::gridCurve::createVector(class Curve *curve, int gridNo)
{

    //create the coord vector
    double *coord = new double[this->maxCurvePoints];

    //use the snapped value of the max coordinate as padding
    double snappedPadding1 = round((this->maxCoord - this->displacement[gridNo][0]) / this->delta);
    double snappedPadding2 = round((this->maxCoord - this->displacement[gridNo][1]) / this->delta);

    for (int i = 0; i < this->maxCurvePoints; i += 2)
    {
        coord[i] = snappedPadding1;
        coord[i + 1] = snappedPadding2;
    }

    double temp1 = 0, temp2 = 0;

    //snap the curves onto the grid
    int pos = 0;
    for (int i = 0; i < curve->getSize(); i++)
    {
        temp1 = (curve->getCoord()[i].first - this->displacement[gridNo][0]) / this->delta;
        temp2 = (curve->getCoord()[i].second - this->displacement[gridNo][1]) / this->delta;
        temp1 = round(temp1);
        temp2 = round(temp2);

        // remove duplicates
        if (i > 0)
        {
            if (coord[2 * pos] != temp1 || coord[(2 * pos) + 1] != temp2)
            {
                pos++;
                coord[2 * pos] = temp1;
                coord[(2 * pos) + 1] = temp2;
            }
        }
        else
        {
            coord[2 * i] = temp1;
            coord[(2 * i) + 1] = temp2;
        }
    }
    class Point *gridCurve = new class Point(curve->getID(), coord, this->maxCurvePoints, curve);

    return gridCurve;
}

template <>
class Curve *gridCurve<class LSH<class Curve *> >::findNN(class Curve *query, double *dist)
{
    class Point *queryGridCurve;
    class Curve *temp, *neighbor = NULL;
    double minDistance = numeric_limits<double>::max();

    //find the nearest curves among the LSH structures

    for (int i = 0; i < L; i++)
    {
        queryGridCurve = createVector(query, i);
        temp = algorithm[i]->findNN(queryGridCurve, dist);
        if (temp == NULL) //if no neighbor found
            continue;
        if (*dist < minDistance)
        {
            minDistance = *dist;
            neighbor = temp;
        }
    }

    cout << "Finished searching for nearest neigbor!" << endl;

    if (neighbor == NULL)
        return NULL;

    return neighbor;
}

template<class T>
void gridCurve<T>::assign(unordered_map<class Curve*,unordered_set<class Curve*> >*centroids,vector<class Curve*> *elements){
    unordered_map<class Curve*, class Curve*> pointToCentroid;

    //parse through the hashtables and assign centroids
    for (auto it: *centroids){
        for(int i = 0 ; i < L ;i++){
            unsigned int hashNo;
            int temp = 0;

            //generate a point to get the neighbors
            class Point * tempPoint = createVector(it.first,i);
            vector<class Point*> neighb = algorithm[i]->returnNeighb(tempPoint,&temp);
            delete tempPoint;
            if(!temp)
                continue;
            for (int j = 0; j < neighb.size();j++){
                if(neighb.at(j)->getCurvePtr()!=it.first){ //if this point is not our centroid
            
                    //if this place is not a centroid
                    if (centroids->find(neighb.at(j)->getCurvePtr())!=centroids->end()){
                        //if its already assigned to a centroid
                        if (pointToCentroid.find(neighb.at(j)->getCurvePtr())!=pointToCentroid.end()){
                            
                            //if the assigned centroid isn't the same as the one we examine
                            if(pointToCentroid.find(neighb.at(j)->getCurvePtr())->second != it.first){
                                //compare the distances and assign the closest
                                if(dtwDist(neighb.at(j)->getCurvePtr(),it.first)<dtwDist(neighb.at(j)->getCurvePtr(),(pointToCentroid.find(neighb.at(j)->getCurvePtr())->second))){
                                    //delete the previous assignment from both the pointToCentroid map and the centroid map
                                    //and assign it anew to both of them

                                    //delete from centroid's cluster
                                    centroids->find(pointToCentroid.find(neighb.at(j)->getCurvePtr())->second)->second.erase(neighb.at(j)->getCurvePtr());
                                    //delete from pointToCentroid
                                    pointToCentroid.erase(neighb.at(j)->getCurvePtr());

                                    //add to the new centroid's cluster
                                    it.second.insert(neighb.at(j)->getCurvePtr());
                                    //add to the pointToCentroid
                                    pointToCentroid.insert(make_pair(neighb.at(j)->getCurvePtr(),it.first));
                                    
                                }
                                //else do nothing it's a duplicate from another hash
                            }
                        }
                        else{ //if it's a new point just assign it
                            pointToCentroid.insert(make_pair(neighb.at(j)->getCurvePtr(),it.first));
                            it.second.insert(neighb.at(j)->getCurvePtr());
                        }
                    }
                }
            }
        
        }
    }

        //parse through the elements and assign the ones that are left
    for(int i = 0; i < elements->size();i++){

        //if its unassigned
        if (pointToCentroid.find(elements->at(i))== pointToCentroid.end()){
            double minDist = numeric_limits<double>::max();
            class Curve * closestCentroid = NULL;

            for(auto it: *centroids){
                double distance = dtwDist(elements->at(i),it.first);
                if (distance<minDist){
                    minDist = distance;
                    closestCentroid = it.first;
                }
            }
            centroids->find(closestCentroid)->second.insert(elements->at(i));
        }

    }
}


template <class T>
class Curve *gridCurve<T>::findNN(class Curve *query, double *distance)
{
    return NULL;
}

template <class T>
gridCurve<T>::~gridCurve()
{
    //delete the structures and free the alocated memory
    for (int i = 0; i < L; i++)
    {
        delete [] displacement[i];
        delete algorithm[i];
        while (!points[i]->empty())
        {
            delete points[i]->back();
            points[i]->pop_back();
        }
        delete points[i];
    }

    delete[] displacement;
    delete[] points;
    delete[] algorithm;

    cout << "Destructor Done!" << endl;
}

template <class T>
double gridCurve<T>::calculateDelta()
{
    //calculate the average distance of consecutive points of the curves
    double curveAvg, totalAvg = 0, tempDist;
    int total = curves->size();

    for (int i = 0; i < curves->size(); i++)
    {
        curveAvg = 0;
        int size = curves->at(i)->getSize();
        if (size < 2)
        {
            total--;
            continue;
        }
        pair<double, double> *coords = curves->at(i)->getCoord();

        for (int j = 1; j < size; j++)
        {
            curveAvg += euclideanDist(coords[j], coords[j - 1]);
        }
        curveAvg /= (double)(size - 1);
        totalAvg += curveAvg;
    }
    totalAvg /= (double)total;

    return totalAvg;
}
template class gridCurve<class LSH<class Curve *> >;
