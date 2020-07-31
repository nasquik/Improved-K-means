#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <iterator>
#include <limits>
#include <random>
#include <utility>
#include <fstream>

#include "LSH.hpp"
#include "dataStructs.hpp"
#include "clustering.hpp"
#include "dynamicTimeWarping.hpp"
#include "metrics.hpp"
#include "objectiveFunctions.hpp"
#include "gridCurve.hpp"

using namespace std;

template <class T>
clustering<T>::clustering(int clusterNo, int L, int k, vector<T *> *items)
{
    kClusters = clusterNo;
    elements = items;
    centroids = new unordered_map<T *, unordered_set<T *> >();
    centroidIndex = new T *[kClusters];
    nHT = L;
    nHF = k;
}

template <class T>
T *clustering<T>::binarySearch(vector<pair<T *, double> > index, double value)
{
    int lowerBound = 0;
    int upperBound = index.size() - 1;
    int half = upperBound / 2;

    //cut the vector in half and compare the value of the middle cell with your given value
    //iterate that in binary search manner 
    while (1)
    {

        //if we have concluded to a single element or found the value
        if (upperBound - lowerBound == 1 || index.at(lowerBound).second == value)
        {
            return index.at(lowerBound).first;
        }
        else if (index.at(upperBound).second == value) //if we have found the value
        {
            return index.at(upperBound).first;
        }

        if (index.at(half).second > value) //decide which portion of the dataset to pick
        {
            upperBound = half;
        }
        else if (index.at(half).second < value)
        {
            lowerBound = half;
        }
        else
        {
            return index.at(half).first;
        }

        half = (upperBound - lowerBound) / 2 + lowerBound; //get the new limits of the subarray
    }
}

template <class T>
void clustering<T>::initialize()
{

    unordered_set<T *> curr_centroids;

    // choose centroids in a uniformly random way

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, elements->size());

    int index = 0;

    for (int i = 0; i < kClusters; i++)
    {

        index = dis(gen);
        curr_centroids.insert(elements->at(index));
        centroidIndex[i] = elements->at(index);
    }

    for (auto itr = curr_centroids.begin(); itr != curr_centroids.end(); ++itr)
    {
        cout << "centroid: " << (*itr)->getID() << endl;
        centroids->insert(make_pair(*itr, unordered_set<T *>()));
    }
}

template <>
void clustering<class Point>::initializePlus()
{
    unordered_set<class Point *> curr_centroids;

    // choose centroids in a uniformly random way

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, elements->size());

    int index = dis(gen);

    curr_centroids.insert(elements->at(index));
    centroidIndex[0] = elements->at(index);

    while (curr_centroids.size() < kClusters)
    {
        //create a vector of pairs
        vector<pair<class Point *, double> > nonCentroids;
        nonCentroids.reserve(elements->size() - curr_centroids.size());
        double sum = 0;

        //iterate through the non-centroid elements
        for (int i = 0; i < elements->size(); i++)
        {
            //if this is a non centroid element
            if (curr_centroids.find(elements->at(i)) == curr_centroids.end())
            {
                //loop through the distances
                double min_distance = numeric_limits<double>::max();

                for (auto itr = curr_centroids.begin(); itr != curr_centroids.end(); ++itr)
                {
                    double dist = manhattanDist(*itr, elements->at(i));
                    if (dist < min_distance)
                        min_distance = dist;
                }

                nonCentroids.push_back(make_pair(elements->at(i), sum));
                sum += min_distance;
            }
        }

        // //select a random number in the range of 0 to max d_i
        std::uniform_real_distribution<double> unif(0, sum);

        double value = unif(gen);
        class Point *item = binarySearch(nonCentroids, value);
        centroidIndex[curr_centroids.size()] = item;
        curr_centroids.insert(item);
    }

    for (int i = 0; i < kClusters; i++)
    {
        cout << "centroid: " << centroidIndex[i]->getID() << endl;
        centroids->insert(make_pair(centroidIndex[i], unordered_set<class Point *>()));
    }
};

template <>
void clustering<class Curve>::initializePlus()
{
    unordered_set<class Curve *> curr_centroids;

    // choose centroids in a uniformly random way

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, elements->size());

    int index = dis(gen);

    curr_centroids.insert(elements->at(index));
    centroidIndex[0]=elements->at(index);

    while (curr_centroids.size() < kClusters)
    {

        //create a vector of pairs
        vector<pair<class Curve *, double> > nonCentroids;
        nonCentroids.reserve(elements->size() - curr_centroids.size());
        double sum = 0;

        //iterate through the non-centroid elements
        for (int i = 0; i < elements->size(); i++)
        {
            //if this is a non centroid element
            if (curr_centroids.find(elements->at(i)) == curr_centroids.end())
            {
                //loop through the distances
                double min_distance = numeric_limits<double>::max();

                for (auto itr = curr_centroids.begin(); itr != curr_centroids.end(); ++itr)
                {
                    double dist = 0;
                    dist = dtwDist(*itr, elements->at(i));
                    if (dist < min_distance)
                        min_distance = dist;
                }

                nonCentroids.push_back(make_pair(elements->at(i), sum));
                sum += min_distance;
            }
        }

        // //select a random number in the range of 0 to max d_i
        std::uniform_real_distribution<double> unif(0, sum);

        double value = unif(gen);
        class Curve *item = binarySearch(nonCentroids, value);
        centroidIndex[curr_centroids.size()] = item;
        curr_centroids.insert(item);
    }

    for (auto itr = curr_centroids.begin(); itr != curr_centroids.end(); ++itr)
    {
        cout << "centroid: " << (*itr)->getID() << endl;
        centroids->insert(make_pair(*itr, unordered_set<class Curve *>()));
    }
};

template <class T>
void clustering<T>::initializePlus(){};

template <>
void clustering<class Point>::rangeSearchAssign()
{

    // create an lsh structure
    class LSH<class Point *> *lsh = new class LSH<class Point *>(nHF, nHT, 3600, elements);

    // assign points to clusters using the lsh class
    lsh->assign(centroids, elements);
    
    for (auto it : *centroids)
    {
        cout << "Centroid: " << it.first->getID() << endl;
        //      << "Points in Cluster: " << endl;

        // for (class Point *x : it.second)
        // {
        //     cout << x->getID() << endl;
        // }
    }

    delete lsh;
}

template <>
void clustering<class Curve>::rangeSearchAssign()
{
    unordered_map<class Curve *, class Curve *> curveToCentroid;

    //create and lsh structure
    class gridCurve<class LSH<class Curve *> > *lshCurve = new class gridCurve<class LSH<class Curve *> >(elements, nHF, nHT, 3600, 5, 10, 10, 50000);
    lshCurve->initializeAlgorithm();
    
    // assign points to clusters using the lsh class
    lshCurve->assign(centroids, elements);

    for (auto it : *centroids)
    {
        cout << "Centroid: " << it.first->getID() << endl;
            //  << "Curves in Cluster: " << endl;

        // for (class Curve *x : it.second)
        // {
        //     cout << x->getID() << endl;
        // }
    }

    delete lshCurve;
}

template <class T>
void clustering<T>::rangeSearchAssign(){};

template <>
void clustering<class Point>::lloydAssign()
{

    double minDist, tempDist;
    class Point *point;
    int index;
    class Point *closest;

    // iterate through all elements
    for (int i = 0; i < elements->size(); i++)
    {
        point = elements->at(i);

        // if point is a centroid, skip it
        if (centroids->find(point) != centroids->end())
            continue;

        // else calculate manhattan distance from all centroids and assign it to the closest
        minDist = numeric_limits<int>::max();

        for (auto it : *centroids)
        {
            tempDist = manhattanDist(point, it.first);
            if (tempDist < minDist)
            {
                minDist = tempDist;
                closest = it.first;
            }
        }
        centroids->at(closest).insert(point);
    }

    for (auto it : *centroids)
    {
        cout << "Centroid: " << it.first->getID() << endl;
        //      << "Points in Cluster: " << endl;

        // for (class Point *x : it.second)
        // {
        //     cout << x->getID() << endl;
        // }
    }
}

template <>
void clustering<class Curve>::lloydAssign()
{

    double minDist, tempDist;
    class Curve *curve;
    int index;
    class Curve *closest;

    // iterate through all elements
    for (int i = 0; i < elements->size(); i++)
    {
        curve = elements->at(i);

        // if point is a centroid, skip it
        if (centroids->find(curve) != centroids->end())
            continue;

        // else calculate manhattan distance from all centroids and assign it to the closest
        minDist = numeric_limits<int>::max();

        for (auto it : *centroids)
        {
            tempDist = dtwDist(curve, it.first);
            if (tempDist < minDist)
            {
                minDist = tempDist;
                closest = it.first;
            }
        }
        centroids->at(closest).insert(curve);
    }

    for (auto it : *centroids)
    {
        cout << "Centroid: " << it.first->getID() << endl;
        //      << "Curves in Cluster: " << endl;

        // for (class Curve *x : it.second)
        // {
        //     cout << x->getID() << endl;
        // }
    }
}

template <class T>
void clustering<T>::lloydAssign()
{}

template <>
int clustering<class Point>::meanVectorUpdate()
{

    //create a new unordered map for the updated centroids
    unordered_map<class Point *, unordered_set<class Point *> > *newMap = new unordered_map<class Point *, unordered_set<class Point *> >();
    //get the dimensions of each point
    int dimensions = elements->at(0)->getSize();
    class Point **newIndex = new class Point *[kClusters];

    int count = 1;
    //create the new mean vectors
    for (int j = 0; j < kClusters; j++)
    {
        auto it = centroids->find(centroidIndex[j]);
        //check the cluster
        double *temp = new double[dimensions];
        for (int i = 0; i < dimensions; i++)
        {
            temp[i] = 0;
        }

        for (auto it2 : it->second)
        {
            double *resulting = it2->getCoord();
            for (int i = 0; i < dimensions; i++)
            {
                temp[i] += resulting[i];
            }
        }

        //now calulated the average
        if ((it)->second.size() > 0)
        {
            for (int i = 0; i < dimensions; i++)
            {
                temp[i] /= it->second.size();
            }
        }

        //create the new vector
        class Point *newPoint = new class Point("fake" + to_string(count), temp, dimensions, 1);
        newMap->insert(make_pair(newPoint, unordered_set<class Point *>()));
        newIndex[j] = newPoint;
        count++;
    }

    // convergence test

    int flag = 1;

    for (int i = 0; i < kClusters; i++)
    {
        if (manhattanDist(centroidIndex[i], newIndex[i]) != 0)
        {
            flag = 0;
            break;
        }
    }

    //delete new structures
    if (flag)
    {
        for (auto it : *newMap)
        {
            if (it.first->Centroid() == 1)
                delete it.first;
        }

        delete newMap;
        delete[] newIndex;
        return 1;
    }

    for (auto it : *centroids)
    {
        if (it.first->Centroid() == 1)
            delete it.first;
    }

    delete centroids;
    delete[] centroidIndex;
    centroids = newMap;
    centroidIndex = newIndex;

    for (auto it : *centroids)
    {
        cout << "Centroid: " << it.first->getID() << endl;
    }
    cout << endl;

    return 0;
}

template <class T>
int clustering<T>::meanVectorUpdate()
{
    return 0;
}

template <>
int clustering<class Point>::pamUpdate()
{
    unordered_map<class Point *, unordered_set<class Point *> > *newClusters = new unordered_map<class Point *, unordered_set<class Point *> >();
    class Point **newIndex = new class Point *[kClusters];

    // iterate through all clusters
    for (int i = 0; i < kClusters; i++)
    {
        auto iterMap = centroids->find(centroidIndex[i]);

        class Point *currCentroid = iterMap->first;
        class Point *oldCentroid = currCentroid;
        unordered_set<class Point *> *nonMedoids = &iterMap->second;

        // calculate objective function for current centroid
        double centroidMinDist = minDistObjective(currCentroid, *nonMedoids);
        double potentialMinDist;

        // calculate objective function for all non-medoid points in cluster
        for (auto potentialCentroid : *nonMedoids)
        {
            potentialMinDist = minDistObjective(potentialCentroid, *nonMedoids);

            if (potentialMinDist < centroidMinDist)
            {
                currCentroid = potentialCentroid;
                centroidMinDist = potentialMinDist;
            }
        }

        // make the point with the min objective function the new centroid
        newClusters->insert(make_pair(currCentroid, unordered_set<class Point *>()));
        newIndex[i] = currCentroid;
        cout << "Old centroid ID: " << oldCentroid->getID() << endl;
        cout << "New centroid ID: " << currCentroid->getID() << endl;
    }

    // convergence test

    int flag = 1;

    for (int i = 0; i < kClusters; i++)
    {
        if (manhattanDist(centroidIndex[i], newIndex[i]) != 0)
        {
            flag = 0;
            break;
        }
    }

    //delete new structures
    if (flag)
    {
        delete[] newIndex;
        delete newClusters;
        return 1;
    }

    delete centroids;
    delete[] centroidIndex;
    centroids = newClusters;
    centroidIndex = newIndex;

    return 0;
}

template <>
int clustering<class Curve>::pamUpdate()
{
    unordered_map<class Curve *, unordered_set<class Curve *> > *newClusters = new unordered_map<class Curve *, unordered_set<class Curve *> >();
    class Curve **newIndex = new class Curve *[kClusters];

    // iterate through all clusters
    for (int i = 0; i < kClusters; i++)
    {
        auto iterMap = centroids->find(centroidIndex[i]);

        class Curve *currCentroid = iterMap->first;
        class Curve *oldCentroid = currCentroid;
        unordered_set<class Curve *> *nonMedoids = &iterMap->second;

        // calculate objective function for current centroid
        double centroidMinDist = minDistObjective(currCentroid, *nonMedoids);
        double potentialMinDist;

        // calculate objective function for all non-medoid points in cluster
        for (auto potentialCentroid : *nonMedoids)
        {
            potentialMinDist = minDistObjective(potentialCentroid, *nonMedoids);

            if (potentialMinDist < centroidMinDist)
            {
                currCentroid = potentialCentroid;
                centroidMinDist = potentialMinDist;
            }
        }

        // make the point with the min objective function the new centroid
        newClusters->insert(make_pair(currCentroid, unordered_set<class Curve *>()));
        newIndex[i] = currCentroid;
        cout << "Old centroid ID: " << oldCentroid->getID() << endl;
        cout << "New centroid ID: " << currCentroid->getID() << endl;
    }

    // convergence test

    int flag = 1;

    for (int i = 0; i < kClusters; i++)
    {
        if (dtwDist(centroidIndex[i], newIndex[i]) != 0)
        {
            flag = 0;
            break;
        }
    }

    //delete new structures
    if (flag)
    {
        delete[] newIndex;
        delete newClusters;
        return 1;
    }

    delete centroids;
    delete[] centroidIndex;
    centroids = newClusters;
    centroidIndex = newIndex;

    return 0;
}

template <class T>
int clustering<T>::pamUpdate()
{
    return 0;
}

template <>
int clustering<class Curve>::dbaUpdate()
{
    int lambda;
    unordered_map<class Curve *, unordered_set<class Curve *> > *newClusters = new unordered_map<class Curve *, unordered_set<class Curve *> >();
    class Curve **newIndex = new class Curve *[kClusters];

    class Curve *seqC;
    class Curve *newCurve;
    vector<pair<double, double> > *A;
    pair<double, double>* coords;
    pair<double, double> *coordsC;
    vector<pair<int, int> > *pairs;
    double sumX, sumY, distance;
    random_device rd;
    mt19937 gen(rd());
    int index, window;

    int count;
    int extCount = 0;

    // iterate through all clusters
    for (int i = 0; i < kClusters; i++)
    {
        auto iter = centroids->find(centroidIndex[i]);

        extCount++;

        // find mean length of all curves in cluster = lambda
        lambda = 0;
        for (auto clusterCurve : iter->second)
        {
            lambda += clusterCurve->getSize();
        }

        lambda /= iter->second.size();
        
        coords = new pair<double, double>[lambda];

        // choose a curve(sequence) s0 from the cluster in a uniformly random way
        uniform_int_distribution<> dis(0, iter->second.size() - 1);

        // choose a lambda-length subsequence of s0 in a uniformly random way
        auto s0 = begin(iter->second);

        do
        {
            index = dis(gen);
            s0 = begin(iter->second);
            advance(s0, index);
        } while ((*s0)->getSize() < lambda);

        window = (*s0)->getSize() - lambda;
        uniform_int_distribution<> dist(0, window);
        index = dist(gen);

        for (int i = index; i < (index + lambda); i++)
        {
            coords[i-index] = (*s0)->getCoord()[i];
        } 

        // set that subsequence as seqC
        count = 0;
        seqC = new class Curve("fake_" + to_string(extCount) + "_0", coords, lambda, 1);
        newCurve = NULL;

        // until the 2 sequences converge or there have been enough iterations
        do
        {
            if (newCurve != NULL)
            {
                delete seqC;
                seqC = newCurve;
            }

            // create vector A of size lambda
            A = new vector<pair<double, double> >[lambda];

            // iterate through cluster
            for (auto seqS : iter->second)
            {
                // find index-pairs od best traversal
                pairs = dtwBestTraversal(seqC, seqS);
                int P1, P2;
                pair<double, double> coordsS;

                for (int i = 0; i < pairs->size(); i++)
                {
                    P1 = pairs->at(i).first;
                    P2 = pairs->at(i).second;
                    A[P1].push_back(seqS->getCoord()[P2]);
                }

                delete pairs;
            }

            // calculate coordinates of mean curve
            coordsC = new pair<double, double>[lambda];

            for (int j = 0; j < lambda; j++)
            {
                sumX = 0.0;
                sumY = 0.0;
                pair<double, double> coord;

                for (int t = 0; t < A[j].size(); t++)
                {
                    sumX += A[j].at(t).first;
                    sumY += A[j].at(t).second;
                }

                if (A[j].size() == 0)
                {
                    coord.first = 0.0;
                    coord.second = 0.0;
                }
                else
                {
                    coord.first = sumX / A[j].size();
                    coord.second = sumY / A[j].size();
                }

                coordsC[j] = coord;
            }

            count++;

            // create new curve that reflects everything calculated above
            newCurve = new class Curve("fake_" + to_string(extCount) + "_" + to_string(count), coordsC, lambda, 1);

            // calculate dtw distance between seqC and the new curve
            distance = dtwDist(seqC, newCurve);
            //cout << "Distance: " << distance << endl;

            delete[] A;

        } while (distance > 0.002 && count < 1000);

        // set mean curve as the new centroid
        newClusters->insert(make_pair(newCurve, unordered_set<class Curve *>()));
        newIndex[i] = newCurve;
        delete seqC;
    }

    // convergence test

    int flag = 1;

    for (int i = 0; i < kClusters; i++)
    {
        if (dtwDist(centroidIndex[i], newIndex[i]) != 0)
        {
            flag = 0;
            break;
        }
    }

    //delete new structures
    if (flag)
    {
        for (auto it : *newClusters)
        {
            if (it.first->Centroid() == 1)
                delete it.first;
        }
        delete[] newIndex;
        return 1;
    }

    for (auto it : *centroids)
    {
        if (it.first->Centroid() == 1)
            delete it.first;
    }

    delete centroids;
    delete[] centroidIndex;
    centroids = newClusters;
    centroidIndex = newIndex;

    return 0;
}

template <class T>
int clustering<T>::dbaUpdate()
{
    return 0;
}

template <>
vector<double> *clustering<class Point>::silhouetteEvaluate()
{

    vector<double> *clusterEvaluations = new vector<double>();
    vector<double> *pointEvaluations;

    int currClusterSize, secClusterSize, overallSize = 0;
    class Point *currCentroid, *tempCentroid, *secCentroid;
    unordered_set<class Point *> *secCluster;
    double minDistCentroid, tempDistCentroid;
    double averageDistCurr, averageDistSec;
    double silhouetteEvaluationPoint, silhouetteEvaluationCluster;
    double overallEvaluation = 0.0;

    // iterate through all clusters
    for (int i = 0; i < kClusters; i++)
    {
        auto iter = centroids->find(centroidIndex[i]);

        currCentroid = iter->first;
        currClusterSize = (iter->second).size();
        overallSize += currClusterSize;
        pointEvaluations = new vector<double>();
        silhouetteEvaluationCluster = 0.0;

        // for every point in the cluster, calculate
        for (auto currClusterPoint : iter->second)
        {
            averageDistCurr = 0.0;
            averageDistSec = 0.0;
            minDistCentroid = numeric_limits<double>::max();

            // average distance from points in the same cluster a(i)

            for (auto pointInCluster : iter->second)
            {
                averageDistCurr += manhattanDist(currClusterPoint, pointInCluster);
            }

            averageDistCurr = averageDistCurr / currClusterSize;

            // 2nd closest centroid

            for (auto tempCluster : *centroids)
            {
                tempCentroid = tempCluster.first;

                if (tempCentroid->getID() != currCentroid->getID())
                {

                    tempDistCentroid = manhattanDist(currClusterPoint, tempCentroid);
                    if (tempDistCentroid < minDistCentroid)
                    {
                        minDistCentroid = tempDistCentroid;
                        secCentroid = tempCentroid;
                    }
                }
            }

            // average distance from points in the second closest cluster b(i);

            auto secCluster = centroids->find(secCentroid);
            secClusterSize = (secCluster->second).size();
            for (auto pointInSecCluster : secCluster->second)
            {
                averageDistSec += manhattanDist(currClusterPoint, pointInSecCluster);
            }
            averageDistSec = averageDistSec / secClusterSize;

            // slhouette evaluation

            silhouetteEvaluationPoint = (averageDistSec - averageDistCurr) / max(averageDistSec, averageDistCurr);
            pointEvaluations->push_back(silhouetteEvaluationPoint);
            overallEvaluation += silhouetteEvaluationPoint;
        }

        // silhouette evaluation for each cluster
        for (int i = 0; i < currClusterSize; i++)
        {
            silhouetteEvaluationCluster += pointEvaluations->at(i);
        }

        delete pointEvaluations;
        silhouetteEvaluationCluster = silhouetteEvaluationCluster / currClusterSize;
        clusterEvaluations->push_back(silhouetteEvaluationCluster);
    }

    // overall evaluation
    overallEvaluation = overallEvaluation / overallSize;
    clusterEvaluations->push_back(overallEvaluation);
    return clusterEvaluations;
}

template <>
vector<double> *clustering<class Curve>::silhouetteEvaluate()
{

    vector<double> *clusterEvaluations = new vector<double>;
    vector<double> *curveEvaluations;

    int currClusterSize, secClusterSize, overallSize = 0;
    class Curve *currCentroid, *tempCentroid, *secCentroid;
    unordered_set<class Curve *> *secCluster;
    double minDistCentroid, tempDistCentroid;
    double averageDistCurr, averageDistSec;
    double silhouetteEvaluationCurve, silhouetteEvaluationCluster;
    double overallEvaluation = 0.0;

    // for every point in the cluster, calculate
    for (int i = 0; i < kClusters; i++)
    {
        auto iter = centroids->find(centroidIndex[i]);

        currCentroid = iter->first;
        currClusterSize = (iter->second).size();
        overallSize += currClusterSize;
        curveEvaluations = new vector<double>;
        silhouetteEvaluationCluster = 0.0;

        for (auto currClusterCurve : iter->second)
        {
            averageDistCurr = 0.0;
            averageDistSec = 0.0;
            minDistCentroid = numeric_limits<double>::max();

            // average distance from curves in same cluster a(i)

            for (auto curveInCluster : iter->second)
            {
                averageDistCurr += dtwDist(currClusterCurve, curveInCluster);
            }

            averageDistCurr = averageDistCurr / currClusterSize;

            // second closest centroid

            for (auto tempCluster : *centroids)
            {
                tempCentroid = tempCluster.first;

                if (tempCentroid->getID() != currCentroid->getID())
                {

                    tempDistCentroid = dtwDist(currClusterCurve, tempCentroid);
                    if (tempDistCentroid < minDistCentroid)
                    {
                        minDistCentroid = tempDistCentroid;
                        secCentroid = tempCentroid;
                    }
                }
            }

            // average distance from curves in second closest cluster b(i);

            auto secCluster = centroids->find(secCentroid);
            secClusterSize = (secCluster->second).size();
            for (auto curveInSecCluster : secCluster->second)
            {
                averageDistSec += dtwDist(currClusterCurve, curveInSecCluster);
            }
            averageDistSec = averageDistSec / secClusterSize;

            // slhouette evaluation

            silhouetteEvaluationCurve = (averageDistSec - averageDistCurr) / max(averageDistSec, averageDistCurr);
            curveEvaluations->push_back(silhouetteEvaluationCurve);
            overallEvaluation += silhouetteEvaluationCurve;
        }

        // silhouette evaluation for each cluster
        for (int i = 0; i < currClusterSize; i++)
        {
            silhouetteEvaluationCluster += curveEvaluations->at(i);
        }

        delete curveEvaluations;
        silhouetteEvaluationCluster = silhouetteEvaluationCluster / currClusterSize;
        clusterEvaluations->push_back(silhouetteEvaluationCluster);
    }

    // overall evaluation
    overallEvaluation = overallEvaluation / overallSize;
    clusterEvaluations->push_back(overallEvaluation);
    return clusterEvaluations;
}

template <class T>
vector<double> *clustering<T>::silhouetteEvaluate()
{
    return NULL;
}

template <class T>
clustering<T>::~clustering()
{

    for (auto it : *centroids)
    {
        if (it.first->Centroid() == 1)
            delete it.first;
    }

    delete centroids;
    delete[] centroidIndex;
}

template <>
void clustering<class Point>::printResults(string fileName, string I, string A, string U, double clusteringTime, char isComplete)
{
    string algo = I + "x" + A + "x" + U + "x";

    //open file to write into
    ofstream myfile;
    myfile.open(fileName);
    myfile << "Algorithm: " << algo << endl;

    //start printing info about each cluster
    int count = 1;
    for (auto it : *centroids)
    {
        myfile << "CLUSTER-" << count << " "
               << "{ size: " << it.second.size() << ", centroid: ";
        count++;
        myfile << it.first->getID() << " }" << endl;
    }
    myfile << "clustering_time: " << clusteringTime << endl;

    //print sihlouette
    vector<double> *result = silhouetteEvaluate();
    myfile << "Silhouette: [";
    for (int j = 0; j < result->size(); j++)
    {
        myfile << result->at(j);
        if (j < result->size() - 1)
        {
            myfile << ", ";
        }
    }
    myfile << "]" << endl;

    delete result;
    if (isComplete)
    {
        count = 1;
        for (auto it : *centroids)
        {
            myfile << "CLUSTER-" << count << " (";

            int printed = 0;
            int size = it.second.size();
            for (auto it2 : it.second)
            {
                myfile << it2->getID();
                printed++;
                if (printed != size)
                {
                    myfile << ", ";
                }
            }
            myfile << ")" << endl;
            count++;
        }
    }
    myfile.close();
}

template <>
void clustering<class Curve>::printResults(string fileName, string I, string A, string U, double clusteringTime, char isComplete)
{
    string algo = I + "x" + A + "x" + U + "x";

    //open file to write into
    ofstream myfile;
    myfile.open(fileName);
    myfile << "Algorithm: " << algo << endl;

    //start printing info about each cluster
    int count = 1;
    for (auto it : *centroids)
    {
        myfile << "CLUSTER-" << count << " "
               << "{ size: " << it.second.size() << ", centroid: ";
        count++;
        //check if the cluster is not a valid curve
        if (U == "M")
        {
            myfile << "[ ";
            pair<double, double> *coords = it.first->getCoord();
            for (int i = 0; i < it.first->getSize(); i++)
            {
                myfile << "(" << coords[i].first << ", " << coords[i].second << ")";
                myfile << " ";
            }
            myfile << "] }" << endl;
        }
        else
        {
            myfile << it.first->getID() << " }" << endl;
        }
    }
    myfile << "clustering_time: " << clusteringTime << endl;

    //print sihlouette
    vector<double> *result = silhouetteEvaluate();
    myfile << "Silhouette: [";
    for (int j = 0; j < result->size(); j++)
    {
        myfile << result->at(j);
        if (j < result->size() - 1)
        {
            myfile << ", ";
        }
    }
    myfile << "]" << endl;

    delete result;

    if (isComplete)
    {
        count = 1;
        for (auto it : *centroids)
        {
            myfile << "CLUSTER-" << count << " (";

            int printed = 0;
            int size = it.second.size();
            for (auto it2 : it.second)
            {
                myfile << it2->getID();
                printed++;
                if (printed != size)
                {
                    myfile << ", ";
                }
            }
            myfile << ")" << endl;
            count++;
        }
    }
    myfile.close();
}

template <class T>
void clustering<T>::printResults(string fileName, string I, string A, string U, double clusteringTime, char isComplete)
{
}

template class clustering<class Point>;
template class clustering<class Curve>;
