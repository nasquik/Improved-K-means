#include <iostream>
#include <limits>
#include <cstring>
#include <cmath>
#include <unordered_set>

#include "dataStructs.hpp"
#include "metrics.hpp"
#include "dynamicTimeWarping.hpp"
#include "hashTable.hpp"
#include "LSH.hpp"

using namespace std;

// Class used for the implementation of the LSH algorithm
template <class T>
LSH<T>::LSH(int k, int L, int w, vector<class Point *> *input)
{
    this->k = k;
    this->L = L;
    this->w = w;

    this->input = input;

    int dim = input->at(0)->getSize();

    unsigned int hashtableSize = input->size() / 8;
    hashTables = new class HashTable<vector<class Point *> >[L];
    for (int i = 0; i < L; i++)
    {
        hashTables[i].initialize(hashtableSize, w, k, dim, 0);
    }
    //insert the points to the hashtables
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < input->size(); j++)
        {
            hashTables[i].insertPoint(input->at(j));
        }
        hashTables[i].printHashTable();
    }
}

template <class T>
LSH<T>::~LSH()
{
    delete[] hashTables;
}

template <class T>
int LSH<T>::getk()
{
    return k;
}

template <class T>
int LSH<T>::getL()
{
    return L;
}

template <class T>
int LSH<T>::getw()
{
    return w;
}

template <>
class Point *LSH<class Point *>::findNN(class Point *query, double *dist)
{

    //find the nearest neighbors
    class Point *b = NULL;
    class Point *p = NULL;
    double manhattanD;
    double distance = numeric_limits<double>::max();
    vector<class Point* > neighbors;
    int count;
    int temp;
    for (int i = 0; i < L; i++)
    {
        unsigned int amplifiedResult = 0;
        neighbors = hashTables[i].getneighbors(query, &amplifiedResult,&temp);
        count = 0;

        for (int j = 0; j < neighbors.size(); j++)
        {
                count++;
                if (count > 10 * L)
                    break;

                p = neighbors.at(j);
                manhattanD = manhattanDist(query, p);

                if (manhattanD < distance)
                {
                    b = p;
                    distance = manhattanD;
                }
            
        }
    }
    *dist = distance;

    return b;
}

template <>
class Curve *LSH<class Curve *>::findNN(class Point *query, double *dist)
{
    class Point *b = NULL;
    class Point *p = NULL;
    double dtwD;
    double distance = numeric_limits<double>::max();
    vector<class Point * > neighbors;
    int count;
    int temp;
    cout << query->getCurvePtr()->getID() << endl;

    for (int i = 0; i < L; i++)
    {
        unsigned int amplifiedResult = 0;
        neighbors = hashTables[i].getneighbors(query, &amplifiedResult,&temp);
        count = 0;

        for (int j = 0; j < neighbors.size(); j++)
        {

                count++;
                if (count > 30 * L)
                    break;

                p = neighbors.at(j);
                dtwD = dtwDist(query->getCurvePtr(), p->getCurvePtr());

                if (dtwD < distance)
                {
                    b = p;
                    distance = dtwD;
                }
            
        }
    }

    *dist = distance;
    if (b != NULL)
        return b->getCurvePtr();
    else
        return NULL;
}

template <class T>
T LSH<T>::findNN(class Point *query, double *dist) {}

// template <>
// vector<pair<class Point *, double> > *LSH<class Point *>::findRadiusNN(class Point *query, double radius)
// {

//     class Point *p = NULL;
//     double manhattanD;
//     vector<pair<class Point *, double> > *radiusNaighbors = new vector<pair<class Point *, double> >;
//     vector<pair<class Point *, unsigned int> > neighbors;
//     int count;

//     for (int i = 0; i < L; i++)
//     {
//         unsigned int amplifiedResult = 0;
//         neighbors = hashTables[i].getneighbors(query, &amplifiedResult);
//         count = 0;

//         for (int j = 0; j < neighbors.size(); j++)
//         {
//             if (amplifiedResult == neighbors.at(j).second)
//             {
//                 count++;
//                 if (count > 10 * L)
//                     break;

//                 p = neighbors.at(j).first;
//                 manhattanD = manhattanDist(query, p);

//                 if (manhattanD <= radius)
//                 {
//                     pair<class Point *, double> neighborRadius;
//                     neighborRadius.first = p;
//                     neighborRadius.second = manhattanD;
//                     radiusNaighbors->push_back(neighborRadius);
//                 }
//             }
//         }
//     }

//     if (!radiusNaighbors->empty())
//         return radiusNaighbors;
//     else
//     {
//         delete radiusNaighbors;
//         return NULL;
//     }
// }

template <class T>
vector<pair<class Point *, double> > *LSH<T>::findRadiusNN(class Point *query, double radius)
{
    return NULL;
}

//////////////////////////////////////////////////////////////////////////
template<>
void LSH<class Point*>::assign(unordered_map<class Point*,unordered_set<class Point*> >*centroids,vector<class Point*> *elements){

    unordered_map<class Point*, class Point*> pointToCentroid;

    //parse through the hashtables and assign the centroids
    for (auto it: *centroids){
        for (int i = 0; i < L;i++){
            unsigned int hashNo;
            int temp = 0;
            vector <class Point *> neighb = hashTables[i].getneighbors(it.first,&hashNo,&temp);
            if (!temp)
                continue;
            //cout << "SIZE" << neighb.size()<<endl;
            for (int j = 0; j < neighb.size();j++){
                if(neighb.at(j)!=it.first){ //if this point is not our centroid
            
                    //if this place is not a centroid
                    if (centroids->find(neighb.at(j))!=centroids->end()){
                        //if its already assigned to a centroid
                        if (pointToCentroid.find(neighb.at(j))!=pointToCentroid.end()){
                            
                            //if the assigned centroid isn't the same as the one we examine
                            if(pointToCentroid.find(neighb.at(j))->second != it.first){
                                //compare the distances and assign the closest
                                if(euclideanDist(neighb.at(j),it.first)<euclideanDist(neighb.at(j),(pointToCentroid.find(neighb.at(j))->second))){
                                    //delete the previous assignment from both the pointToCentroid map and the centroid map
                                    //and assign it anew to both of them

                                    //delete from centroid's cluster
                                    centroids->find(pointToCentroid.find(neighb.at(j))->second)->second.erase(neighb.at(j));
                                    //delete from pointToCentroid
                                    pointToCentroid.erase(neighb.at(j));

                                    //add to the new centroid's cluster
                                    it.second.insert(neighb.at(j));
                                    //add to the pointToCentroid
                                    pointToCentroid.insert(make_pair(neighb.at(j),it.first));
                                    
                                }
                                //else do nothing it's a duplicate from another hash
                            }
                        }
                        else{ //if it's a new point just assign it
                            pointToCentroid.insert(make_pair(neighb.at(j),it.first));
                            it.second.insert(neighb.at(j));
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
            class Point * closestCentroid = NULL;

            for(auto it: *centroids){
                double distance = euclideanDist(elements->at(i),it.first);
                if (distance<minDist){
                    minDist = distance;
                    closestCentroid = it.first;
                }
            }
            centroids->find(closestCentroid)->second.insert(elements->at(i));
        }

    }

};


template <class T>
void LSH<T>::assign(unordered_map<T,unordered_set<T> >*,vector<T>*){};


template <>
vector<class Point*> & LSH<class Curve*>::returnNeighb(class Point* givenPoint,int * valid){
    unsigned int temp;


    return hashTables[0].getneighbors(givenPoint,&temp,valid);
}

template<class T>
vector<class Point*> &LSH<T>::returnNeighb(class Point* givenPoint,int * valid)
{};

template class LSH<class Curve *>;
template class LSH<class Point *>;