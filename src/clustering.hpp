#ifndef CLUSTERING_HPP
#define CLUSTERING_HPP

#include <vector>
#include <utility>
#include <string>
#include <unordered_set>
#include <unordered_map>

using std::pair;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

template <class T>
class clustering
{
    int kClusters;
    int nHT;
    int nHF;

    // T ** centroids; //An index of centroids
    // vector<T*> *clusters;
    unordered_map<T *, unordered_set<T *> > *centroids;
    vector<T *> *elements;
    T **centroidIndex;

public:
    clustering(int, int, int, vector<T *> *);
    void initialize();
    void initializePlus();
    T *binarySearch(vector<pair<T *, double> >, double);
    void rangeSearchAssign();
    void lloydAssign();
    int pamUpdate();
    int meanVectorUpdate();
    int dbaUpdate();
    void printResults(string, string, string, string, double, char);
    vector<double> *silhouetteEvaluate();
    ~clustering();
};

#endif