#ifndef LSH_HPP
#define LSH_HPP

#include <vector>
using std::vector;
#include <unordered_map>
using std::unordered_map;
#include <unordered_set>
using std::unordered_set;

#include "hashTable.hpp"

template <class T>
class LSH
{
    int k;
    int w;
    int L;
    vector<class Point *> *input;
    class HashTable<vector<class Point *> > *hashTables;

public:
    LSH(int, int, int, vector<class Point *> *);
    ~LSH();
    int getk();
    int getw();
    int getL();
    T findNN(class Point *, double *);
    vector <class Point *> &returnNeighb(class Point*,int *);
    vector<pair<class Point *, double> > *findRadiusNN(class Point *, double);
    void assign(unordered_map<T,unordered_set<T> >*,vector<T>*);
};

#endif