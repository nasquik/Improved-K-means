#ifndef _HASHTABLE_
#define _HASHTABLE_

#include <iostream>
#include <unordered_map>
#include "dataStructs.hpp"


#include <vector>
using std::vector;
#include <utility>
using std::pair;
#include <random>
using std::mt19937;
using std::random_device;
using std::unordered_map;
#include <cstring>

template <class T>
class HashTable
{
    unsigned int bucketSize;
    int numFunct;
    int w;
    int dimensions;
    int maxPoints;
    mt19937 generator;
    unsigned int M;
    unordered_map <unsigned int,vector<class Point*> > *buckets;
    double **sVectors;
    unsigned int *mArray;

public:
    HashTable();
    void initialize(size_t, int, int, int, int);
    void printHashTable();
    ~HashTable();

    //hashfunction operations
    unsigned int amplifiedHashFunctionPoint(class Point *x);
    unsigned int hashFunctionPoint(class Point *x, int functionNo);
    int insertPoint(class Point *);
    int findPoint(class Point *);
    vector<class Point*> &getneighbors(class Point *x, unsigned int *amplifiedResult,int* valid)
    {
        *amplifiedResult = amplifiedHashFunctionPoint(x);
        //std::cout << "Point: " << x->getID() << " Bucket: " << *amplifiedResult % bucketSize << std::endl;
        if (buckets->find(*amplifiedResult)==buckets->end()){
            *valid = 0;
            return buckets->begin()->second;
        }
        *valid = 1;
        return buckets->at(*amplifiedResult);//buckets[*amplifiedResult % bucketSize];
    }
    vector<class Point*> getDigit(class Point *x, unsigned int *amplifiedResult)
    {
        *amplifiedResult = amplifiedHashFunctionPoint(x);
        return buckets->at(*amplifiedResult);
    }
};

unsigned int modular_expo(unsigned int base, unsigned int exponent, unsigned int modulus);

#endif