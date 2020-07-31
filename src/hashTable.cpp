#include <iostream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "dataStructs.hpp"
#include "hashTable.hpp"

using namespace std;

// template <>
// void HashTable<map<unsigned int, char> >::initialize(size_t size, int givenw, int givenNumFunct, int givendim, int givenMaxPoints)
// {
//     numFunct = givenNumFunct;
//     w = givenw;
//     dimensions = givendim;
//     maxPoints = givenMaxPoints;

//     // create buckets
//     buckets = new map<unsigned int, char>;

//     // initialize M

//     unsigned int temp = 32 / numFunct;
//     M = 1;
//     M = M << temp;

//     // create si arrays (fill with uniformly random numbers from 0 to w)

//     sVectors = new double *[numFunct];
//     for (int i = 0; i < numFunct; i++)
//     {
//         sVectors[i] = new double[dimensions];
//     }

//     uniform_real_distribution<> dis(0.0, (double)w);

//     for (int i = 0; i < numFunct; i++)
//     {
//         for (int j = 0; j < dimensions; j++)
//         {
//             sVectors[i][j] = dis(generator);
//         }
//     }

//     // create array to keep powers of m

//     mArray = new unsigned int[dimensions];
//     unsigned int m = INT32_MAX - 4;
//     for (int i = 0; i < dimensions; i++)
//     {
//         mArray[i] = modular_expo(m, i, M);
//     }
// }

template <class T>
void HashTable<T>::initialize(size_t size, int givenw, int givenNumFunct, int givendim, int givenMaxPoints)
{
    numFunct = givenNumFunct;
    w = givenw;
    dimensions = givendim;
    maxPoints = givenMaxPoints;

    // create buckets

    bucketSize = size;
    buckets = new unordered_map<unsigned int, vector<class Point*> >();

    // initialize M

    unsigned int temp = 32 / numFunct;
    M = 1;
    M = M << temp;

    // create si arrays (fill with uniformly random numbers from 0 to w)

    sVectors = new double *[numFunct];
    for (int i = 0; i < numFunct; i++)
    {
        sVectors[i] = new double[dimensions];
    }

    uniform_real_distribution<> dis(0.0, (double)w);

    for (int i = 0; i < numFunct; i++)
    {
        for (int j = 0; j < dimensions; j++)
        {
            sVectors[i][j] = dis(generator);
        }
    }

    // create array to keep powers of m

    mArray = new unsigned int[dimensions];
    unsigned int m = INT32_MAX - 4;
    for (int i = 0; i < dimensions; i++)
    {
        mArray[i] = modular_expo(m, i, M);
    }
}

template <>
void HashTable<vector<class Point *> >::printHashTable()
{

    // for (int i = 0; i < bucketSize; i++)
    // {
    //     cout << buckets[i].size() << " ";
    // }
    // cout << endl;
    // for(auto& it: *buckets){
    //     cout << it.second.size() << " ";
    // }
    // cout <<endl;
}

template <class T>
void HashTable<T>::printHashTable()
{
}

template <class T>
HashTable<T>::HashTable() : generator((std::random_device())())
{
}

template <>
HashTable<map<unsigned int, char> >::~HashTable()
{
    //clear the structures and the allocated memory
    buckets->clear();
    delete buckets;
    for (int i = 0; i < numFunct; i++)
    {
        delete[] sVectors[i];
    }
    delete[] sVectors;
    delete[] mArray;
}

template <class T>
HashTable<T>::~HashTable()
{
    delete buckets;
    for (int i = 0; i < numFunct; i++)
    {
        delete[] sVectors[i];
    }
    delete[] sVectors;
    delete[] mArray;
}

template <>
int HashTable<vector<class Point*> >::insertPoint(class Point *point)
{
    unsigned int result = amplifiedHashFunctionPoint(point);
    //insert point to the dataset to the position indicated by the hashFunction
    unsigned int position = result;
    if (buckets->find(position)==buckets->end())
        buckets->insert(make_pair(position,vector<class Point*>()));
    
    buckets->at(position).push_back(point);

    return 0;
}

// template <>
// int HashTable<map<unsigned int, char> >::insertPoint(class Point *point)
// {
//     unsigned int result = amplifiedHashFunctionPoint(point);
//     char random;
//     map<unsigned int, char>::iterator it;

//     it = buckets->find(result);

//     if (it == buckets->end())
//     {
//         uniform_int_distribution<> dis(0, 1);
//         random = dis(generator);
//         buckets->insert(pair<unsigned int, unsigned int>(result, random));
//         return random;
//     }
//     else
//     {
//         return it->second;
//     }
// }

template <class T>
int HashTable<T>::insertPoint(class Point *point)
{
    return 0;
};

template <class T>
unsigned int HashTable<T>::amplifiedHashFunctionPoint(class Point *x)
{

    int shiftAmount = 32 / numFunct;
    unsigned int temp, result = 0;
    for (int i = 0; i < numFunct; i++)
    {
        temp = hashFunctionPoint(x, i);
        temp = temp << shiftAmount * i; //prepare the 32/numFunct binary digits for concatenation
        result += temp;
    }
    return result;
}

template <class T>
unsigned int HashTable<T>::hashFunctionPoint(class Point *x, int functionNo)
{

    int *a = new int[dimensions];
    double *s = sVectors[functionNo];
    unsigned int m = UINT32_MAX - 4;



    for (int i = 0; i < dimensions; i++)
    {
        a[i] = (int)floor((x->getCoord()[i] - s[i]) / w);
    }



    unsigned int result = 0;
    unsigned int step1;

    for (int i = 0; i < dimensions; i++)
    {
        step1 = (a[dimensions - 1 - i] % M + M) % M;
        result += (step1 * mArray[i]) % M;
    }
    result = result % M;
    delete[] a;

    return result;
}

unsigned int modular_expo(unsigned int base, unsigned int exponent, unsigned int modulus)
{
    //function that calculates modulos of exponentials that produce very big values
    if (modulus == 1)
        return 0;
    if (exponent == 0)
        if (modulus == 1)
            return 1;
        else
            return 0;

    unsigned long long c = 1;
    for (unsigned int i = 1; i <= exponent; i++)
    {
        c = (c * base) % modulus;
    }
    return c;
}

template class HashTable<vector<class Point *> >;
template class HashTable<class Curve >;