#include <iostream>
#include <string>
#include <vector>

#include "dataStructs.hpp"
#include "fileReading.hpp"
#include "metrics.hpp"
#include "dynamicTimeWarping.hpp"

using namespace std;

double dtwDist(class Curve *P, class Curve *Q)
{   
    int length = P->getSize();
    int width = Q->getSize();

    if(length == 0 || width == 0){
        return -1;
    }

    double **C = new double *[length];

    for (int i = 0; i < length; i++)
        C[i] = new double[width];

    double minimum, result = 0.0;

    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (i == 0 && j == 0)
                C[i][j] = euclideanDist(P->getCoord()[i], Q->getCoord()[j]);
            else if (j > 0 && i == 0)
                C[i][j] = C[i][j - 1] + euclideanDist(P->getCoord()[i], Q->getCoord()[j]);
            else if (i > 0 && j == 0)
                C[i][j] = C[i - 1][j] + euclideanDist(P->getCoord()[i], Q->getCoord()[j]);
            else
            {
                minimum = min(min(C[i - 1][j], C[i - 1][j - 1]), C[i][j - 1]);
                C[i][j] = minimum + euclideanDist(P->getCoord()[i], Q->getCoord()[j]);
                minimum = 0.0;
            }
        }
    }

    result = C[length - 1][width - 1];




    for (int i = 0; i < length; i++)
        delete[] C[i];

    delete[] C;

    return result;
}

vector<pair<int, int> >* dtwBestTraversal(class Curve* P, class Curve* Q)
{   
    int length = P->getSize();
    int width = Q->getSize();

    double **C = new double *[length];

    for (int i = 0; i < length; i++)
        C[i] = new double[width];

    double minimum;
    vector<pair<int, int> >* result;

    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < width; j++)
        {
            if (i == 0 && j == 0)
                C[i][j] = euclideanDist(P->getCoord()[i], Q->getCoord()[j]);
            else if (j > 0 && i == 0)
                C[i][j] = C[i][j - 1] + euclideanDist(P->getCoord()[i], Q->getCoord()[j]);
            else if (i > 0 && j == 0)
                C[i][j] = C[i - 1][j] + euclideanDist(P->getCoord()[i], Q->getCoord()[j]);
            else
            {
                minimum = min(min(C[i - 1][j], C[i - 1][j - 1]), C[i][j - 1]);
                C[i][j] = minimum + euclideanDist(P->getCoord()[i], Q->getCoord()[j]);
                minimum = 0.0;
            }
        }
    }

    result = optimalPath(C, length, width);

    for (int i = 0; i < length; i++)
        delete[] C[i];

    delete[] C;

    return result;
}

vector <pair<int, int> > * optimalPath(double** C, int length, int width){

    vector <pair<int, int> > *path = new vector <pair<int, int> >();
    path->insert(path->begin(),make_pair(length-1, width-1));

    int i = length - 1;
    int j = width - 1;

    while( i != 0 || j != 0){
        if (i == 0 && j == 0){
            //ignore
        }
        else if (i == 0 && j > 0){
            path->insert(path->begin(),make_pair(i,j-1));
            j--;
        }
        else if (i>0 && j == 0){
            path->insert(path->begin(),make_pair(i-1,j));
            i--;
        }
        else{
            //compare the upper triangle
            int pos;

            if(C[i-1][j-1]< C[i][j-1]){
                if(C[i-1][j-1]<C[i-1][j]){
                    pos = 2;
                }
                else{
                    pos = 1;
                }
            }
            else{
                if(C[i][j-1]<C[i-1][j]){
                    pos = 3;
                }
                else{
                    pos = 1;
                }
            }

            if(pos == 1){
                path->insert(path->begin(),make_pair(i-1,j));
                i --;
            }
            else if(pos == 2){
                path->insert(path->begin(),make_pair(i-1,j-1));
                i--;
                j--;
            }
            else if(pos == 3){
                path->insert(path->begin(),make_pair(i,j-1));
                j--;
            }
        }
    }
    return path;
}