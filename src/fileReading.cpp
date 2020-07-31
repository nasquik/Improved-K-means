#include <iostream>
#include <fstream>
#include <sstream>

#include "dataStructs.hpp"
#include "fileReading.hpp"

using namespace std;

vector<class Point*>* Reading::readPoints(string FileLocation)
{
    string line;

    //count the lines of the file

    int count = 0;
    int columns = 0;
    string x;
    ifstream file(FileLocation);

    getline(file, line);


    if (getline(file, line))
    {
        count++;
        //count the columns
        istringstream buffer(line);
        while (buffer >> x)
        {
            columns++;
        }
        columns--; //exclude the Point name
    }

    while (getline(file, line))
        count++;

    cout << "Number of lines in file: " << count << endl;
    cout << "Number of dimensions per point: " << columns << endl;

    //construct the table of Points

    ifstream infile(FileLocation);

    getline(infile, line);

    vector<class Point*>* table = new vector<class Point *>;
    table->reserve(count);

    string ID, buff;
    double *coord;
    
    while (getline(infile, line))
    {
        istringstream iss(line);

        iss >> ID;

        coord = new double[columns];

        for (int i = 0; i < columns; i++)
        {
            iss >> coord[i];
        }

        class Point *newPoint = new class Point(ID, coord, columns);
        table->push_back(newPoint);
    }

    return table;
}

vector<class Curve*>* Reading::readCurves(string FileLocation)
{
    string line;

    //count the lines of the file

    int count = 0;
    int columns = 0;
    string x;
    ifstream file(FileLocation);

    getline(file, line);

    if (getline(file, line))
    {
        count++;
        //count the columns
        istringstream buffer(line);
        while (buffer >> x)
        {
            columns++;
        }
        columns--; //exclude the Curve name
    }

    while (getline(file, line))
        count++;

    cout << "Number of lines in file: " << count << endl;
    cout << "Number of dimensions per point: " << columns << endl;

    vector<class Curve*>* table = new vector<class Curve *>;
    table->reserve(count);

    string ID, buff;
    int pointsNumber;
    


    ifstream infile(FileLocation);
    
    getline(infile, line);

    while (getline(infile, line))
    {
        istringstream iss(line);
        iss >> ID;
        iss >> pointsNumber;
        pair<double, double> *coord = new pair<double, double>[pointsNumber];
        pair<double, double> mypair;
        
        for (size_t i = 0; i < pointsNumber; i++)
        {
            //read first part of the string
            iss >> buff;
            sscanf(buff.c_str(), "(%lf,", &coord[i].first);
            iss >> buff;
            sscanf(buff.c_str(), "%lf)", &coord[i].second);
        }

        class Curve *curve = new class Curve(ID, coord, pointsNumber, 0);
        table->push_back(curve);
    }
    return table;
}

int readConfig(string confFile, int *kClusters, int *nGrids, int *nHT, int *nHashFunctions)
{

    string cline;
    ifstream cfile(confFile);

    if(cfile.fail())
        return 1;

    while (getline(cfile, cline))
    {
        istringstream buffer(cline);
        string command;
        buffer >> command;

        if (command == "number_of_clusters:")
        {
            buffer >> *kClusters;
        }
        else if (command == "number_of_grids:")
        {
            buffer >> *nGrids;
        }
        else if (command == "number_of_vector_hash_tables:")
        {
            buffer >> *nHT;
        }
        else if (command == "number_of_vector_hash_functions:")
        {
            buffer >> *nHashFunctions;
        }
        else{
            cout << "Unexpected Input!" << endl;
            return -1;
        }
    }

    cout << "Cluster: " << *kClusters << endl;
    cout << "Grids: " << *nGrids << endl;
    cout << "nHT: " << *nHT << endl;
    cout << "HashFunctions: " << *nHashFunctions << endl;

    return 0;
}