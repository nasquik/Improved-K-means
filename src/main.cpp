#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>

using namespace std;

#include "dataStructs.hpp"
#include "fileReading.hpp"
#include "clustering.hpp"

// $./cluster –i <input file> –c <configuration file> -o <output file> -
// complete <optional>

int setUp(string *, string *, string *);

int main(int argc, char *argv[])
{

    if (argc > 9)
    {
        cout << "Invalid Input!\n";
        return -1;
    }

    string inputFile, confFile, outputFile;
    int complete = 0;
    int grids = 2, L = 3;
    clock_t time_count;

    string optionBuffer, parameterBuffer;

    for (int i = 1; i <= (argc - 1); i += 2)
    {
        optionBuffer = string(argv[i]);
        if ((i + 1) < argc)
            parameterBuffer = string(argv[i + 1]);

        if (optionBuffer == "-i")
            inputFile = parameterBuffer;
        else if (optionBuffer == "-c")
            confFile = parameterBuffer;
        else if (optionBuffer == "-o")
            outputFile = parameterBuffer;
        else if (optionBuffer == "-complete")
            complete = 1;
    }

    int w = 4500;

    // read configuration file //

    int nGrids = 2, nHT = 3, kClusters, nHashFunctions;
    readConfig(confFile, &kClusters, &nGrids, &nHT, &nHashFunctions);

    // determine the methods that are to be used //

    string initialization, assignment, update;
    setUp(&initialization, &assignment, &update);

    // determine type of data that is to be handled //

    class Reading reader;

    string line, type;
    ifstream file(inputFile);

    getline(file, line);

    istringstream iss(line);
    iss >> type;

    // execute the algorithm //

    if (type == "vectors")
    {
        vector<class Point *> *inputTable;
        inputTable = reader.readPoints(inputFile);

        time_count = clock();
        class clustering<class Point> *clusteringHandler = new clustering<class Point>(kClusters, nHT, nHashFunctions, inputTable);
        if (initialization == "S")
            clusteringHandler->initialize();
        else
            clusteringHandler->initializePlus();

        for (int i = 0; i < 30; i++)
        {
            cout << "ITERATION: " << i+1 <<endl;
            if (i != 0)
            {
                if (update == "P"){
                   if (clusteringHandler->pamUpdate()){
                        break;
                    } 
                }   
                else
                {
                    if (clusteringHandler->meanVectorUpdate())
                        break;
                }
            }

            if (assignment == "L")
                clusteringHandler->lloydAssign();
            else
                clusteringHandler->rangeSearchAssign();
        }
        time_count = clock() - time_count;
        clusteringHandler->printResults(outputFile, initialization, assignment, update, (float)time_count / CLOCKS_PER_SEC , complete);
        delete clusteringHandler;

        while (!inputTable->empty())
        {
            delete inputTable->back();
            inputTable->pop_back();
        }
        delete (inputTable);
    }
    else if (type == "curves")
    {
        vector<class Curve *> *inputTable;
        inputTable = reader.readCurves(inputFile);

        time_count = clock();
        class clustering<class Curve> *clusteringHandler = new clustering<class Curve>(kClusters, nGrids, nHashFunctions, inputTable);

        if (initialization == "S")
            clusteringHandler->initialize();
        else
            clusteringHandler->initializePlus();

        for (int i = 0; i < 30; i++)
        {
            cout << "ITERATION: " << i+1 <<endl;
            if (i != 0)
            {
                if (update == "P"){
                    if (clusteringHandler->pamUpdate())
                        break;
                }
                else
                {
                    if (clusteringHandler->dbaUpdate())
                        break;
                }
            }

            if (assignment == "L")
                clusteringHandler->lloydAssign();
            else
                clusteringHandler->rangeSearchAssign();
        }
        time_count = clock() - time_count;

        clusteringHandler->printResults(outputFile, initialization, assignment, update, (float)time_count / CLOCKS_PER_SEC , complete);
        delete clusteringHandler;

        while (!inputTable->empty())
        {
            delete inputTable->back();
            inputTable->pop_back();
        }
        delete (inputTable);
    }
    else
    {
        cout << "Error reading type of input file!" << endl;
        return -1;
    }
}

int setUp(string *initialization, string *assignment, string *update)
{

    bool flag = false;

    while (flag == false)
    {
        cout << "Choose initialization method: (S)imple OR (K)-means++?" << endl;
        cin >> *initialization;
        if (*initialization == "S")
            flag = true;
        else if (*initialization == "K")
            flag = true;
        else
            cout << "Incorrect Input! Please try again!" << endl;
    }

    flag = false;

    while (flag == false)
    {
        cout << "Choose assignment method: (L)loyd's OR (R)ange-Search with LSH?" << endl;
        cin >> *assignment;
        if (*assignment == "L")
            flag = true;
        else if (*assignment == "R")
            flag = true;
        else
            cout << "Incorrect Input! Please try again!" << endl;
    }

    flag = false;

    while (flag == false)
    {
        cout << "Choose update method: (P)AM OR (M)eanVector?" << endl;
        cin >> *update;
        if (*update == "P")
            flag = true;
        else if (*update == "M")
            flag = true;
        else
            cout << "Incorrect Input! Please try again!" << endl;
    }

    return 0;
}