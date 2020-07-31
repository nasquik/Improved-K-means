#include <iostream>
#include <gtest/gtest.h>
#include "fileReading.hpp"
#include <string>
#include <utility>

#include "fileReading.hpp"
#include "dataStructs.hpp"
#include "dynamicTimeWarping.hpp"

using std::string;
using std::vector;
using std::pair;
using std::make_pair;

TEST(ReadInput, readCorrectPoints){
    class Reading reader;
    vector<class Point *> *inputTable;
    inputTable = reader.readPoints("./Ex2_Datasets/DataVectors_5_500x100.csv");
    int x = inputTable->size();
    EXPECT_EQ(456,x);
    //free memory
    while (!inputTable->empty())
        {
            delete inputTable->back();
            inputTable->pop_back();
        }
    delete (inputTable);

}

TEST(ReadInput, readInvalidFile){
    class Reading reader;
    vector<class Point *> *inputTable;
    inputTable = reader.readPoints("./Ex2_Datasets/DataVectors");
    int x = inputTable->size();
    EXPECT_EQ(0,x);
    //free memory
    while (!inputTable->empty())
        {
            delete inputTable->back();
            inputTable->pop_back();
        }
    delete (inputTable);

}

TEST(metrics, dtwDistanceZeroPoint){
    class Curve* curve1, *curve2;
    pair<double, double>  *index2;
    curve1 = new class Curve("emptyCurve",NULL,0,0);
    index2 = new pair<double, double>[2];
    index2[0] = make_pair(1,3);
    index2[1] = make_pair(1,4);
    curve2 = new class Curve("nonEmpty",index2,2,0);
    EXPECT_DOUBLE_EQ(-1.0, dtwDist(curve1,curve2));

    delete curve1;
    delete curve2;


}

TEST(metrics, dtwDistance){
    class Curve* curve1, *curve2;
    pair<double, double> *index1, *index2;
    index1 = new pair<double, double>[4];
    index1[0] = make_pair(0.3,3.5);
    index1[1] = make_pair(3,4);
    index1[2] = make_pair(4,9);
    index1[3] = make_pair(2,2);
    curve1 = new class Curve("curve1",index1,4,0);
    index2 = new pair<double, double>[2];
    index2[0] = make_pair(1,3);
    index2[1] = make_pair(1,4);

    curve2 = new class Curve("curve2",index2,2,0);
    EXPECT_DOUBLE_EQ(10.927252399049353, dtwDist(curve1,curve2));

    delete curve1;
    delete curve2;

}

TEST(readConfFileTest, ExpectedInput)
{ 
    int kClusters, nGrids, nHT, nHashFunctions;
    ASSERT_EQ(0, readConfig("cluster.conf", &kClusters, &nGrids, &nHT, &nHashFunctions));
}
 
TEST(readConfFileTest, FileNotExist)
{
    int kClusters, nGrids, nHT, nHashFunctions;
    ASSERT_EQ(1, readConfig("cluster.con", &kClusters, &nGrids, &nHT, &nHashFunctions));
}

TEST(readConfFileTest, UnexpectedInput)
{ 
    int kClusters, nGrids, nHT, nHashFunctions;
    ASSERT_EQ(-1, readConfig("fakeCluster.conf", &kClusters, &nGrids, &nHT, &nHashFunctions));
}

TEST(metrics, dtwDistanceEqual){
    class Curve* curve1;
    pair<double, double> *index1;
    index1 = new pair<double, double>[4];
    index1[0] = make_pair(0.3,3.5);
    index1[1] = make_pair(3,4);
    index1[2] = make_pair(4,9);
    index1[3] = make_pair(2,2);
    curve1 = new class Curve("identical",index1,4,0);

    EXPECT_DOUBLE_EQ(0, dtwDist(curve1,curve1));

    delete curve1;
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv); 
    return RUN_ALL_TESTS();
}