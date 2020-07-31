CXX = g++
FLAGS = -g -c #-std=c++14
#CFLAGS=-I
DEPS = fileReading.hpp dataStructs.hpp clustering.hpp metrics.hpp LSH.hpp dynamicTimeWarping.hpp hashTable.hpp objectiveFunctions.hpp bruteForce.hpp gridCurve.hpp
OBJ = fileReading.o dataStructs.o clustering.o metrics.o LSH.o dynamicTimeWarping.o hashTable.o objectiveFunctions.o bruteForce.o gridCurve.o
OUT = cluster mainTest

all: $(OUT)

cluster: main.o $(OBJ)
	$(CXX) -g main.o $(OBJ) -o cluster

mainTest: mainTest.cpp mainTest.o fileReading.o dataStructs.o dynamicTimeWarping.o metrics.o
	$(CXX) -g mainTest.o fileReading.o dataStructs.o dynamicTimeWarping.o metrics.o -o mainTest -L /usr/lib -l gtest -l pthread

main.o: main.cpp $(OBJ) $(DEPS)
	$(CXX) $(FLAGS) main.cpp -o main.o

mainTest.o: mainTest.cpp fileReading.cpp fileReading.hpp dataStructs.hpp
	$(CXX) $(FLAGS) -Wall -I /usr/local/include mainTest.cpp -o mainTest.o

clustering.o: clustering.cpp clustering.hpp
	$(CXX) $(FLAGS) clustering.cpp -o clustering.o

fileReading.o: fileReading.cpp fileReading.hpp
	$(CXX) $(FLAGS) fileReading.cpp -o fileReading.o

dataStructs.o: dataStructs.cpp dataStructs.hpp
	$(CXX) $(FLAGS) dataStructs.cpp -o dataStructs.o

metrics.o: metrics.cpp metrics.hpp
	$(CXX) $(FLAGS) metrics.cpp -o metrics.o

LSH.o: LSH.cpp LSH.hpp
	$(CXX) $(FLAGS) LSH.cpp -o LSH.o

hashTable.o: hashTable.cpp hashTable.hpp
	$(CXX) $(FLAGS) hashTable.cpp -o hashTable.o

dynamicTimeWarping.o: dynamicTimeWarping.cpp dynamicTimeWarping.hpp
	$(CXX) $(FLAGS) dynamicTimeWarping.cpp -o dynamicTimeWarping.o

objectiveFunctions.o: objectiveFunctions.cpp objectiveFunctions.hpp
	$(CXX) $(FLAGS) objectiveFunctions.cpp -o objectiveFunctions.o

bruteForce.o: bruteForce.cpp bruteForce.hpp
	$(CXX) $(FLAGS) bruteForce.cpp -o bruteForce.o

gridCurve.o: gridCurve.cpp gridCurve.hpp
	$(CXX) $(FLAGS) gridCurve.cpp -o gridCurve.o

clean:
	rm $(OUT) *.o
