CXX = g++
CXXFLAGS = --std=c++11

.PHONY: all clean

all: grader

# Specific part

SOURCES_SPECIFIC = Dataset.cpp PositionScoringMatrix.cpp ConfusionMatrix.cpp
OBJECTS_SPECIFIC = Dataset.o PositionScoringMatrix.o ConfusionMatrix.o

Dataset.o: Dataset.cpp Dataset.hpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<
PositionScoringMatrix.o: PositionScoringMatrix.cpp PositionScoringMatrix.hpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<
ConfusionMatrix.o: ConfusionMatrix.cpp ConfusionMatrix.hpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<

# Common part

SOURCES_COMMON = main.cpp
OBJECTS_COMMON = main.o

grader: $(OBJECTS_COMMON) $(OBJECTS_SPECIFIC)
	$(CXX) $(CXXFLAGS) -o grader $(OBJECTS_COMMON) $(OBJECTS_SPECIFIC)

main.o: main.cpp 
	$(CXX) -c $(CXXFLAGS) -o main.o main.cpp

clean:
	del /F grader.exe *.o