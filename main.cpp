#include <iostream>
#include "Dataset.hpp"
#include "PositionScoringMatrix.hpp"
#include "ConfusionMatrix.hpp"

void evalPositionScoringMatrix(Dataset test_dataset, PositionScoringMatrix myMatrix, float thresholde){
    ConfusionMatrix myCM;

    for(int i=0; i<test_dataset.GetNbrSequences(); i++){
        int cleavage = test_dataset.GetCleavage(i);
        std::string sequence = test_dataset.GetSequence(i);
        for(int j=1 + myMatrix.Getp(); j + myMatrix.Getq() - 1 < sequence.size(); j++){
            std::string test_data = sequence.substr(j - myMatrix.Getp(), myMatrix.Getp() + myMatrix.Getq());
            myCM.AddPrediction((j == cleavage), myMatrix.prediction(test_data, thresholde));
        }
    }

    myCM.PrintEvaluation();
}

int main(int argc, char* argv[]) {
    Dataset GS_13("./GRAM+SIG_13.red.txt");
    std::cout << GS_13.GetSequence(0) << std::endl;
    std::cout << GS_13.GetCleavage(0) << std::endl;
    Dataset GS_13_train(GS_13, true, 0.8);
    Dataset GS_13_test(GS_13, false, 0.8);
    std::cout << GS_13_train.GetSequence(0) << std::endl;
    std::cout << GS_13_train.GetCleavage(0) << std::endl;
    std::cout << GS_13_test.GetSequence(0) << std::endl;
    std::cout << GS_13_test.GetCleavage(0) << std::endl;

    PositionScoringMatrix myMatrix(GS_13_train, 13, 2, 0.001);
    float min_score;
    float avg_score;
    float max_score;
    myMatrix.getTrianScore(GS_13_train, min_score, avg_score, max_score);
    std::cout << min_score << avg_score << max_score << std::endl;

    std::string test_data = GS_13_test.GetSequence(0).substr(GS_13_test.GetCleavage(0) - 13, 15);
    std::cout << myMatrix.prediction(test_data, avg_score) << std::endl;

    test_data = GS_13_test.GetSequence(0).substr(GS_13_test.GetCleavage(0) - 14, 15);
    std::cout << myMatrix.prediction(test_data, avg_score) << std::endl;

    evalPositionScoringMatrix(GS_13_test, myMatrix, avg_score);
}