#include <iostream>
#include <fstream>
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

void toSVM(Dataset dataset, int p, int q, const char* output_file){

    std::ofstream outfile;
    outfile.open(output_file, std::ios::out | std::ios::trunc );

    for(int i=0; i<dataset.GetNbrSequences(); i++){
        int cleavage = dataset.GetCleavage(i);
        std::string sequence = dataset.GetSequence(i);
        for(int j=1 + p; j + q - 1 < sequence.size(); j++){
            std::string line;
            if(j == cleavage) line = "1 ";
            else line = "-1 ";
            for(int k=0; k< p+q; k++){
                int a = (sequence[j-p + k] - 'A') + k*26 + 1;
                line += std::to_string(a);
                line += ":1.00 ";
            }
            outfile << line << std::endl;
        }
    }

    outfile.close();
}

int main(int argc, char* argv[]) {
    Dataset GS_13("./data/GRAM+SIG_13.red.txt");
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

    toSVM(GS_13_train, 13, 2, "./data/GRAM+SIG_13_train.svm");
    toSVM(GS_13_test, 13, 2, "./data/GRAM+SIG_13_test.svm");
}