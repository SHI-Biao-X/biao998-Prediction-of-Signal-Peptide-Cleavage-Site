#include <iostream>
#include <fstream>
#include <array>
#include <cstdio>
#include <memory>
#include <cmath>
#include "Dataset.hpp"
#include "PositionScoringMatrix.hpp"
#include "ConfusionMatrix.hpp"

double evalPositionScoringMatrix(Dataset test_dataset, PositionScoringMatrix myMatrix, float thresholde){
    ConfusionMatrix myCM;

    for(int i=0; i<test_dataset.GetNbrSequences(); i++){
        int cleavage = test_dataset.GetCleavage(i);
        std::string sequence = test_dataset.GetSequence(i);
        for(int j=1 + myMatrix.Getp(); j + myMatrix.Getq() - 1 < sequence.size(); j++){
            std::string test_data = sequence.substr(j - myMatrix.Getp(), myMatrix.Getp() + myMatrix.Getq());
            myCM.AddPrediction((j == cleavage), myMatrix.prediction(test_data, thresholde));
        }
    }

    // myCM.PrintEvaluation();
    return myCM.f_score();
}

void encodingToSVM(Dataset dataset, int p, int q, const char* output_file){

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
                line += ":1 ";
            }
            outfile << line << std::endl;
        }
    }

    outfile.close();
}

std::string exec(const char* in) {
    std::array<char, 128> buffer;
    std::string out;
    std::unique_ptr<FILE, decltype(&_pclose)> pipe(_popen(in, "r"), _pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        out += buffer.data();
    }
    return out;
}

void findPAndQ(Dataset train_dataset, Dataset test_dataset, std::vector<int> L_p, std::vector<int> L_q) {
    int n_p = L_p.size();
    int n_q = L_q.size();

    float alpha = 0.001;
    
    double max_f_score = 0;
    int max_p = -1;
    int max_q = -1;

    for (int i = 0; i < n_p; i++)
    {
        for (int j = 0; j < n_q; j++)
        {
            PositionScoringMatrix myMatrix(&train_dataset, L_p[i], L_q[j], alpha);
            
            float min_score;
            float avg_score;
            float max_score;
            myMatrix.getTrianScore(min_score, avg_score, max_score);

            double f_score;
            f_score = evalPositionScoringMatrix(test_dataset, myMatrix, avg_score);

            if(f_score > max_f_score){
                max_f_score = f_score;
                max_p = L_p[i];
                max_q = L_q[j];
            }

            std::cout << "p = " << L_p[i] << ", q = " << L_q[j] << ", Fscore = " << f_score << std::endl;
        }
    }

    std::cout << "max p = " << max_p << ", max q = " << max_q << ", max Fscore = " << max_f_score << std::endl;
}

void findCandGamma(const char* train_file, int kernel_type, std::vector<double> L_c, std::vector<double> L_g = {}){
    double max_f_score = 0;
    double max_c = 1;
    double max_g = 0;

    if(kernel_type == 0 || kernel_type == 4 || kernel_type == 6)
    {
        if(!L_g.empty()) {
            std::cout << "Error: L_g is not empty" << std::endl;
            return;
        }
        for(int i = 0 ; i < L_c.size(); i++){
            std::string cmdline = "";
            cmdline += ".\\libsvm-3.23\\svm-train -t ";
            cmdline += std::to_string(kernel_type);
            cmdline += " -c ";
            cmdline += std::to_string(L_c[i]);
            cmdline += " -v 5 ";
            cmdline += train_file;
            
            std::string output = exec(cmdline.c_str());
            std::size_t found1 = output.find("Cross");
            found1 = output.find("=", found1);
            found1 += 2;
            std::size_t found2 = output.find("%", found1);

            float f_score = std::stod(output.substr(found1, found2 - found1))/100;
            if(f_score > max_f_score){
                max_f_score = f_score;
                max_c = L_c[i];
            }

            std::cout << "c = " << L_c[i] << ", Fscore = " << f_score << std::endl;
        }

        std::cout << "max c = " << max_c << ", max Fscore = " << max_f_score << std::endl;
        return;
    }

    else if(kernel_type == 2 || kernel_type == 7)
    {
        if(L_g.empty()) {
            std::cout << "Error: L_g is empty" << std::endl;
            return;
        }
        for(int i = 0 ; i < L_c.size(); i++){
            for(int j = 0 ; j < L_g.size(); j++){
                std::string cmdline = "";
                cmdline += ".\\libsvm-3.23\\svm-train -t ";
                cmdline += std::to_string(kernel_type);
                cmdline += " -c ";
                cmdline += std::to_string(L_c[i]);
                cmdline += " -g ";
                cmdline += std::to_string(L_g[j]);
                cmdline += " -v 5 ";
                cmdline += train_file;
                
                // std::cout << cmdline << std::endl;
                std::string output = exec(cmdline.c_str());
                std::size_t found1 = output.find("Cross");
                found1 = output.find("=", found1);
                found1 += 2;
                std::size_t found2 = output.find("%", found1);

                float f_score = std::stod(output.substr(found1, found2 - found1))/100;
                if(f_score > max_f_score){
                    max_f_score = f_score;
                    max_c = L_c[i];
                    max_g = L_g[j];
                }
                
                std::cout << "c = " << L_c[i] << ", gamma = " << L_g[j] << ", Fscore = " << f_score << std::endl;
            }
        }

        std::cout << "max c = " << max_c << ", max g = " << max_g <<", max Fscore = " << max_f_score << std::endl;
        return;
    }
    
    else
    {
        std::cout << "Error: kernel type is not valid" << std::endl;
        return;
    }
}

int main(int argc, char* argv[]) {
    
    Dataset GS_13("./data/GRAM+SIG_13.red");
    // std::cout << GS_13.GetSequence(0) << std::endl;
    // std::cout << GS_13.GetCleavage(0) << std::endl;
    Dataset GS_13_train(GS_13, true, 0.8);
    Dataset GS_13_test(GS_13, false, 0.8);

    //Fonction test
    /*
    std::cout << GS_13_train.GetSequence(0) << std::endl;
    std::cout << GS_13_train.GetCleavage(0) << std::endl;
    std::cout << GS_13_test.GetSequence(0) << std::endl;
    std::cout << GS_13_test.GetCleavage(0) << std::endl;

    PositionScoringMatrix myMatrix(&GS_13_train, 13, 2, 0.001);
    float min_score;
    float avg_score;
    float max_score;
    myMatrix.getTrianScore(min_score, avg_score, max_score);
    std::cout << min_score << avg_score << max_score << std::endl;

    std::string test_data = GS_13_test.GetSequence(0).substr(GS_13_test.GetCleavage(0) - 13, 15);
    std::cout << myMatrix.prediction(test_data, avg_score) << std::endl;

    test_data = GS_13_test.GetSequence(0).substr(GS_13_test.GetCleavage(0) - 14, 15);
    std::cout << myMatrix.prediction(test_data, avg_score) << std::endl;

    std::cout << evalPositionScoringMatrix(GS_13_test, myMatrix, avg_score) << std::endl;

    encodingToSVM(GS_13_train, 13, 2, "./data/GRAM+SIG_13_train.svm");
    encodingToSVM(GS_13_test, 13, 2, "./data/GRAM+SIG_13_test.svm");

    myMatrix.matrixToSVM("./data/GRAM+SIG_13_kernels_train.svm");
    myMatrix.matrixTestToSVM(GS_13_test, "./data/GRAM+SIG_13_kernels_test.svm");
    */

    // find best p and q
/*
    std::vector<int> L_p({8,9,10,11,12,13});
    std::vector<int> L_q({1,2,3,4,5,6,7});

    findPAndQ(GS_13_train, GS_13_test, L_p, L_q);
    // max p = 8, max q = 2, max Fscore = 0.36
*/
/*
    std::vector<int> L_p({4,5,6,7,8,9});
    std::vector<int> L_q({1,2,3,4});

    findPAndQ(GS_13_train, GS_13_test, L_p, L_q);
    // max p = 6, max q = 2, max Fscore = 0.421053
*/

    // prepare data for p = 8, q = 2
/*
    encodingToSVM(GS_13_train, 8, 2, "./data/p8q2/GRAM+SIG_13_train.svm");
    encodingToSVM(GS_13_test, 8, 2, "./data/p8q2/GRAM+SIG_13_test.svm");

    PositionScoringMatrix myMatrix(&GS_13_train, 8, 2, 0.001);

    myMatrix.matrixToSVM("./data/p8q2/GRAM+SIG_13_kernels_train.svm");
    myMatrix.matrixTestToSVM(GS_13_test, "./data/p8q2/GRAM+SIG_13_kernels_test.svm");
*/

    // prepare data for p = 6, q = 2
/*
    encodingToSVM(GS_13_train, 6, 2, "./data/p6q2/GRAM+SIG_13_train.svm");
    encodingToSVM(GS_13_test, 6, 2, "./data/p6q2/GRAM+SIG_13_test.svm");

    PositionScoringMatrix myMatrix(&GS_13_train, 6, 2, 0.001);

    myMatrix.matrixToSVM("./data/p6q2/GRAM+SIG_13_kernels_train.svm");
    myMatrix.matrixTestToSVM(GS_13_test, "./data/p6q2/GRAM+SIG_13_kernels_test.svm");
*/


    // p = 8, q = 2
    //rbf
/*
    std::vector<double> L_c({pow(2,-1),pow(2,0),pow(2,1),pow(2,2)});
    std::vector<double> L_g({pow(2,-13),pow(2,-9),pow(2,-5),pow(2,-1)});
    findCandGamma("./data/p8q2/GRAM+SIG_13_train.svm", 2, L_c, L_g);
    // max c = 4, max g = 0.03125, max Fscore = 0.312057
*/
/*
    std::vector<double> L_c({pow(2,1),pow(2,2),pow(2,3),pow(2,4)});
    std::vector<double> L_g({pow(2,-6),pow(2,-5),pow(2,-4),pow(2,-3)});
    findCandGamma("./data/p8q2/GRAM+SIG_13_train.svm", 2, L_c, L_g);
    // max c = 8, max g = 0.03125, max Fscore = 0.434783
*/

    // proba kernel
/*
    std::vector<double> L_c({pow(2,-1),pow(2,0),pow(2,1),pow(2,2)});
    findCandGamma("./data/p8q2/GRAM+SIG_13_kernels_train.svm", 4, L_c);
    // max c = 0.5, max Fscore = 0.383178
*/
/*
    std::vector<double> L_c({pow(2,-4),pow(2,-3),pow(2,-2),pow(2,-1)});
    findCandGamma("./data/p8q2/GRAM+SIG_13_kernels_train.svm", 4, L_c);
    // max c = 0.125, max Fscore = 0.404494
*/


    // p = 6, q = 2
    //rbf
/*
    std::vector<double> L_c({pow(2,-1),pow(2,0),pow(2,1),pow(2,2)});
    std::vector<double> L_g({pow(2,-13),pow(2,-9),pow(2,-5),pow(2,-1)});
    findCandGamma("./data/p6q2/GRAM+SIG_13_train.svm", 2, L_c, L_g);
    // max c = 4, max g = 0.03125, max Fscore = 0.186047
*/
/*
    std::vector<double> L_c({pow(2,1),pow(2,2),pow(2,3),pow(2,4)});
    std::vector<double> L_g({pow(2,-6),pow(2,-5),pow(2,-4),pow(2,-3)});
    findCandGamma("./data/p6q2/GRAM+SIG_13_train.svm", 2, L_c, L_g);
    // max c = 8, max g = 0.0625, max Fscore = 0.453488
*/

    //similarity score matrix

    std::vector<double> L_c({pow(2,1),pow(2,2),pow(2,3),pow(2,4)});
    findCandGamma("./data/p6q2/GRAM+SIG_13_train.svm", 1, L_c);
    // max c = 8, max g = 0.0625, max Fscore = 0.453488

    // proba kernel
/*
    std::vector<double> L_c({pow(2,-1),pow(2,0),pow(2,1),pow(2,2)});
    findCandGamma("./data/p6q2/GRAM+SIG_13_kernels_train.svm", 4, L_c);
    // max c = 0.5, max Fscore = 0.427861
*/

/*    std::vector<double> L_c({pow(2,-4),pow(2,-3),pow(2,-2),pow(2,-1)});
    findCandGamma("./data/p6q2/GRAM+SIG_13_kernels_train.svm", 4, L_c);
    // max c = 0.5, max Fscore = 0.427861
*/



}