#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include "PositionScoringMatrix.hpp"

void cutDataset(Dataset* dataset, int p, int q, int& n_intances, std::vector<std::vector<int>>* intances, std::vector<int>* labels){
    n_intances = 0;
    intances->clear();
    labels->clear();
    
    for(int i=0; i<dataset->GetNbrSequences(); i++){
        int cleavage = dataset->GetCleavage(i);
        std::string sequence = dataset->GetSequence(i);
        for(int j=1 + p; j + q - 1 < sequence.size(); j++){
            n_intances ++;

            if(j == cleavage) labels->push_back(1);
            else labels->push_back(-1);

            std::vector<int> intance(p+q);
            for(int k=0; k< p+q; k++){
                int a = sequence[j-p + k] - 'A';
                intance[k] = a;
            }
            intances->push_back(intance);
        }
    }
    // std::cout << n_intances << std::endl;
}

PositionScoringMatrix::PositionScoringMatrix(Dataset* train_dataset, int p, int q, float alpha){
    m_p = p;
    m_q = q;

    m_train_dataset = train_dataset;

    int N = train_dataset->GetNbrSequences();

    float g[26];
    for(int i=0; i<26; i++){
        g[i] = 0;
    }

    int** c = new int*[26];
    for(int i=0; i<26; i++){
        c[i] = new int[p+q];
        for(int j=0; j<p+q; j++){
            c[i][j] = 0;
        }
    }

    std::vector<std::vector<float>> f;
    f.resize(26);
    for(int i = 0 ; i < 26 ; ++i)
    {
        f[i].resize(p+q);
    }

    m_s.resize(26);
    for(int i = 0 ; i < 26 ; ++i)
    {
        m_s[i].resize(p+q);
    }

    int total = 0;
    for(int i=0; i<N; i++){
        int cleavage = train_dataset->GetCleavage(i);
        std::string sequence = train_dataset->GetSequence(i);
        for(int j=1; j<sequence.size(); j++){
            int a = sequence[j] - 'A';
            g[a] += 1;
            if(j >= cleavage - p && j < cleavage + q){
                c[a][j - cleavage + p] += 1;
            }
            total ++;
        }
    }

    for(int i=0; i<26; i++){
        g[i] = (g[i] + alpha) / (total + alpha);
    }

    for(int i=0; i<26; i++){
        for(int j=0; j<p+q; j++){
            f[i][j] = (c[i][j] + alpha) / (N + alpha);
            m_s[i][j] = log(f[i][j]) + log(g[i]);
        }
    }

    cutDataset(m_train_dataset, m_p, m_q, m_n_intances, &m_train_intances, &m_train_labels);
}

PositionScoringMatrix::~PositionScoringMatrix() {
}

void PositionScoringMatrix::getTrianScore(float& min_score, float& avg_score, float& max_score){
    int N = m_train_dataset->GetNbrSequences();

    int cleavage = m_train_dataset->GetCleavage(0);
    std::string sequence = m_train_dataset->GetSequence(0).substr(cleavage - m_p, m_p + m_q);
    float score = 0;
    for(int i=0; i<m_p + m_q; i++){
        int a = sequence[i] - 'A';
        score += m_s[a][i];
    }
    min_score = score;
    avg_score = score;
    max_score = score;

    for(int j=1; j<N; j++){
        cleavage = m_train_dataset->GetCleavage(j);
        sequence = m_train_dataset->GetSequence(j).substr(cleavage - m_p, m_p + m_q);
        score = 0;
        for(int i=0; i<m_p + m_q; i++){
            int a = sequence[i] - 'A';
            score += m_s[a][i];
        }
        if(score < min_score) min_score = score;
        if(score > max_score) max_score = score;
        avg_score += score;
    }

    avg_score = avg_score/N;
}

int PositionScoringMatrix::prediction(std::string test_data, float thresholde){
    if(test_data.size() != m_p + m_q){
        std::cout << "Prediction Error: different sizes" << std::endl;
        return -1;
    }

    float score = 0;
    for(int i=0; i<m_p + m_q; i++){
            int a = test_data[i] - 'A';
            score += m_s[a][i];
    }

    return score > thresholde;
}

int PositionScoringMatrix::Getp() const{
    return m_p;
}

int PositionScoringMatrix::Getq() const{
    return m_q;
}

void PositionScoringMatrix::matrixToSVM(const char* output_file){

    std::ofstream outfile;
    outfile.open(output_file, std::ios::out | std::ios::trunc );
    
    for(int i=0; i<m_n_intances; i++){
        std::string line;
        line += std::to_string(m_train_labels[i]);
        line += " 0:";
        line += std::to_string(i+1);
        line += " ";

        for(int j=0; j<m_n_intances; j++){
            float kernel = 0;

            for(int k=0; k< m_p+m_q; k++){
                int x = m_train_intances[i][k];
                int y = m_train_intances[j][k];

                if(x != y){
                    kernel += (m_s[x][k] + m_s[y][k]);
                    // std::cout << kernel << std::endl;
                }
                else{
                    kernel += (m_s[x][k] + log(1+exp(m_s[x][k])));
                    // std::cout << kernel << std::endl;
                }
            }
            // std::cout << kernel << std::endl;
            // kernel = exp(kernel);
            // std::cout << kernel << std::endl;

            line += std::to_string(j+1);
            line += ":";
            line += std::to_string(kernel);
            line += " ";
        }

        outfile << line << std::endl;
    }

    outfile.close();

}

void PositionScoringMatrix::matrixTestToSVM(Dataset test_dataset, const char* output_file){
    int n_intances;
    std::vector<std::vector<int>> intances;
    std::vector<int> labels;
    cutDataset(&test_dataset, m_p, m_q, n_intances, &intances, &labels);

    std::ofstream outfile;
    outfile.open(output_file, std::ios::out | std::ios::trunc );
    
    for(int i=0; i<n_intances; i++){
        std::string line;
        line += std::to_string(labels[i]);
        line += " 0:";
        line += std::to_string(i+1);
        line += " ";

        for(int j=0; j<m_n_intances; j++){
            float kernel = 0;

            for(int k=0; k< m_p+m_q; k++){
                int x = intances[i][k];
                int y = m_train_intances[j][k];

                if(x != y){
                    kernel += (m_s[x][k] + m_s[y][k]);
                    // std::cout << kernel << std::endl;
                }
                else{
                    kernel += (m_s[x][k] + log(1+exp(m_s[x][k])));
                    // std::cout << kernel << std::endl;
                }
            }
            // std::cout << kernel << std::endl;
            // kernel = exp(kernel);
            // std::cout << kernel << std::endl;

            line += std::to_string(j+1);
            line += ":";
            line += std::to_string(kernel);
            line += " ";
        }

        outfile << line << std::endl;
    }
    

    outfile.close();
}