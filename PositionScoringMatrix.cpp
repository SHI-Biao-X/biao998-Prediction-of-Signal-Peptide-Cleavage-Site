#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "PositionScoringMatrix.hpp"

PositionScoringMatrix::PositionScoringMatrix(Dataset train_dataset, int p, int q, float alpha){
    m_p = p;
    m_q = q;

    int N = train_dataset.GetNbrSequences();

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

    float** f = new float*[26];
    for(int i=0; i<26; i++){
        f[i] = new float[p+q];
    }

    m_s.resize(26);
    for(int i = 0 ; i < 26 ; ++i)
    {
        m_s[i].resize(p+q);
    }

    int total = 0;
    for(int i=0; i<N; i++){
        int cleavage = train_dataset.GetCleavage(i);
        std::string sequence = train_dataset.GetSequence(i);
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
}

void PositionScoringMatrix::getTrianScore(Dataset train_dataset, float& min_score, float& avg_score, float& max_score){
    int N = train_dataset.GetNbrSequences();

    int cleavage = train_dataset.GetCleavage(0);
    std::string sequence = train_dataset.GetSequence(0).substr(cleavage - m_p, m_p + m_q);
    float score = 0;
    for(int i=0; i<m_p + m_q; i++){
        int a = sequence[i] - 'A';
        score += m_s[a][i];
    }
    min_score = score;
    avg_score = score;
    max_score = score;

    for(int j=1; j<N; j++){
        cleavage = train_dataset.GetCleavage(j);
        sequence = train_dataset.GetSequence(j).substr(cleavage - m_p, m_p + m_q);
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