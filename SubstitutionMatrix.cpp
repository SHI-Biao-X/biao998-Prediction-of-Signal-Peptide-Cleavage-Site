#include <vector>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <fstream>
#include "SubstitutionMatrix.hpp"

using namespace std;

SubstitutionMatrix::SubstitutionMatrix (const char* file) {
    nbr_alphabet = 26;

    //alphabet.resize(nbr_alphabet);
    s_matrix.resize(nbr_alphabet);
    alphabet = new int[nbr_alphabet];

    ifstream input;
    input.open(file);
    if (!input.is_open()){
        std::cout << "Error opening file !";
        exit(0);
    }
    string line;
    for(int i = 0; i < 6; i++) {    // to pass the first 6 lines that are irrelevant to the matrix
        getline(input, line);
    }

    getline(input, line);
    char *alphab = (char *)line.c_str();
    char *letter = strtok(alphab, " ");
    int ord;
    int position = 0;
    while (letter != NULL) {
        //cout << letter << ' ';
        ord = (int)(*letter - 'A');
        if (ord < nbr_alphabet && 0 <= ord)
            alphab[ord] = position;
        letter = strtok(NULL, " ");
        position++;
    }

    int row = 0;
    while (getline(input, line)) {
        char *r = (char *)line.c_str();
        char *amino = strtok(r, " ");
        //cout << amino << ' ';
        char *score = strtok(NULL, " ");
        int s;
        while (score != NULL) {
            s = (int)atof(score);
            //cout << s << ' ';
            s_matrix[row].push_back(s);
            score = strtok(NULL, " ");
        }
        row++;
        //cout << endl;
    }
}

int SubstitutionMatrix::GetNbrAlphabet() const {
    return nbr_alphabet;
}

int SubstitutionMatrix::GetLetterPosition(int i) const {
    return alphabet[i];
}
int SubstitutionMatrix::GetScore(int i, int j) const {
    return s_matrix[i][j];
}