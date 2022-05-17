#include <iostream>
#include <vector>

using namespace std;


#ifndef SUBSTITUTIONMTRIX_HPP
#define SUBSTITUTIONMTRIX_HPP

class SubstitutionMatrix {
    public :
        /**
			The constructor needs the path of the file as a string.
		*/ 
        SubstitutionMatrix(const char* file);

        /**
			Standard destructor
		*/
        //~SubstitutionMatrix();

        int GetNbrAlphabet() const;

        int GetScore(int i, int j) const;

        int GetLetter(int i) const;

    private : 
        /**
            The number of letters
        */
        int nbr_alphabet;

        /**
            The substitution matrix
        */
        vector<vector<int> > s_matrix;

        /**
            The order of the letters in the alphabet
        */
        vector<int> alphabet;

};

#endif