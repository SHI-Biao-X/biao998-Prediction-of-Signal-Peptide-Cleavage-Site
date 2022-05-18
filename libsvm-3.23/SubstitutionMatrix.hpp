#include <iostream>
#include <vector>


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

        int GetLetterPosition(int i) const;

    private : 
        /**
            The number of letters
        */
        int nbr_alphabet;

        /**
            The substitution matrix
        */
        std::vector<std::vector<int> > s_matrix;

        /**
            The order of the letters in the alphabet
        */
        int *alphabet;

};

#endif