#include <iostream>
#include <string>
#include <vector>

#ifndef DATASET_HPP
#define DATASET_HPP

class Dataset {
    public:
		/**
			The constructor needs the path of the file as a string.
		*/
		Dataset(const char* file);

		Dataset(const Dataset& data, bool is_train, float train_rate);

		/**
			Standard destructor
		*/
		~Dataset();

		/**
          The getter to the number of sequences.
        */
    	int GetNbrSequences() const;

		/**
        Returns a copy of a sequence.
        @param i Sequence number to get.
      	*/
    	const std::string & GetSequence(int i) const;

		const int GetCleavage(int i) const;

		const std::vector<std::string> & GetSequences() const;

		const std::vector<int> & GetCleavages() const;

	private:
        /**
          The number of sequences.
        */
		int m_nseq;

		/**
          The sequence is stored as a vector of string.
        */
        std::vector<std::string> m_sequences;

		/**
          The cleavage is stored as a vector of int.
        */
        std::vector<int> m_cleavages;
};
#endif //DATASET_HPP