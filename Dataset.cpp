#include <iostream>
#include <fstream>
#include "Dataset.hpp"

Dataset::Dataset(const char* file) {
	m_nseq = 0;

	std::ifstream fin(file);
	
	if (fin.fail()) {
		std::cout<<"Cannot read from file "<<file<<" !"<<std::endl;
		exit(1);
	}
	
    std::string line; 
    int row = 1;

	while (getline(fin, line)) {
        
        switch(row){
        case 1:
            row = 2;
            break;
        case 2:
            m_sequences.push_back(line);
            row = 3;
		    m_nseq ++;
            break;
        case 3:
            for (int i = 0; i < line.size(); i++) {
                if(line[i] == 'C') m_cleavages.push_back(i);
            }
            row = 1;
            break;
        }
	}
	
	fin.close();
}

Dataset::Dataset(const Dataset& data, bool is_train, float train_rate){
    if(is_train){
        m_nseq = int(data.GetNbrSequences() * train_rate);
        m_cleavages.assign(data.GetCleavages().begin(), data.GetCleavages().begin() + m_nseq);
        m_sequences.assign(data.GetSequences().begin(), data.GetSequences().begin() + m_nseq);        
    }
    else{
        m_nseq = data.GetNbrSequences() - int(data.GetNbrSequences() * train_rate);
        m_sequences.assign(data.GetSequences().begin() + m_nseq, data.GetSequences().end());
        m_cleavages.assign(data.GetCleavages().begin() + m_nseq, data.GetCleavages().end());
    }
}

Dataset::~Dataset() {
}

int Dataset::GetNbrSequences() const {
	return m_nseq;
}

const std::string& Dataset::GetSequence(int i) const {
	return m_sequences[i];
}

const int Dataset::GetCleavage(int i) const {
	return m_cleavages[i];
}

const std::vector<std::string> & Dataset::GetSequences() const {
    return m_sequences;
}

const std::vector<int> & Dataset::GetCleavages() const{
    return m_cleavages;
}