#include <vector>
#include <string>
#include "Dataset.hpp"

#ifndef POSITIONSCORINGMATRIX_HPP
#define POSITIONSCORINGMATRIX_HPP

class PositionScoringMatrix {
    public:
		PositionScoringMatrix(Dataset train_dataset, int p, int q, float alpha);

		int prediction(std::string test_data, float thresholde);

		void getTrianScore(Dataset train_dataset, float& min_score, float& avg_score, float& max_score);

		int Getp() const;

		int Getq() const;

	private:

		int m_p;
		int m_q;

        std::vector<std::vector<float>> m_s;

};
#endif //POSITIONSCORINGMATRIX_HPP