#include <vector>
#include <string>
#include "Dataset.hpp"

#ifndef POSITIONSCORINGMATRIX_HPP
#define POSITIONSCORINGMATRIX_HPP

class PositionScoringMatrix {
    public:
		PositionScoringMatrix(Dataset* train_dataset, int p, int q, float alpha);

		~PositionScoringMatrix();

		int prediction(std::string test_data, float thresholde);

		void getTrianScore(float& min_score, float& avg_score, float& max_score);

		int Getp() const;

		int Getq() const;

		void matrixToSVM(const char* output_file);

		void matrixTestToSVM(Dataset test_dataset, const char* output_file);

	private:

		int m_p;
		int m_q;

        std::vector<std::vector<float>> m_s;

		Dataset* m_train_dataset;

		int m_n_intances;
		std::vector<std::vector<int>> m_train_intances;
		std::vector<int> m_train_labels;

};
#endif //POSITIONSCORINGMATRIX_HPP