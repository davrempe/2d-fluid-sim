#include "SimUtil.h"

#include <iostream>

namespace SimUtil {

	
	template<typename T>
	T** initMat2D(int m, int n) {
		T** mat = new T*[m];
		for (int i = 0; i < m; i++) {
			mat[i] = new T[n];
		}

		return mat;
	}

	template<typename T>
	void deleteMat2D(int m, int n, T ** mat) {
		for (int i = 0; i < m; i++) {
			delete[] mat[i];
		}
		delete[] mat;
	}

	template <typename T>
	void printMat2D(int m, int n, T** mat) {
		for (int i = 0; i < m; i++) {
			std::cout << "[";
			for (int j = 0; j < n - 1; j++) {
				std::cout << mat[i][j] << ", ";
			}
			std::cout << mat[i][n - 1] << "]" << std::endl;
		}
	}

	template<typename T>
	T*** initMat3D(int m, int n, int l) {
		T*** mat = new T**[m];
		for (int i = 0; i < m; i++) {
			mat[i] = new T*[n];
			for (int j = 0; j < l; j++) {
				mat[i][j] = new T[l];
			}
		}

		return mat;
	}

	template<typename T>
	void deleteMat3D(int m, int n, int l, T *** mat) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				delete[] mat[i][j];
			}
			delete[] mat[i];
		}
		delete[] mat;
	}

}

// explicit instantiation of template functions for compilation
template int** SimUtil::initMat2D<int>(int, int);
template float** SimUtil::initMat2D<float>(int, int);
template void SimUtil::deleteMat2D<int>(int, int, int**);
template void SimUtil::deleteMat2D<float>(int, int, float**);
template void SimUtil::printMat2D<int>(int, int, int**);
template void SimUtil::printMat2D<float>(int, int, float**);
template int*** SimUtil::initMat3D<int>(int, int, int);
template float*** SimUtil::initMat3D<float>(int, int, int);
template void SimUtil::deleteMat3D<int>(int, int, int, int***);
template void SimUtil::deleteMat3D<float>(int, int, int, float***);