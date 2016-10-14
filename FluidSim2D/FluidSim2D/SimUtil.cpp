#include "SimUtil.h"

namespace SimUtil {

	
	template<typename T>
	void initMat2D(int x, int y, T ** mat) {
		mat = new T*[x];
		for (int i = 0; i < x; i++) {
			mat[i] = new T[y];
		}
	}

	template<typename T>
	void deleteMat2D(int x, int y, T ** mat) {
		for (int i = 0; i < x; i++) {
			delete[] mat[i];
		}
		delete[] mat;
	}

	template<typename T>
	void initMat3D(int x, int y, int z, T *** mat) {
		mat = new T**[x];
		for (int i = 0; i < x; i++) {
			mat[i] = new T*[y];
			for (int j = 0; j < z; j++) {
				mat[i][j] = new T[z];
			}
		}
	}

	template<typename T>
	void deleteMat3D(int x, int y, int z, T *** mat) {
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < y; j++) {
				delete[] mat[i][j];
			}
			delete[] mat[i];
		}
		delete[] mat;
	}
}