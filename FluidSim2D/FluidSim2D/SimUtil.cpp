#include "SimUtil.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

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
	T** initGrid2D(int x, int y) {
		return initMat2D<T>(x, y);
	}

	template<typename T>
	void deleteGrid2D(int x, int y, T ** grid) {
		deleteMat2D<T>(x, y, grid);
	}

	template <typename T>
	void printGrid2D(int x, int y, T** grid) {
		std::cout << "====================================================================================\n";
		for (int i = y - 1; i >= 0; i--) {
			std::cout << "[";
			for (int j = 0; j < x - 1; j++) {
				std::cout << grid[j][i] << ", ";
			}
			std::cout << grid[x - 1][i] << "]" << std::endl;
		}
		std::cout << "====================================================================================\n";
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

	template<typename T>
	T*** initGrid3D(int x, int y, int z) {
		return initMat3D<T>(x, y, z);
	}

	template<typename T>
	void deleteGrid3D(int x, int y, int z, T *** grid) {
		deleteMat3D<T>(x, y, z, grid);
	}

	void readInGeom(int x, int y, std::string geomFileName, Mat2Di grid) {
		// open the geometry file
		std::ifstream geomFile(geomFileName);
		if (geomFile.is_open()) {
			std::string lineStr;
			// parse file based on given dimensions, will error if file does not match these
			// fills grid so that [0][0] is at bottom left corner of simulation
			for (int i = y - 1; i >= 0; i--) {
				std::getline(geomFile, lineStr);
				for (int j = 0; j < x; j++) {
					switch (lineStr[j]) {
					case 'f':
						grid[j][i] = SimUtil::FLUID;
						break;
					case 's':
						grid[j][i] = SimUtil::SOLID;
						break;
					case 'a':
						grid[j][i] = SimUtil::AIR;
						break;
					}
				}
			}
			geomFile.close();
		}
	}

	Vec2 getGridCellPosition(float i, float j, float dx) {
		return Vec2(i*dx + 0.5f*dx, j*dx + 0.5f*dx);
	}

	int* getGridCellIndex(Vec2 pos, float dx) {
		int index[2] = { (int)(pos.x / dx), (int)(pos.y / dx) };
		return index;
	}

	Vec2 add(Vec2 vec1, Vec2 vec2) {
		return Vec2(vec1.x + vec2.x, vec1.y + vec2.y);
	}

	Vec3 add(Vec3 vec1, Vec3 vec2) {
		return Vec3(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z);
	}

	Vec2 sub(Vec2 vec1, Vec2 vec2) {
		return Vec2(vec1.x - vec2.x, vec1.y - vec2.y);
	}
	
	Vec3 sub(Vec3 vec1, Vec3 vec2) {
		return Vec3(vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z);
	}

	Vec2 scale(Vec2 vec1, float scalar) {
		return Vec2(scalar * vec1.x, scalar * vec1.y);
	}

	Vec3 scale(Vec3 vec1, float scalar) {
		return Vec3(scalar * vec1.x, scalar * vec1.y, scalar * vec1.z);
	}

	float norm(Vec2 vec) {
		return sqrt(vec.x * vec.x + vec.y * vec.y);
	}

	float norm(Vec3 vec) {
		return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
	}

	template <typename T>
	double dot(T** grid1, T** grid2, int x, int y) {
		double dotProd = 0.0;
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < y; j++) {
				dotProd += grid1[i][j] * grid2[i][j];
			}
		}

		return dotProd;
	}

	template <typename T>
	T max(T** grid1, int x, int y) {
		T maxVal = std::numeric_limits<T>::lowest();
		for (int i = 0; i < x; i++) {
			for (int j = 0; j < y; j++) {
				if (grid1[i][j] > maxVal) {
					maxVal = grid1[i][j];
				}
			}
		}

		return maxVal;
	}

}

// explicit instantiation of template functions for compilation
template int** SimUtil::initMat2D<int>(int, int);
template float** SimUtil::initMat2D<float>(int, int);
template double** SimUtil::initMat2D<double>(int, int);
template void SimUtil::deleteMat2D<int>(int, int, int**);
template void SimUtil::deleteMat2D<float>(int, int, float**);
template void SimUtil::deleteMat2D<double>(int, int, double**);
template void SimUtil::printMat2D<int>(int, int, int**);
template void SimUtil::printMat2D<float>(int, int, float**);
template int** SimUtil::initGrid2D<int>(int, int);
template float** SimUtil::initGrid2D<float>(int, int);
template double** SimUtil::initGrid2D<double>(int, int);
template void SimUtil::deleteGrid2D<int>(int, int, int**);
template void SimUtil::deleteGrid2D<float>(int, int, float**);
template void SimUtil::deleteGrid2D<double>(int, int, double**);
template void SimUtil::printGrid2D<int>(int, int, int**);
template void SimUtil::printGrid2D<float>(int, int, float**);
template void SimUtil::printGrid2D<double>(int, int, double**);
template int*** SimUtil::initMat3D<int>(int, int, int);
template float*** SimUtil::initMat3D<float>(int, int, int);
template double*** SimUtil::initMat3D<double>(int, int, int);
template void SimUtil::deleteMat3D<int>(int, int, int, int***);
template void SimUtil::deleteMat3D<float>(int, int, int, float***);
template void SimUtil::deleteMat3D<double>(int, int, int, double***);
template int*** SimUtil::initGrid3D<int>(int, int, int);
template float*** SimUtil::initGrid3D<float>(int, int, int);
template double*** SimUtil::initGrid3D<double>(int, int, int);
template void SimUtil::deleteGrid3D<int>(int, int, int, int***);
template void SimUtil::deleteGrid3D<float>(int, int, int, float***);
template void SimUtil::deleteGrid3D<double>(int, int, int, double***);

template double SimUtil::dot(int**, int**, int, int);
template double SimUtil::dot(float**, float**, int, int);
template double SimUtil::dot(double**, double**, int, int);

template int SimUtil::max(int**, int, int);
template float SimUtil::max(float**, int, int);
template double SimUtil::max(double**, int, int);