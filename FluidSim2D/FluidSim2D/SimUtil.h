#ifndef SIM_UTIL_H
#define SIM_UTIL_H

namespace SimUtil {
	//----------------------------------------------------------------------
	// Constants
	//----------------------------------------------------------------------

	const int SOLID = 0;
	const int FLUID = 1;
	const int AIR = 2;

	//----------------------------------------------------------------------
	// Data Structures
	//----------------------------------------------------------------------

	struct Vec2 {
		float x, y;
		Vec2(float x, float y) : x(x), y(y) {}
	};

	// particle struct

	//----------------------------------------------------------------------
	// Functions
	//----------------------------------------------------------------------

	/*
	Initializes a 2D matrix.
	Args:
	x - num rows
	y - num cols
	*/
	template <typename T> void initMat2D(int, int, T**);
	/*
	Deletes a 2D matrix.
	Args:
	x - num rows
	y - num cols
	*/
	template <typename T> void deleteMat2D(int, int, T**);
	/*
	Initializes a 3D matrix.
	Args:
	x - num rows
	y - num cols
	z - depth
	*/
	template <typename T> void initMat3D(int, int, int, T***);
	/*
	Deletes a 3D matrix.
	Args:
	x - num rows
	y - num cols
	z - depth
	*/
	template <typename T> void deleteMat3D(int, int, int, T***);

	// vec operations
}

#endif //SIM_UTIL_H

