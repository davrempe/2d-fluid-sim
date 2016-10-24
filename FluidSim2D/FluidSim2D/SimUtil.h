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

	struct Vec3 {
		float x, y, z;
		Vec3(float x, float y, float z) : x(x), y(y), z(z) {}
	};

	struct Particle2D {
		Vec2 pos;
		Vec2 vel;
		Particle2D(Vec2 pos, Vec2 vel) : pos(pos), vel(vel) {}
	};

	struct Particle3D {
		Vec3 pos;
		Vec3 vel;
		Particle3D(Vec3 pos, Vec3 vel) : pos(pos), vel(vel) {}
	};

	//----------------------------------------------------------------------
	// Functions
	//----------------------------------------------------------------------

	/*
	Initializes a 2D matrix.
	Args:
	m - num rows
	n - num cols
	Returns:
	T** - the dynamic array
	*/
	template <typename T> T** initMat2D(int, int);
	/*
	Deletes a 2D matrix.
	Args:
	m - num rows
	n - num cols
	mat - matrix to delete
	*/
	template <typename T> void deleteMat2D(int, int, T**);
	/*
	Prints the given matrix to stdout.
	Args:
	m - num rows
	n - num cols
	T** - matrix to print
	*/
	template <typename T> void printMat2D(int, int, T**);
	/*
	Initializes a 3D matrix.
	Args:
	m - num rows
	n - num cols
	l - depth
	Returns:
	T*** - the dynamic array
	*/
	template <typename T> T*** initMat3D(int, int, int);
	/*
	Deletes a 3D matrix.
	Args:
	m - num rows
	n - num cols
	l - depth
	mat - matrix to delete
	*/
	template <typename T> void deleteMat3D(int, int, int, T***);

	// vec operations
}

#endif //SIM_UTIL_H

