#ifndef SIM_UTIL_H
#define SIM_UTIL_H

#include <string>

namespace SimUtil {
	typedef int** Mat2Di;
	typedef float** Mat2Df;
	typedef double** Mat2Dd;
	typedef int*** Mat3Di;
	typedef float*** Mat3Df;
	typedef double*** Mat3Dd;

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
	Initializes a 2D grid. This is done in column major so the grid can lie on
	a cartesian coordinate system so that (0,0) is at the bottom left corner of grid cell [0][0]. And
	generally grid cells can be accessed [x][y].
	Args:
	x - x dimension width (in number of cells)
	y - y dimension height (in number of cells)
	Returns:
	T** - the dynamic array
	*/
	template <typename T> T** initGrid2D(int, int);
	/*
	Deletes a 2D grid.
	Args:
	x - x dimension width
	y - y dimension height
	grid - grid to delete
	*/
	template <typename T> void deleteGrid2D(int, int, T**);
	/*
	Prints the given grid to stdout.
	Args:
	x - x dimension width
	n - y dimension height
	T** - grid to print
	*/
	template <typename T> void printGrid2D(int, int, T**);

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

	/*
	Initializes a 3D grid. This is done in column major so the grid can lie on
	a cartesian coordinate system so that (0,0,0) is at the bottom left back corner of grid cell [0][0][0]. And
	generally grid cells can be accessed [x][y][z].
	Args:
	x - x dimension width (in grid cells)
	y - y dimension height (in grid cells)
	z - z dimension depth (in grid cells)
	Returns:
	T** - the dynamic array
	*/
	template <typename T> T*** initGrid3D(int, int, int);
	/*
	Deletes a 3D grid.
	Args:
	x - x dimension width
	y - y dimension height
	z - z dimension depth
	grid - grid to delete
	*/
	template <typename T> void deleteGrid3D(int, int, int, T***);

	/*
	Builds initial grid of dimensions (x, y) that contains the initial
	geometry for the system to simulate. Cell [0][0] in the grid is at the bottom
	left corner of the input geometry, so it's treated as if the input grid was initialized
	using initGrid2D. x is positive right, y is positive up. 
	It reads the initial geometry from the specified input file parameter.
	Args:
	x, y - grid dimensions in number of cells
	geomFile - the file containing the geometry
	grid - the 2D array to put the initial grid in
	*/
	void readInGeom(int, int, std::string, SimUtil::Mat2Di);

	/*
	Finds the physical location of the cell with index [x][y]
	in a grid where [0][0] is the center of the bottom left cell, based on the given dx. 
	Integer indices are treated in the center of cells while fractional indices may lie anywhere
	in a grid cell.
	Args:
	i - x index of cell
	j - y index of cell
	dx - single cell dimension
	Returns:
	Vec2 (x, y) containing physical location from bottom left corner of grid.
	*/
	Vec2 getGridCellPosition(float, float, float);
	/*
	Returns array of size 2 that contains the integer grid cell with index [i][j] at its center that contains the given position.
	Args:
	pos - a Vec2 (x,y) coordinate containing the position to use based on the origin at the bottom left of the simulation.
	dx - single cell dimension
	*/
	int* getGridCellIndex(Vec2 pos, float);

	// vec operations

	/*
	Calculates the sum of vec1 and vec2 (vec1 + vec2) and returns a new vector containing this subtraction.
	Args
	vec1 - first vector
	vec2 - second vector
	*/
	Vec2 add(Vec2, Vec2);
	/*
	Calculates the sum of vec1 and vec2 (vec1 + vec2) and returns a new vector containing this subtraction.
	Args
	vec1 - first vector
	vec2 - second vector
	*/
	Vec3 add(Vec3, Vec3);
	/*
	Calculates the difference of vec1 and vec2 (vec1 - vec2) and returns a new vector containing this subtraction.
	Args
	vec1 - first vector
	vec2 - second vector
	*/
	Vec2 sub(Vec2, Vec2);
	/*
	Calculates the difference of vec1 and vec2 (vec1 - vec2) and returns a new vector containing this subtraction.
	Args
	vec1 - first vector
	vec2 - second vector
	*/
	Vec3 sub(Vec3, Vec3);
	/*
	Scales the given vector by a scalar and returns the new scaled vector. 
	Args
	vec1 - the vector
	scalar - the value to scale by
	*/
	Vec2 scale(Vec2, float);
	/*
	Scales the given vector by a scalar and returns the new scaled vector.
	Args
	vec1 - the vector
	scalar - the value to scale by
	*/
	Vec3 scale(Vec3, float);
	/*
	Calculates Euclidean norm of the given vector.
	Args
	vec - the vector to calculate norm of.
	*/
	float norm(Vec2);
	/*
	Calculates Euclidean norm of the given vector.
	Args
	vec - the vector to calculate norm of.
	*/
	float norm(Vec3);
	/*
	Calculates the dot product of two vectors being stored
	as 2D grids (grid1 dot grid2).
	Args:
	grid1 - the first vector (2D grid)
	grid2 - the second vector (2D grid)
	x/y - grid dimensions
	*/
	template <typename T> double dot(T**, T**, int, int);
	/*
	Finds the maximum value in a grid.
	Args:
	grid1 - the grid to look in
	x/y - grid dimensions
	*/
	template <typename T> T max(T**, int, int);
}

#endif //SIM_UTIL_H

