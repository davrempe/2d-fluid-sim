#ifndef SIM_UTIL_H
#define SIM_UTIL_H

//----------------------------------------------------------------------
// Constants
//----------------------------------------------------------------------

enum class CellType {
	SOLID,
	FLUID,
	AIR
};

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

// vec operations

#endif //SIM_UTIL_H

