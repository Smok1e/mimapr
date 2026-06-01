#pragma once

#include <cstdint>

//========================================

struct Node
{
	enum class Type: uint8_t
	{
		Internal,
		External,
		Boundary
	} type = Type::External;

	enum class BoundaryConditon: uint8_t
	{
		None,
		Dirichlet, // T
		Neumann    // dT/dn
	} boundary_condition = BoundaryConditon::None;

	double x = 0;
	double y = 0;
	double T = 0;

	int id = -1;
};

//========================================