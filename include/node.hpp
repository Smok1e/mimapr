#pragma once

#include <cstdint>

//========================================

struct Node
{
	enum class Type: uint8_t
	{
		Internal,
		External
	} type = Type::External;

	enum class BoundaryCondition: uint8_t
	{
		None,
		Dirichlet,
		Neumann,
		Neuton
	} bc = BoundaryCondition::None;

	double x = 0;
	double y = 0;
	double T = 0;
	double bc_value = 0;

	int id = -1;
};

//========================================
