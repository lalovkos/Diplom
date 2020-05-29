#include "Point.h"

std::string Point::ToString(){
	return "{" + std::to_string(this->x) + "," + std::to_string(this->y) + "," + std::to_string(this->z) + "}";
}