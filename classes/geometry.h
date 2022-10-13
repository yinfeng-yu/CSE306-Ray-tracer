#ifndef GEOMETRY
#define GEOMETRY

#include "ray.h" 

enum class Type {
    diffuse,
    mirror,
    transparent,
    hollow,
    light,
};

class Geometry {
public:
	Vector albedo;
	Type type;
	Vector velocity;

	virtual bool intersect(const Ray& ray, double& t) = 0;
};

#endif