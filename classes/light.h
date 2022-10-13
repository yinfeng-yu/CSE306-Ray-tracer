#ifndef LIGHT
#define LIGHT

#include "sphere.h"

class Light: public Sphere {
public:
    double intensity;

    Light() {};
    Light(Vector center, double radius, double intensity = 20000000000) : intensity (intensity) {
        this->center = center;
        this->radius = radius;
        this->albedo = Vector(1, 1, 1);
        this->type = Type::light;
    }
    
};

#endif