#ifndef RAY
#define RAY

#include "vector.h"

using namespace std;

class Ray {
public:
    Vector origin;
    Vector unit;
    double time;

    Ray(Vector origin, Vector unit, double time = 0.) : origin(origin), unit(unit), time(time) { };
};

#endif