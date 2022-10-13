#ifndef SPHERE
#define SPHERE

#include "geometry.h"
#include "../utils/utils.h"

using namespace std;

class Sphere: public Geometry {
public: 
    Vector center;
    double radius;

    Sphere() {
        center = Vector();
        radius = 1.;
    };
    Sphere(Vector center, double radius, Vector albedo, Type type, Vector velocity = Vector(0., 0., 0.)) : center(center), radius(radius) {
        this->albedo = albedo;
        this->type = type;
        this->velocity = velocity;
    }
    
    bool intersect(const Ray& ray, double& t) {
        Vector O = ray.origin;
        Vector U = ray.unit;

#ifdef BLUR
        Vector C = this->center + ray.time * velocity;
#else
        Vector C = this->center;
#endif
        double R = this->radius;
        double delta = (dot(U, O - C) * dot(U, O - C)) - (norm(O - C) * norm(O - C)) + (R * R);

        if (delta < 0) {
            return false;
        }
        else {
            double t_1 = dot(U, C - O) - sqrt(delta);
            double t_2 = dot(U, C - O) + sqrt(delta);

            if (t_2 < 0) {
                return false;
            }
            else {
                t = t_1 >= 0 ? t_1 : t_2;
                return true;
            }
        }
    }

};

#endif