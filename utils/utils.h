#ifndef UTILS
#define UTILS

#include <random>
#include <math.h>
#include <stdio.h>
#include "../classes/vector.h"
#include "scene_setup.h"

using namespace std;
static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);

void boxMuller(double stdev, double& x, double& y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
    y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
}

// returns a vector that is randomly placed around N
Vector random_cos(const Vector &N) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);

    double min_N = min(min(abs(N[0]), abs(N[1])), abs(N[2]));
    Vector T1;
    for (int i = 0; i < 3; i++) {
        if (abs(N[i]) == min_N) {
            switch (i) {
                case 0:
                    T1 = Vector(0, -N[2], N[1]);
                    break;
                case 1:
                    T1 = Vector(-N[2], 0, N[0]);
                    break;
                case 2:
                    T1 = Vector(-N[1], N[0], 0);
                    break;
                default:
                    break;
            }
            break;
        }
    }
    Vector T2 = cross(N, normalize(T1));
    return T1 * x + T2 * y + N * z;
}

Vector random_point_on_light_sphere(const Vector& x, const Vector& C, const double& R) {
    Vector D = normalize(x - C);
    Vector V = random_cos(D);
    return R * V + C;
}

#endif