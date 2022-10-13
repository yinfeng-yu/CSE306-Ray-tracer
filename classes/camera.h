#ifndef CAMERA
#define CAMERA

#include "geometry.h"
#include "../utils/utils.h"

using namespace std;

class Camera {
public:
    Vector center;
    double focus_length;

    Camera() {}
    Camera(Vector center, double f) {
        this->center = center;
        this->focus_length = f;
    }
    
    Ray make_ray(int W, int H, int i, int j, double focus_distance = 50., double spread = 0.2, double shutter_time = 0.2) {
        double Qx = center[0];
        double Qy = center[1];
        double Qz = center[2];
        double offset_x, offset_y;

        Vector ray_dir(Qx + j + 0.5 - W / 2, Qy + i + 0.5 - H / 2, Qz - focus_length);

        boxMuller(1, offset_x, offset_y);
        ray_dir[0] += offset_x * spread;
        ray_dir[1] += offset_y * spread;

        Vector u = normalize(ray_dir - center);

#ifdef REGULAR
        return Ray(center, u);
#endif
#ifdef DOF
        Vector P = get_point_in_focus(focus_distance, u);
        Vector Q_prime = center + random_point_on_unit_disk(); // new ray origin from a random position on the unit disk around the camera center
        Vector u_prime = normalize(P - Q_prime);
        return Ray(Q_prime, u_prime);
#endif
#ifdef BLUR
        return Ray(center, u, uniform(engine) * shutter_time);
#endif
        
    }

    Vector get_point_in_focus(double D, Vector u) {
        return center + (D / abs(u[2])) * u;
    }

    Vector random_point_on_unit_disk() {
        double r = sqrt(uniform(engine));
        double theta = 2 * M_PI * uniform(engine);
        return Vector(r * cos(theta), r * sin(theta), 0.);
    }


};

#endif