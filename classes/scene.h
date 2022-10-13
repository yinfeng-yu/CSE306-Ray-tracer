#ifndef SCENE
#define SCENE

#include <string>
#include <iostream>
#include "vector.h"
#include "light.h"

using namespace std;
#define EPSILON 0.001

class Scene {
public: 
    Light* light;
    vector<Sphere*> collection;

    Scene() {
#ifdef REGULAR
        light = new Light(Vector(-10, 25, -10), 5, 20000000000);
#endif
#ifdef DOF
        light = new Light(Vector(-10, 20, 40), 5, 20000000000);
#endif
#ifdef BLUR
        light = new Light(Vector(-10, 20, 40), 5, 20000000000);
#endif

        static Sphere* collection[] = {

#ifdef REGULAR
            new Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1), Type::transparent),                    // center sphere
            new Sphere(Vector(-20, 0, 0), 10, Vector(1, 1, 1), Type::mirror),                       // left sphere
            new Sphere(Vector(20, 0, 0), 10, Vector(1, 1, 1), Type::transparent),                   // right sphere
            new Sphere(Vector(20, 0, 0), 9.5, Vector(1, 1, 1), Type::hollow),                       // right_inner_hollow
#endif
#ifdef DOF
            new Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1), Type::diffuse),                        // center sphere
            new Sphere(Vector(-10, 0, 25), 10, Vector(1, 1, 1), Type::transparent),                 // front sphere
            new Sphere(Vector(10, 0, -25), 10, Vector(1, 1, 1), Type::mirror),                      // back sphere
#endif
#ifdef BLUR
            new Sphere(Vector(0, 0, 0), 10, Vector(1, 1, 1), Type::diffuse, Vector(0., 50., 0.)),   // center sphere
#endif
            new Sphere(Vector(0,  1000, 0), 940, Vector(1, 0, 0), Type::diffuse),                   // red wall
            new Sphere(Vector(0, 0, -1000), 940, Vector(0, 1, 0), Type::diffuse),                   // green wall
            new Sphere(Vector(0, -1000, 0), 990, Vector(0, 0, 1), Type::diffuse),                   // blu wall
            new Sphere(Vector(0, 0,  1000), 940, Vector(1, 0, 1), Type::diffuse),                   // pink wall
            new Sphere(Vector( 1000, 0, 0), 940, Vector(0, 1, 1), Type::diffuse),                   // cyan wall
            new Sphere(Vector(-1000, 0, 0), 940, Vector(1, 1, 0), Type::diffuse),                   // yellow wall

            light,
        };

        this->collection = vector<Sphere *>(collection, collection + sizeof(collection) / sizeof(collection[0]));
    }

    // finds the nearest sphere which intersects the ray
    bool intersect(const Ray& ray, double& t_near, int& id) {
        double t = 0.;
        bool res = false;

        for (int i = 0; i < collection.size(); i++) {
            Geometry* cur = this->collection[i];
            if (cur->intersect(ray, t)) {
                if (t < t_near) {
                    t_near = t;
                    id = i;
                    res = true;
                }
            }
        }
        return res;
    }

    Vector get_color(const Ray& ray, const int ray_depth, const bool last_bounce_diffuse) {
        // terminates recursion at some point
        if (ray_depth < 0) return Vector(); 
    
        double t_near = INFINITY;
        int id = 0;

        if (intersect(ray, t_near, id)) {
            const Vector albedo = collection[id]->albedo;
            Vector P = ray.origin + t_near * ray.unit;
            Vector N = normalize(P - collection[id]->center);
            double I = light->intensity;

            // distance from P to the light
            double R = norm(light->center - P); 

            if (collection[id]->type == Type::light) {
                // if this is an indirect diffuse bounce
                if (last_bounce_diffuse) { 
                    // if we hit a light source by chance via an indirect diffuse bounce, return 0 to avoid counting it twice
                    return Vector(0., 0., 0.);
                }
                else {
                    return Vector(1., 1., 1.) * I / (4 * M_PI * M_PI * R * R);
                }
            }
            else if (collection[id]->type == Type::diffuse) {
                // handle diffuse surfaces
                Vector Lo(0., 0., 0.);

                // add direct lighting
                Vector x_prime = random_point_on_light_sphere(P, light->center, light->radius);

                Vector N_prime = normalize(x_prime - light->center);
                Vector omega_i = normalize(x_prime - P);
                double pdf = dot(N_prime, normalize(P - light->center)) / (M_PI * R * R);

                Ray shadow_ray(P + N * EPSILON, omega_i); // computes the visibility term by launching a ray of direction omega_i

                double visibility = 0;
                double shadow_test_t_near = INFINITY;
                int shadow_test_id = 0;

                if (intersect(shadow_ray, shadow_test_t_near, shadow_test_id)) {
                    visibility = collection[shadow_test_id]->type == Type::light ? 1. : 0.;
                }

                Lo = ((I * albedo) / (4 * M_PI * M_PI * M_PI * R * R)) * visibility * max(dot(N, omega_i), 0.)
                    * max(dot(N_prime, -omega_i), 0.) / (squared_norm(x_prime - P) * pdf);

                // add indirect lighting
                Ray random_ray = Ray(P + N * EPSILON, normalize(random_cos(N))); // randomly sample ray using random_cos

                // albedo * get_color(random_ray, ray_depth - 1, true) seems to be not thread-safe (some undetermined behaviours)
                const Vector indirect_color = get_color(random_ray, ray_depth - 1, true);

                Lo += albedo * indirect_color;

                return Lo;
            }

            else if (collection[id]->type == Type::mirror) {
                Vector omega_i = ray.unit;
                Vector omega_r = omega_i - N * (2 * dot(N, omega_i));
                Ray refl = Ray(P + N * EPSILON, normalize(omega_r));
                return get_color(refl, ray_depth - 1, false);
            }

            else if (collection[id]->type == Type::transparent) {
                double n1 = 1.;
                double n2 = 1.5;
                Vector omega_i = ray.unit;

                // exiting the sphere
                if (dot(omega_i, N) > 0) {
                    n1 = 1.5;
                    n2 = 1.;

                    N = N * -1;
                }
                // entering the sphere
                else {
                    n1 = 1.;
                    n2 = 1.5;
                }

                // Fresnel law
                double k0 = pow(((n1 - n2) / (n1 + n2)), 2);
                double R = k0 + (1 - k0) * pow((1 - abs(dot(N, omega_i))), 5);

                // (uniform) random number between 0 and 1
                double u = (double)rand() / (RAND_MAX); 

                // launch a reflection ray
                if (u < R) { 
                    Vector omega_r = omega_i - N * (2 * dot(N, omega_i));
                    Ray refl = Ray(P + N * EPSILON, normalize(omega_r));
                    return get_color(refl, ray_depth - 1, false);
                }

                // launch a refraction ray
                else {
                    double radicand = 1 - (((n1 * n1) / (n2 * n2)) * (1 - (dot(omega_i, N) * dot(omega_i, N))));

                    // Snell_Decartes law
                    // radicand < 0, only happens when n1 > n2 (exiting the sphere) and lead to total reflection
                    if (radicand < 0) {
                        // N (normal) is already inverted
                        Vector omega_r = omega_i - N * (2 * dot(N, omega_i));
                        Ray refl = Ray(P + N * EPSILON, normalize(omega_r));
                        return get_color(refl, ray_depth - 1, false);
                    }

                    // radicand > 0, refraction
                    else {
                        Vector omega_t_T = (omega_i - N * dot(omega_i, N)) * (n1 / n2);
                        Vector omega_t_N = -1 * N * sqrt(radicand);
                        Vector omega_t = omega_t_N + omega_t_T;
                        Ray refr = Ray(P - N * EPSILON, normalize(omega_t));
                        return get_color(refr, ray_depth - 1, false);
                    }
                }
            }

            else if (collection[id]->type == Type::hollow) {
                double n1 = 1.;
                double n2 = 1.5;
                Vector omega_i = ray.unit;

                // (Hack) to simulate hollow, invert the normal
                N = -N;

                // exiting the sphere
                if (dot(omega_i, N) > 0) {
                    n1 = 1.5;
                    n2 = 1.;
                    N = -N;
                }
                // entering the sphere
                else {
                    n1 = 1.;
                    n2 = 1.5;
                }

                // Fresnel law
                double k0 = pow(((n1 - n2) / (n1 + n2)), 2);
                double R = k0 + (1 - k0) * pow((1 - abs(dot(N, omega_i))), 5);

                // (uniform) random number between 0 and 1
                double u = (double)rand() / (RAND_MAX); 
                
                // launch a reflection ray
                if (u < R) { 
                    // N (normal) is already inverted
                    // to reflect we need to revert it back (to original)
                    N = -N;
                    Vector omega_r = omega_i - N * (2 * dot(N, omega_i));
                    Ray refl = Ray(P + N * EPSILON, normalize(omega_r));
                    return get_color(refl, ray_depth - 1, false);
                }

                // launch a refraction ray
                else {
                    double radicand = 1 - (((n1 * n1) / (n2 * n2)) * (1 - (dot(omega_i, N) * dot(omega_i, N))));

                    // Snell_Decartes law
                    // radicand < 0, only happens when n1 > n2 (exiting the sphere) and lead to total reflection
                    if (radicand < 0) {
                        // N (normal) is already inverted
                        // to reflect we need to revert it back (to original)
                        N = -N;
                        Vector omega_r = omega_i - N * (2 * dot(N, omega_i));
                        Ray refl = Ray(P + N * EPSILON, normalize(omega_r));
                        return get_color(refl, ray_depth - 1, false);
                    }

                    // radicand > 0, refraction
                    else {
                        Vector omega_t_T = (omega_i - N * dot(N, omega_i)) * (n1 / n2);
                        Vector omega_t_N = -1 * N * sqrt(radicand);
                        Vector omega_t = omega_t_N + omega_t_T;
                        Ray refr = Ray(P - N * EPSILON, normalize(omega_t));
                        return get_color(refr, ray_depth - 1, false);
                    }
                }
            }
            
        }
        
    }

};

#endif