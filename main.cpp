#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STBIW_WINDOWS_UTF8
#include "./include/include.h"
#include <chrono>
#include <iostream>
#include <string>
using namespace chrono;

void render(Image& image, Camera& camera, Scene& scene, int num_path, int max_ray_depth, double gamma = 2.2) {
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < image.height; ++i) {
        std::cout << "\rLines scanned: " << i << ' ' << std::flush;

        for (int j = 0; j < image.width; ++j) {
            Vector pixel_color;
            bool printed = false;
            for (int k = 0; k < num_path; ++k) {

                Ray ray = camera.make_ray(image.width, image.height, i, j);

                Vector color = scene.get_color(ray, max_ray_depth, false);
                pixel_color += color;
            }
            pixel_color = pixel_color / num_path;

            image.pixels[(image.height - i - 1) * image.width * 3 + j * 3 + 0] = min(255., pow(pixel_color[0], 1 / gamma));
            image.pixels[(image.height - i - 1) * image.width * 3 + j * 3 + 1] = min(255., pow(pixel_color[1], 1 / gamma));
            image.pixels[(image.height - i - 1) * image.width * 3 + j * 3 + 2] = min(255., pow(pixel_color[2], 1 / gamma));
        }
    }
}

void save_image(const Image& image, const int num_path) {
    std::string name_str = "images/";
    name_str.append(to_string(image.width));
    name_str.append("x");
    name_str.append(to_string(image.height));
    name_str.append("_num_path=");
    name_str.append(to_string(num_path));
#ifdef REGULAR
    name_str.append("_regular");
#endif
#ifdef DOF
    name_str.append("_dof");
#endif
#ifdef BLUR
    name_str.append("_blur");
#endif
    name_str.append(".png");

    stbi_write_png(name_str.c_str(), image.width, image.height, 3, &image.pixels[0], 0);
}

int main(int argc, char** argv) {
    // start timer
    auto start_time = system_clock::now();

    // first define the scene, variables, ... 
    printf("Initializing.....\n");

    // image width and height
    int W = IMAGE_WIDTH;
    int H = IMAGE_HEIGHT;

    // camera focal distance
    double f = FOCAL_DISTANCE;

    // number of rays per pixel
    int num_path = NUM_PATH;
    // maximum ray depth (how many bounces allowed)
    int max_ray_depth = MAX_PATH_DEPTH;

    // the image
    Image image;
    image.width = W;
    image.height = H;
    image.pixels = vector<uint8_t>(W * H * 3);

    // the camera
    Camera camera;
    camera.center = Vector(0, 0, 55);
    camera.focus_length = f;

    // the scene
    Scene scene = Scene();
    
    printf("Ray-tracing with %d rays per pixel.....\n", num_path);

    // then scan all pixels
    render(image, camera, scene, num_path, max_ray_depth, 2.2);

    printf("\nDone.\n");

    // save the image 
    save_image(image, num_path);

    // time stuff
    auto end_time = system_clock::now();
    auto duration = duration_cast<microseconds>(end_time - start_time);
    std::cout << "Total rendering time: "
        << double(duration.count()) * microseconds::period::num / microseconds::period::den
        << (double(duration.count()) > 1. ? " seconds" : " second") << std::endl;

    // everything's good
    return 0;
}
