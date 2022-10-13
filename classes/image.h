#ifndef IMAGE
#define IMAGE

#include "vector.h"
#include <vector>

using namespace std;

struct Image {
    int width;
    int height;
    vector<uint8_t> pixels;
};

#endif