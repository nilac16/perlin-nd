#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#ifdef _OPENMP
#include <omp.h>
#endif //_OPENMP

#include "include/perlin.h"

#define N_OCTAVES 4


void average_char_arrays(size_t N,
                         unsigned char arr1[static restrict N],
                         const unsigned char arr2[static restrict N]);


int main(void)
{
    const size_t density[] = { 16, 16 };
    const size_t resolution[] = { 2048, 2048 };
    const size_t rank = (sizeof resolution) / (sizeof *resolution);
    const unsigned int seed = time(NULL);
    const size_t N = resolution[0] * resolution[1];

    /* unsigned char *img = perlin_noise(rank, density, resolution, seed, 1); */

    unsigned char *img = malloc(N * (sizeof *img));
    if (!img) {
        fputs("Image failed to generate (possibly out of memory)\n", stderr);
        perror("Errno");
        return -1;
    }
    memset(img, 0, N * (sizeof *img));

    unsigned char *oct;
    for (size_t i = N_OCTAVES; i >= 1; i--) {
        oct = perlin_noise(rank, density, resolution, seed, i);
        if (!oct) {
            fputs("Image failed to generate (possibly out of memory)\n", stderr);
            perror("Errno");
            return -1;
        }
        average_char_arrays(N, img, oct);
        free(oct);
    }

    stbi_write_jpg("output.jpg", resolution[0], resolution[1], 1, img, 80);

    /*  This bit is for writing animations of 3D noise

    unsigned char *img = perlin_noise(rank, density, resolution, seed, 1);
    char fname[] = "output/output   .jpg", num[32];
    const size_t dim_offset = resolution[0] * resolution[1];
    for (size_t i = 0; i < resolution[2]; i++) {
        sprintf(num, "%3li", i);
        fname[13] = (i < 100) ? '0' : num[0];
        fname[14] = (i < 10) ? '0' : num[1];
        fname[15] = num[2];
        stbi_write_jpg(fname, resolution[0], resolution[1], 1, img + dim_offset * i, 80);
    } */

    free(img);
    return 0;
}


void average_char_arrays(size_t N,
                         unsigned char arr1[static restrict N],
                         const unsigned char arr2[static restrict N])
{
    for (size_t i = 0; i < N; i++) {
        arr1[i] = arr1[i] / 2 + arr2[i] / 2;
    }
}
