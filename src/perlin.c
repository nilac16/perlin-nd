#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "include/perlin.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __GNUC__
#define gnu_attribute(...) __attribute__((__VA_ARGS__))
#else
#define gnu_attribute(...)
#endif


struct perlin_hypergrid {
    // dimensionality of the image
    size_t rank;

    // number of lattice points in each dimension
    // total number of lattice vectors is the product of this
    // this and the rank are required to initialize a hypergrid
    size_t *restrict density;

    // final resolution of the image
    size_t *restrict resolution;

    // coordinates of each vector
    // indexed with an integer, offset by rank as its stride
    double *restrict locations;

    // vertex vectors, stride is the rank
    double *restrict vectors;
};


gnu_attribute(const)
inline static size_t pow2l(size_t n)
{
    return (size_t)1 << n;
}

gnu_attribute(always_inline)
inline static double random_unity()
{
    return (2 * (double)rand() / RAND_MAX) - 1;
}

gnu_attribute(nonnull)
static void random_unit_vector(size_t dim, double r[static dim])
{
    double inner = 0.;
    for (size_t i = 0; i < dim; i++) {
        r[i] = random_unity();
        inner += r[i] * r[i];
    }
    if (inner > 1.0) {
        random_unit_vector(dim, r);
    } else {
        inner = sqrt(inner);
        for (size_t i = 0; i < dim; i++) {
            r[i] /= inner;
        }
    }
}

gnu_attribute(nonnull, pure)
/// finds the product of every component in @p density
static size_t integrate_density(size_t dim, const size_t density[static dim])
{
    size_t N = 1;
    for (size_t i = 0; i < dim; i++) {
        N *= density[i];
    }
    return N;
}

gnu_attribute(malloc)
struct perlin_hypergrid *new_perlin_hypergrid(size_t rank, const size_t density[static rank])
{
    struct perlin_hypergrid *grid = malloc(sizeof *grid);
    if (!grid) {
        return grid;
    }
    grid->rank = rank;
    grid->density = malloc(rank * (sizeof *grid->resolution));
    if (!grid->density) {
        free(grid);
        return NULL;
    }
    memcpy(grid->density, density, rank * (sizeof *grid->density));
    grid->resolution = malloc(rank * (sizeof *grid->resolution));
    if (!grid->resolution) {
        free(grid->density);
        free(grid);
        return NULL;
    }
    const size_t N = integrate_density(rank, density);
    printf("Hypergrid will have %lu lattice points\n", N);
    printf("Allocating space for %lu doubles\n", N * rank);
    grid->locations = malloc(N * rank * (sizeof *grid->locations));
    if (!grid->locations) {
        free(grid->density);
        free(grid->resolution);
        free(grid);
        return NULL;
    }
    grid->vectors = malloc(N * rank * (sizeof *grid->vectors));
    if (!grid->vectors) {
        free(grid->density);
        free(grid->locations);
        free(grid->resolution);
        free(grid);
        return NULL;
    }
    return grid;
}

gnu_attribute(nonnull)
void free_perlin_hypergrid(struct perlin_hypergrid *hg)
{
    free(hg->vectors);
    free(hg->locations);
    free(hg->density);
    free(hg->resolution);
    free(hg);
}

gnu_attribute(nonnull)
void set_perlin_hypergrid_resolution(struct perlin_hypergrid *hg,
                                     const size_t resolution[static hg->rank])
{
    memcpy(hg->resolution, resolution, hg->rank * (sizeof *hg->resolution));
/// STRATAGEM:
///     Find the 0 coordinate along dimension 0
///     Along dimension 1, copy the entire set of dimension 0 coordinates
///     and set the dimension 1 coordinate
///     Then repeat for dimension 2 with the combined dim 0 and dim 1
///     coordinates
/// You only need to keep track of the current dimensional offset, which
/// increases as the dimension multiplier whenever a loop is completed
    double cell_res, coord;
    size_t d_mult = 1;
    for (size_t dim = 0; dim < hg->rank; dim++) {
                       /* I got fenceposted here vvv */
        cell_res = (double)(hg->resolution[dim] - 1) / (hg->density[dim] - 1);
        coord = 0.;
        size_t stride = d_mult * hg->rank;
        size_t dim_limit = stride * hg->density[dim];
        for (size_t i = 0; i < dim_limit; i += stride) {
            if (i) {
                memcpy(hg->locations + i, hg->locations, stride * (sizeof *hg->locations));
            }
            for (size_t j = 0; j < stride; j += hg->rank) {
                hg->locations[i + j + dim] = coord;
            }
            coord += cell_res;
        }
        d_mult *= hg->density[dim];
    }
}

gnu_attribute(nonnull)
void randomize_perlin_vectors(struct perlin_hypergrid *hg, size_t seed)
{
    srand(seed);
    const size_t N = hg->rank * integrate_density(hg->rank, hg->density);
    for (size_t i = 0; i < N; i += hg->rank) {
        random_unit_vector(hg->rank, hg->vectors + i);
    }
}

gnu_attribute(nonnull, pure)
static size_t find_hypercube_origin(const struct perlin_hypergrid *hg,
                                    const size_t coords[static hg->rank])
{
    // delinearize in-place
    // in each dimension, find out which cell the point is in by dividing the
    // point's coordinate by the pixels per hypercube
    // add that times the dimension multiplier to the running index
    size_t idx = 0, d_mult = 1, hcube_offset;
    for (size_t dim = 0; dim < hg->rank; dim++) {
        hcube_offset = (coords[dim] * (hg->density[dim] - 1)) / hg->resolution[dim];
        idx += d_mult * hcube_offset;
        d_mult *= hg->density[dim];
    }
    return idx;
}

gnu_attribute(nonnull)
// first index contains the hypercube origin
static void locate_hypercube_vertices(const struct perlin_hypergrid *hg,
                                      size_t vertices[static (1 << hg->rank)])
{
    // ok, you figured this out once, you can do it again
    size_t d_mult = 1, n = 1;
    for (size_t dim = 0; dim < hg->rank; dim++) {
        memcpy(vertices + n, vertices, n * (sizeof *vertices));
        for (size_t i = 0; i < n; i++) {
            vertices[n + i] += d_mult;
        }
        d_mult *= hg->density[dim];
        n *= 2;
    }
}

gnu_attribute(nonnull, pure)
static const double *get_hypergrid_lattice_vector(const struct perlin_hypergrid *hg, size_t i)
{
    return hg->vectors + i * hg->rank;
}

gnu_attribute(nonnull, pure)
static const double *get_hypergrid_lattice_point(const struct perlin_hypergrid *hg, size_t i)
{
    return hg->locations + i * hg->rank;
}

gnu_attribute(nonnull, pure)
/// dot the input coordinates with the vector at @p i in @p hg
static double dot_perlin_vector(const struct perlin_hypergrid *hg, size_t h,
                                const size_t coords[static hg->rank])
{
    const double *const vect = get_hypergrid_lattice_vector(hg, h);
    const double *const origin = get_hypergrid_lattice_point(hg, h);
    double x = 0.;
    for (size_t i = 0; i < hg->rank; i++) {
        x += vect[i] * ((double)coords[i] - origin[i]);
    }
    //#endif
    return x;
}

gnu_attribute(const, unused)
/// Input shall be the dot products along some dimension as @p p0 and @p p1, 
/// and the RELATIVE coordinate in that dimension of the field point @p x
///
/// \param p0
///     Dot product at lower bound
/// \param p1
///     Dot product at upper bound
/// \param x
///     Normalized distance between the endpoints
/// \returns Interpolated dot product along this dimension
///
static double interpol8(double p0, double p1, double x)
{
    return p0 + (p1 - p0) * x * x * x * (10 - 15 * x + 6 * x * x);
}

gnu_attribute(nonnull, pure)
static double dot_perlin_hypercube(const struct perlin_hypergrid *hg, 
                                   const size_t coords[static restrict hg->rank],
                                   size_t vertices[static restrict (1 << hg->rank)],
                                   double dots[static (1 << hg->rank)])
{
    *vertices = find_hypercube_origin(hg, coords);
    // now dot with hyper's vector, and all of its cohypercubical vertices
    locate_hypercube_vertices(hg, vertices);
    const size_t N = pow2l(hg->rank);
    for (size_t i = 0; i < N; i++) {
        dots[i] = dot_perlin_vector(hg, vertices[i], coords);
    }
    // now "interpolate" along each dimension
    for (size_t dim = 0; dim < hg->rank; dim++) {
        size_t di = pow2l(dim);
        double relative = (double)coords[dim] - hg->locations[hg->rank * vertices[0] + dim];
        relative /= hg->locations[hg->rank * vertices[di] + dim] - hg->locations[hg->rank * vertices[0] + dim];
        for (size_t i = 0; i < N; i += 2 * di) {
            dots[i] = interpol8(dots[i], dots[i + di], relative);
        }
    }
    return *dots;
}

gnu_attribute(nonnull)
/// This function is a bottleneck
static void delinearize_coordinate(int rank, size_t coords[static restrict rank],
                                   const size_t resolution[static restrict rank])
{
    size_t idx = coords[0];
    coords[0] = 1;
    for (int i = 1; i < rank; i++) {
        // put the dimension multiple in coords[i]
        coords[i] = coords[i - 1] * resolution[i - 1];
    }
    size_t d_mult;
    for (int i = rank - 1; i >= 0; i--) {
        d_mult = coords[i];
        coords[i] = idx / d_mult;
        idx -= coords[i] * d_mult;
    }
}

gnu_attribute(nonnull, malloc)
float *generate_perlin_noise(const struct perlin_hypergrid *hg, size_t *n_pix)
{
    const size_t N = *n_pix = integrate_density(hg->rank, hg->resolution);
    float *img = malloc(N * (sizeof *img));
    if (!img) {
        *n_pix = 0;
        return img;
    }
    size_t *local_coords, *local_vertices;
    size_t N_locals = hg->rank, N_vertices = pow2l(hg->rank);
    size_t thread_coord_offset = 0, thread_vertex_offset = 0;
    double *local_dots;
    #ifdef _OPENMP
    N_locals *= omp_get_max_threads();
    N_vertices *= omp_get_max_threads();
    printf("Allocating space for %lu possible threads, %lu dimensions: %lu\n", (size_t)omp_get_max_threads(), hg->rank, N_locals);
    #endif //_OPENMP
    local_coords = malloc(N_locals * (sizeof *local_coords));
    if (!local_coords) {
        *n_pix = 0;
        free(img);
        return NULL;
    }
    local_vertices = malloc(N_vertices * (sizeof *local_vertices));
    if (!local_vertices) {
        *n_pix = 0;
        free(local_coords);
        free(img);
        return NULL;
    }
    local_dots = malloc(N_vertices * (sizeof *local_dots));
    if (!local_dots) {
        *n_pix = 0;
        free(local_vertices);
        free(local_coords);
        free(img);
        return NULL;
    }
    #ifdef _OPENMP
    #pragma omp parallel private(thread_coord_offset, thread_vertex_offset)
    {
    thread_coord_offset = hg->rank * omp_get_thread_num();
    thread_vertex_offset = pow2l(hg->rank) * omp_get_thread_num();
    #pragma omp for
    #endif //_OPENMP
    for (size_t i = 0; i < N; i++) {
        local_coords[thread_coord_offset] = i;
        delinearize_coordinate(hg->rank, local_coords + thread_coord_offset,
                                         hg->resolution);
        img[i] = dot_perlin_hypercube(hg, local_coords + thread_coord_offset,
                                          local_vertices + thread_vertex_offset,
                                          local_dots + thread_vertex_offset);
    }
    #ifdef _OPENMP
    }
    #endif //_OPENMP
    free(local_dots);
    free(local_vertices);
    free(local_coords);
    return img;
}

gnu_attribute(nonnull)
static void floating_extrema(size_t N,
                             const float arr[static restrict N],
                             float *restrict min, float *restrict max)
{
    if (!N) {
        return;
    }
    *min = *max = *arr;
    for (size_t i = 1; i < N; i++) {
        if (arr[i] > *max) {
            *max = arr[i];
        }
        if (arr[i] < *min) {
            *min = arr[i];
        }
    }
}

gnu_attribute(nonnull, malloc)
static unsigned char *convert_noise_to_image(size_t N, const float img[static N])
{
    unsigned char *out = malloc(N * (sizeof *out));
    if (!out) {
        return out;
    }
    float min = 0., max = 0., value;
    floating_extrema(N, img, &min, &max);
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif //_OPENMP
    for (size_t i = 0; i < N; i++) {
        //value = 255. * (img[i] - min) / (max - min);
        value = interpol8(0., 255., (img[i] - min) / (max - min));
        out[i] = value;
    }
    return out;
}

gnu_attribute(nonnull, malloc)
unsigned char *perlin_noise(const size_t dimensions,
                            const size_t density[static dimensions],
                            const size_t resolution[static dimensions],
                            size_t seed, int octave)
{
    size_t density_actual[dimensions];
    for (size_t i = 0; i < dimensions; i++) {
        density_actual[i] = density[i] * ((size_t)1 << octave);
    }
    struct perlin_hypergrid *g = new_perlin_hypergrid(dimensions, density_actual);
    if (!g) {
        return NULL;
    }
    set_perlin_hypergrid_resolution(g, resolution);
    randomize_perlin_vectors(g, seed);
    size_t N;
    float *noise = generate_perlin_noise(g, &N);
    free_perlin_hypergrid(g);
    unsigned char *img = convert_noise_to_image(N, noise);
    free(noise);
    return img;
}

/* gnu_attribute(nonnull)
static void linearize_coordinates(int rank, size_t coords[static restrict rank],
                                  const size_t resolution[static restrict rank])
{
    size_t d_mult = 1, idx = 0;
    for (int dim = 0; ; dim++) {
        idx += d_mult * coords[dim];
        if (dim >= rank - 1) {
            // avoid that extra multiplication
            break;
        }
        d_mult *= resolution[dim];
    }
    coords[0] = idx;
} */
