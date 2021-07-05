#pragma once

#ifndef PERLIN_H
#define PERLIN_H


/// Generates an 8-bit image of Perlin noise given the passed parameters
///
/// \param dimensions
///     The dimensionality of the image
/// \param density
///     The number of hypergrid lattice points in each dimension. Must be of size @p dimensions
/// \param resolution
///     The resolution of the image in each dimension. Must also be of size @p dimensions
/// \param seed
///     The seed for the libstdc number generator used to create the lattice vectors
/// \param octave
///     Constant multiple applied to @p density. Useful for creating fractal noise
/// \returns An 8-bit amplitude image of Perlin noise, or NULL on failure
///
unsigned char *perlin_noise(size_t dimensions,
                            const size_t density[static dimensions],
                            const size_t resolution[static dimensions],
                            size_t seed, int octave);

/**
 *  Everything below this point can be used to generate a raw image of floats. */

struct perlin_hypergrid;


/// Creates a new Perlin hypergrid with the minimum necessary values to
/// make it unique
///
/// \param dimensions
///     The dimensionality of the hypergrid, and the size of the density
///     and resolution arrays passed to its initializers
/// \param density
///     The number of lattice points in each dimension
/// \returns A new hypergrid object used for generating the noise, or NULL
///     on allocation failure. Its resolution and vectors are uninitialized!
///
struct perlin_hypergrid *
new_perlin_hypergrid(size_t dimensions, const size_t density[static dimensions]);


void free_perlin_hypergrid(struct perlin_hypergrid *grid);


/// Sets the resolution of the image that will be produced by the passed
/// hypergrid
///
/// \param grid
///     Hypergrid to be changed
/// \param resolution
///     Resolution of the resulting image. Must be of size grid->rank
///
void set_perlin_hypergrid_resolution(struct perlin_hypergrid *grid,
                                     const size_t resolution[]);


/// Randomizes each vector using the Monte Carlo scheme/rejection sampling
///
/// \param grid
///     Hypergrid to be changed
/// \param seed
///     Seed for the random number generator. The libstdc generator is used
///
void randomize_perlin_vectors(struct perlin_hypergrid *grid, size_t seed);


/// Generates an image of floats from the passed, FULLY INITIALIZED hypergrid
///
/// \param grid
///     Hypergrid to use for image generation
/// \param N
///     Size of the returned image array. This is set to zero if allocation fails
/// \returns An array of floats representing the raw image, or NULL on failure
///
float *generate_perlin_noise(const struct perlin_hypergrid *grid, size_t *N);


#endif //PERLIN_H
