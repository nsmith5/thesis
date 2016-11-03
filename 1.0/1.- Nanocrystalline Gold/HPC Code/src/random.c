#include <complex.h>
#include <stdint.h>
#include <math.h>

#define PI 2.0*acos(0)

/* Mersenne Twister Random Number Generator */

/* Mersenne Twister Constants */
enum
{
    N = 624,
    M = 397,
    R = 31,
    A = 0x9908B0DF,

    F = 1812433253,

    U = 11,

    S = 7,
    B = 0x9D2C5680,

    T = 15,
    C = 0xEFC60000,

    L = 18,

    MASK_LOWER = (1ull << R) - 1,
    MASK_UPPER = (1ull << R)
};

static uint32_t  mt[N];
static uint16_t  index;

void seed (const uint32_t  seed)
{
    /* Seed RNG */

    uint32_t  i;

    mt[0] = seed;

    for ( i = 1; i < N; i++ )
    {
        mt[i] = (F * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i);
    }

    index = N;
}

static void twist()
{
    uint32_t  i, x, xA;

    for ( i = 0; i < N; i++ )
    {
        x = (mt[i] & MASK_UPPER) + (mt[(i + 1) % N] & MASK_LOWER);

        xA = x >> 1;

        if ( x & 0x1 )
            xA ^= A;

        mt[i] = mt[(i + M) % N] ^ xA;
    }

    index = 0;
}

// Obtain a 32-bit random number
uint32_t extract_uint32()
{
    uint32_t  y;
    int i = index;

    if ( index >= N )
    {
        twist();
        i = index;
    }

    y = mt[i];
    index = i + 1;

    y ^= (mt[i] >> U);
    y ^= (y << S) & B;
    y ^= (y << T) & C;
    y ^= (y >> L);

    return y;
}

double rand_uniform (void)
{
    /* Generate random number between 0 and 1 */
    return (double) extract_uint32() / (double) 0xffffffff;
}

double complex rand_normal_complex ()
{
    /* Generate a random complex number where the   */
    /* real and imaginary part are normally         */
    /* distributed.                                 */

    double U1, U2;
    double Z1, Z2;

    U1 = rand_uniform ();
    U2 = rand_uniform ();

    Z1 = sqrt(-2.0*log(U1))*cos(2.0*PI*U2);
    Z2 = sqrt(-2.0*log(U1))*sin(2.0*PI*U2);

    return Z1 + I*Z2;
}
