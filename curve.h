#include <stdio.h>
#include <gmp.h>
#include "point.h"

struct curve
{
    mpz_t A, B, C, p, q;
};

void curve_init(struct curve *curve);
void curve_ladder(const struct curve *curve, struct point *point, const mpz_t pow);
void curve_clear(struct curve *curve);
int point_on_curve(const struct curve *curve, struct point *point);

