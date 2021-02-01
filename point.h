#pragma once
#include <stdio.h>
#include <gmp.h>


struct point
{
    mpz_t X, Z;
};

void neutral(struct point *point);
void point_init(struct point *point);
void point_add(struct point *q, const struct point *r, const struct point *p1, const mpz_t p);
void point_double(struct point *point, const mpz_t c, const mpz_t p);
void point_free(struct point *point);
