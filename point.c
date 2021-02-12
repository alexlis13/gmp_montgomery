#include <stdio.h>
#include "point.h"
#include "options.h"

void neutral(struct point *point)
{
    mpz_set_str(point->X, "1", 10);
    mpz_set_str(point->Z, "0", 10);
}
/*
* Инициализирует точку и кривую параметрами из стандарта
*/
void point_init(struct point *point)
{
    mpz_inits(point->X, point->Z, 0);
    //
    mpz_t p, v, e, sum, dif;
    mpz_inits(p, v, e, sum, dif, 0);
    mpz_set_str(p, P, 16);
    mpz_set_str(v, V, 16);
    mpz_set_str(e, "1", 10);
    mpz_set_str(point->Z, "1", 10);

    //x = (1 + v) / (1 - v)
    mpz_sub(dif, e, v);
    mpz_add(e, e, v);
    mpz_invert(sum, dif, p);
    mpz_mul(point->X,sum,e);
    mpz_mod(point->X, point->X, p);

    //clear
    mpz_clears(p, v, e, sum, dif, 0);
}

/*
* Складывает точку q с точкой r и результат помещает в q, используя алгоритм xADD
*/
void point_add(struct point *q, const struct point *r, const struct point *p1, const mpz_t p)
{
    mpz_t sq, sr, dq, dr, pow;
    mpz_inits(sq, sr, dq, dr, pow, 0);
    mpz_set_str(pow, "2", 10);

    mpz_sub(dq, q->X, q->Z);
    mpz_sub(dr, r->X, r->Z);
    mpz_add(sq, q->X, q->Z);
    mpz_add(sr, r->X, r->Z);
    mpz_mul(dq, dq,sr);
    mpz_mul(sq, sq, dr);

    //X = (((q.X - q.Z) * (r.X + r.Z) + (q.X + q.Z) * (r.X - r.Z))^2) * Z1
    mpz_add(q->X, dq, sq);
    mpz_powm(q->X, q->X, pow, p);
    mpz_mul(q->X, q->X, p1->Z);
    mpz_mod(q->X, q->X, p);

    //Z = (((q.X - q.Z) * (r.X + r.Z) - (q.X + q.Z) * (r.X - r.Z))^2) * X1
    mpz_sub(q->Z, dq, sq);
    mpz_powm(q->Z, q->Z, pow, p);
    mpz_mul(q->Z, q->Z, p1->X);
    mpz_mod(q->Z, q->Z, p);

    //clear
    mpz_clears(sq, sr, dq, dr, pow, 0);
}

/*
* Удваивает точку point, используя алгоритм xDBL
*/
void point_double(struct point *point, const mpz_t c, const mpz_t p)
{
    mpz_t sp, dp, mul, pow;
    mpz_inits(sp, dp, mul, pow, 0);
    mpz_set_str(pow, "2", 10);

    mpz_add(sp, point->X, point->Z);
    mpz_sub(dp, point->X, point->Z);
    mpz_powm(sp, sp, pow, p);
    mpz_powm(dp, dp, pow, p);

    //X = (point.X + point.Z)^2 * (point.X - point.Z)^2
    mpz_mul(point->X, sp, dp);
    mpz_mod(point->X, point->X, p);

    //Z = ((point.X + point.Z)^2 - (point.X - point.Z)^2) ((A - 2) / 4) * ((point.X + point.Z)^2 - (point.X - point.Z)^2) + (point.X + point.Z)^2)
    //Z = (sp - dp)(C * (sp - dp) + sp)
    mpz_sub(dp, sp, dp);
    mpz_mul(mul, c, dp);
    mpz_add(sp, sp, mul);
    mpz_mul(point->Z, dp, sp);
    mpz_mod(point->Z, point->Z, p);

    //clear
    mpz_clears(sp, dp, mul, pow, 0);
}

void point_free(struct point *point)
{
    mpz_clears(point->X, point->Z, 0);
}
