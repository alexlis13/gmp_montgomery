#include <stdio.h>
#include <gmp.h>
#include "curve.h"
#include "options.h"
#include "point.h"

void curve_init(struct curve *curve)
{
    mpz_inits(curve->A,curve->B, curve->C, curve->p, curve->q, 0);
    //
    mpz_t e, a, d, sum, dif;
    mpz_inits(e, a, d, sum, dif, 0);
    mpz_set_str(curve->p, P, 16);
    mpz_set_str(curve->q, Q, 16);
    mpz_set_str(e, E, 16);
    mpz_set_str(d, D, 16);

    // A = 2 * (e + d) / (e - d)
    mpz_add(sum, e, d);
    mpz_sub(dif, e, d);
    mpz_set_str(curve->A, "2", 10);
    mpz_mul(curve->A, curve->A, sum);
    mpz_invert(dif,dif,curve->p);
    mpz_mul(curve->A, curve->A, dif);
    mpz_mod(curve->A, curve->A,curve->p);

    // B = 4 / (e - d)
    mpz_set_str(curve->B, "4", 10);
    mpz_mul(curve->B,curve->B,dif);
    mpz_mod(curve->B, curve->B, curve->p);

    // C = (A - 2) / 4
    mpz_set_str(curve->C, "2", 10);
    mpz_sub (curve->C, curve->A, curve->C);
    mpz_set_str(a, "4", 10);
    mpz_invert(a,a,curve->p);
    mpz_mul(curve->C, curve->C, a);
    mpz_mod(curve->C, curve->C, curve->p);

    //clear
    mpz_clears(a, d, e, sum, dif, 0);
}
/*
* 	Выполняет возведение точки point в степень k на кривой по алгоритму «Лесенка Монтгомери»
*/

void curve_ladder(const struct curve *curve, struct point *point, const mpz_t pow)
{
     size_t size = mpz_sizeinbase (pow, 2); //Size of pow in bites

     struct point r, q;
     point_init(&r);
     point_init(&q);
     neutral(&q);
     mpz_set(r.X, point->X);
     mpz_set(r.Z, point->Z);

     for(int i = size - 1; i>=0;i--)
     {
        if(mpz_tstbit(pow,i))
        {
            point_add(&q, &r, point, curve->p);
            point_double(&r, curve->C, curve->p);
        }
        else
        {
            point_add(&r, &q, point, curve->p);
            point_double(&q, curve->C, curve->p);
        }
     }

     mpz_set(point->X, q.X);
     mpz_set(point->Z, q.Z);

     //x = X * Z^(-1) mod p
     mpz_t i;
     mpz_init(i);
     if (mpz_sgn(point->Z)!=0)
     {
         mpz_invert(i, point->Z, curve->p);
         mpz_mul(point->X, point->X, i);
         mpz_mod(point->X, point->X, curve->p);
         mpz_mul(point->Z, point->Z, i);
         mpz_mod(point->Z, point->Z, curve->p);
     }
     else
     {
         mpz_set_str(point->X, "1", 10);
     }

     //clear
     mpz_clear(i);
     point_free(&q);
     point_free(&r);
}

void curve_clear(struct curve *curve)
{
    mpz_clears(curve->A, curve->B, curve->C, curve->p, curve->q, 0);
}
/*
* 	Проверяет, находится ли точка point на кривой
*/
int point_on_curve(const struct curve *curve, struct point *point)
{
    mpz_t i, n, pow, mul;
    mpz_inits(i, n, pow, mul, 0);
    mpz_set_str(pow, "2", 10);

    // Y^2 = n = ((X^2 + A*X)X + X)/B
    mpz_powm(n, point->X, pow, curve->p);
    mpz_mul(mul, curve->A, point->X);
    mpz_add(n, n ,mul);
    mpz_mul(n, n, point->X);
    mpz_add(n, n, point->X);
    mpz_invert(i, curve->B, curve->p);
    mpz_mul(n, n, i);
    return mpz_jacobi(n, curve->p);
}


