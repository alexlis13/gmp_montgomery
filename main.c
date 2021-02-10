#include <stdio.h>
#include <gmp.h>
#include "point.h"
#include "curve.h"
#include "options.h"
#include <stdlib.h>

int main()
{
    printf("Initialisation of variables...\n");
    struct point Point;
    point_init(&Point);
    struct curve Curve;
    curve_init(&Curve);
    if(point_on_curve(&Curve, &Point)!=1)
    {
        printf("Point is not on curve!\n");
        exit(1);
    }
    gmp_randstate_t state;
    gmp_randinit_default(state);
    mpz_t pow, res, e, rand, test;
    mpz_inits(pow, res, e, rand, test, 0);
    mpz_set_str(e,"1", 10);
    mpz_set_str(rand, MAX, 10);
    mpz_urandomm(pow, state, rand);
    gmp_printf("Start with Point(%Zd, %Zd),  Power(%Zd)\n", Point.X, Point.Z, pow);

    curve_ladder(&Curve, &Point, pow);
    gmp_printf("Result is: Point(%Zd, %Zd)\n", Point.X, Point.Z);

    printf("Start testing...\n");
    printf("Test 1:\n");
    if(point_on_curve(&Curve, &Point)==1)
        printf("Test 1 true!\n");
    else
        printf("Test 1 false!\n");


    printf("Test 2:\n");
    mpz_set(test, Point.X);
    mpz_add(pow, Curve.q, e);
    curve_ladder(&Curve, &Point, pow);
    if(mpz_cmp(Point.X,test)!=0)
        printf("Test 2 false!\n");
    else
        printf("Test 2 true!\n");
    gmp_printf("Result is: Point(%Zd, %Zd)\n", Point.X, Point.Z);

    printf("Test 3:\n");
    curve_ladder(&Curve, &Point, Curve.q);
    if(mpz_cmp(Point.X, e)!=0)
        printf("Test 3 false!\n");
    else
        printf("Test 3 true!\n");
    gmp_printf("Result is: Point(%Zd, %Zd)\n", Point.X, Point.Z);

    //clear
    printf("Clear all variables...\n");
    mpz_clears(pow, res, e, rand, test, 0);
    gmp_randclear(state);
    curve_clear(&Curve);
    point_free(&Point);
    printf("End!\n");

    return 0;
}

