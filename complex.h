#ifndef COMPLEX_H
#define COMPLEX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef enum {
    polar,
    coordinate
}Form;

typedef struct {
    double real;
    double imag;
}Coordinate;

typedef struct {
    double mag;
    double ang;
}Polar;

typedef struct {
    Form form;
    union data {
        Polar polar;
        Coordinate coord;
    }data;
}Complex;

void* malloc_trace (size_t size);
void cPrint (Complex* z);

Complex* Init (double v1, double v2, Form form);
Complex* convert_to_polar (Complex* z);
Complex* convert_to_coordinate (Complex* z);
Complex* cLog (Complex* z1);
Complex* cSin (Complex* z1);
Complex* cCos (Complex* z1);
Complex* cPow (Complex* z1, double n);
Complex* cAdd (Complex* z1, Complex* z2);
Complex* cSub (Complex* z1, Complex* z2);
Complex* cMult (Complex* z1, Complex* z2);
Complex* cDiv (Complex* z1, Complex* z2);
Complex* cScale (Complex* z1, double scale);
Complex* conjugate (Complex* z);

#endif