#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "complex.h"

static size_t total_allocated_memory;

void* malloc_trace (size_t size) {
    void* ptr = malloc(size);
    if (ptr != NULL) {
        total_allocated_memory += size;
    }
    return ptr;
}

void cPrint (Complex* z) {
    if (z->form == coordinate) {
        if (z->data.coord.imag < 0) {
            printf("%f  %fi", z->data.coord.real, z->data.coord.imag);
        } else {
            printf("%f + %fi", z->data.coord.real, z->data.coord.imag);
        }
    } else {
        printf("%f∠%f", z->data.polar.mag, z->data.polar.ang);
    }
}

double sinh (double x) {
    return (exp(x) - exp(-x)) / 2;
}

double cosh (double x) {
    return (exp(x) + exp(-x)) / 2;
}

Complex* Init (double v1, double v2, Form form) {

    Complex* z = (Complex*)malloc_trace(sizeof(Complex));

    switch (form) {
        case (polar) :
            z->form = form;
            z->data.polar.ang = v1;
            z->data.polar.mag = v2;
        
        case (coordinate) :
            z->form = form;
            z->data.coord.real = v1;
            z->data.coord.imag = v2;
    }

    return z;
}

Complex* convert_to_polar (Complex* z) {
    assert(z->form == coordinate);

    double r = sqrt(pow(z->data.coord.real, 2) + pow(z->data.coord.imag, 2));
    double angle = atan2(z->data.coord.imag, z->data.coord.real);
    Complex* newZ = Init(angle, r, polar);

    free(z);

    return newZ;
}

Complex* convert_to_coordinate (Complex* z) {
    assert(z->form == polar);

    double real = z->data.polar.mag * cos(z->data.polar.ang);
    double imag = z->data.polar.mag * sin(z->data.polar.ang);;
    Complex* newZ = Init(real, imag, coordinate);

    free(z);

    return newZ;
}

Complex* cAdd (Complex* z1, Complex* z2) {
    assert(z1->form == coordinate);
    assert(z2->form == coordinate);

    double real = z1->data.coord.real + z2->data.coord.real;
    double imag = z1->data.coord.imag + z2->data.coord.imag;
    Complex* z3 = Init(real, imag, coordinate);

    return z3;
}

Complex* cSub (Complex* z1, Complex* z2) {
    assert(z1->form == coordinate);
    assert(z2->form == coordinate);

    double real = z1->data.coord.real - z2->data.coord.real;
    double imag = z1->data.coord.imag - z2->data.coord.imag;
    Complex* z3 = Init(real, imag, coordinate);

    return z3;
}

Complex* cMult (Complex* z1, Complex* z2) {
    assert(z1->form == z2->form);

    double v1;
    double v2;

    switch (z1->form) {

        case (polar) :
            v1 = z1->data.polar.ang + z2->data.polar.ang;
            v2 = z1->data.polar.mag * z2->data.polar.mag;

        case (coordinate) :
            v1 = (z1->data.coord.real * z2->data.coord.real) - (z1->data.coord.imag * z2->data.coord.imag);
            v2 = (z1->data.coord.real * z2->data.coord.imag) + (z1->data.coord.real * z2->data.coord.imag);
    }

    Complex* z3 = Init(v1, v2, z1->form);

    return z3;
}

Complex* cDiv (Complex* z1, Complex* z2) {
    assert(z1->form == z2->form);

    double v1;
    double v2;

    switch (z1->form) {

        case (polar) :
            assert(z2->data.polar.mag > 0);

            v1 = z1->data.polar.ang - z2->data.polar.ang;
            v2 = z1->data.polar.mag / z2->data.polar.mag;

        case (coordinate) :
            assert (pow(z2->data.coord.real,2) + pow(z2->data.coord.imag,2) > 0);

            v1 = ((z1->data.coord.real * z2->data.coord.real) + 
                  (z1->data.coord.imag * z2->data.coord.imag)) / 
                  (pow(z2->data.coord.real,2) + pow(z2->data.coord.imag,2));

            v2 = ((z1->data.coord.imag * z2->data.coord.real) - 
                  (z1->data.coord.real * z2->data.coord.imag)) / 
                  (pow(z2->data.coord.real,2) + pow(z2->data.coord.imag,2));
    }

    Complex* z3 = Init(v1, v2, z1->form);

    return z3;
}

Complex* cPow (Complex* z, double power) {

    Form og = z->form;

    if (og == coordinate) {
        z = convert_to_polar(z);
    }

    z->data.polar.mag = pow(z->data.polar.mag, power);
    z->data.polar.ang = z->data.polar.ang * power;

    if (og == coordinate && z->form == polar) {
        z = convert_to_coordinate(z);
    }

    return z;
}

Complex* cLog (Complex* z) {
    if (z->form == coordinate) {
        z = convert_to_polar(z);
    }

    double r = log(z->data.polar.mag);
    double a = z->data.polar.ang;

    Complex* nz = Init(r, a, polar);

    return nz;
}

Complex* cSin (Complex* z) {

    if (z->form == polar) {
        z = convert_to_coordinate(z);
    }

    double real = sin(z->data.coord.real) * cosh(z->data.coord.imag);
    double imag = cos(z->data.coord.real) * sinh(z->data.coord.imag);

    Complex* nz = Init(real, imag, coordinate);

    return z;
}

Complex* cCos (Complex* z) {

    if (z->form == polar) {
        z = convert_to_coordinate(z);
    }

    double real = cos(z->data.coord.real) * cosh(z->data.coord.imag);
    double imag = sin(z->data.coord.real) * sinh(z->data.coord.imag);

    Complex* nz = Init(real, imag, coordinate);

    return z;
}

Complex* conjugate (Complex* z) {
    assert (z->form == coordinate);
    
    Complex* z1 = Init(z->data.coord.real, -1 * z->data.coord.imag, coordinate);
    return z1;
}

// gcc complex.c -o exec