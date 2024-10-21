#ifndef CTENSOR_H
#define CTENSOR_H

#include <complex.h>

typedef enum {
    TRUE,
    FALSE
}Bool;

typedef enum {
    FREE_0,
    FREE_1,
    FREE_2,
    FREE_BOTH,
}FreeFlag;

typedef enum {
    OP_NULL,
    OP_ADD,
    OP_SUBTRACT,
    OP_MULTIPLY,
    OP_MATMUL,
    OP_SIN,
    OP_COS,
    OP_LOG,
    OP_TRANSPOSE,
    OP_SCALAR,
    OP_HERMITIAN,
    OP_ELE_POW
}Operations;

typedef struct {
    double complex* array;
    int     shape[2];
    int     stride[2];
}ComplexMatrix;

typedef struct Tensor {
    Bool           requires_grad;
    double         power;
    double         scalar;
    ComplexMatrix* gradient;
    ComplexMatrix* tensor_matrix;
    Operations     creation_operation;
    struct Tensor* parents[2];
}ComplexTensor;

ComplexMatrix c_empty (int shape[2]);
ComplexMatrix c_zero_matrix (int shape[2]);
ComplexMatrix c_ones_matrix (int shape[2]);
ComplexMatrix c_random_matrix (int shape[2]);
ComplexMatrix c_gaussian_matrix (int shape[2], double mean, double std);

ComplexMatrix* c_sin (ComplexMatrix* z1);
ComplexMatrix* c_cos (ComplexMatrix* z1);
ComplexMatrix* c_log (ComplexMatrix* z1);
ComplexMatrix* c_scalar (ComplexMatrix* z1, double complex* z);
ComplexMatrix* c_add_matrix (ComplexMatrix* z1, ComplexMatrix* z2);
ComplexMatrix* c_sub_matrix (ComplexMatrix* z1, ComplexMatrix* z2);
ComplexMatrix* c_mult_matrix (ComplexMatrix* z1, ComplexMatrix* z2);
ComplexMatrix* c_matmul_matrix (ComplexMatrix* z1, ComplexMatrix* z2);

ComplexTensor Zeros (int shape[2], Bool requires_grad);
ComplexTensor Ones (int shape[2], Bool requires_grad);
ComplexTensor Random (int shape[2], Bool requires_grad);
ComplexTensor Gaussian (int shape[2], double mean, double std, Bool requires_grad);

ComplexTensor* Sin (ComplexTensor* z1);
ComplexTensor* Cos (ComplexTensor* z1);
ComplexTensor* Log (ComplexTensor* z1);
ComplexTensor* Scalar (ComplexTensor* z1, double complex* z);
ComplexTensor* Add (ComplexTensor* z1, ComplexTensor* z2);
ComplexTensor* Sub (ComplexTensor* z1, ComplexTensor* z2);
ComplexTensor* Mult (ComplexTensor* z1, ComplexTensor* z2);
ComplexTensor* MatMul (ComplexTensor* z1, ComplexTensor* z2);

void FFT (ComplexTensor* z);
void iFFT (ComplexTensor* z);
void Backward (ComplexTensor* z);

#endif