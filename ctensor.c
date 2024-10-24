#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "complex.h"
#include "ctensor.h"

#define PI 3.14159265358979323846

// int is_power_of_two(int n) {
//     return (n > 0) && ((n & (n - 1)) == 0);
// }

// void free_complexmatrix (ComplexMatrix* matrix) {
//     if (matrix != NULL) {
//         free(matrix->array);
//         free(matrix);
//     }
// }

// void cPrintMatrix (ComplexMatrix* z) {
//     int rows = z->shape[0];
//     int cols = z->shape[1];

//     for (int i = 0; i < rows; i++) {
//         for (int j = 0; j < cols; j++) {
//             cPrint(&z->array[(i * z->stride[0]) + (j * z->stride[1])]);
//         }
//         printf("\n");
//     }
//     printf("\n");
// }

// void cPrintTensor (ComplexTensor* z, int grad) {
//     if (grad) {
//         cPrintMatrix(z->gradient);
//     } else {
//         cPrintMatrix(z->tensor_matrix);
//     }
// }

// double normal_random (double mean, double std) {
//     double u1, u2, z1;
//     do {
//         u1 = (double)rand() / RAND_MAX;
//         u2 = (double)rand() / RAND_MAX;
//     }while (u1 <= 0 || u2 <= 0); 

//     z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

//     return mean + std * z1;
// }

// Complex* fft_1d(Complex* array, int size, int inverse) {

//     if (size == 1) {
//         Complex* result = malloc(sizeof(Complex));
//         *result = *array;
//         return result;
//     }

//     if (!is_power_of_two(size)) {
//         return NULL;
//     }

//     Complex* even = (Complex*)malloc(sizeof(Complex) * (size / 2));
//     Complex* odd = (Complex*)malloc(sizeof(Complex) * (size / 2));

//     for (int i = 0; i < size / 2; i++) {
//         even[i] = array[2 * i];
//         odd[i] = array[(2 * i) + 1];
//     }
    
//     Complex* even_fft = fft_1d(even, size / 2, inverse);
//     Complex* odd_fft = fft_1d(odd, size / 2, inverse);

//     free(even);
//     free(odd);

//     Complex* result = malloc(sizeof(Complex) * size);
//     for (int k = 0; k < size / 2; k++) {

//         double angle = (inverse ? 2 : -2) * PI * k / size;
//         Complex* twiddle_factor = Init(cos(angle), sin(angle), coordinate);

//         Complex* product = cMult(twiddle_factor, &odd_fft[k]);
//         result[k] = *cAdd(&even_fft[k], product);
//         result[k + size / 2] = *cSub(&even_fft[k], product);

//         if (inverse) {
//             result[k] = *cScale(&result[k], 0.5);
//             result[k + size / 2] = *cScale(&result[k + size / 2], 0.5);
//         }
//     }

//     free(even_fft);
//     free(odd_fft);

//     return result;
// }

// void fft_2d(ComplexMatrix* z, int inverse) {

//     int rows = z->shape[0];
//     int cols = z->shape[1];

//     if (!is_power_of_two(rows) || !is_power_of_two(cols)) {
//         fprintf(stderr, "Both dimensions must be powers of 2\n");
//         return;
//     }

//     for (int i = 0; i < rows; i++) {
//         Complex* row_fft = fft_1d(&z->array[i * z->stride[0]], cols, inverse);
//         for (int j = 0; j < cols; j++) {
//             z->array[(i * z->stride[0]) + (j * z->stride[1])] = row_fft[j];
//         }
//         free(row_fft);
//     }

//     Complex* col = malloc(sizeof(Complex) * rows);
//     for (int j = 0; j < cols; j++) {
//         for (int i = 0; i < rows; i++) {
//             col[i] = z->array[(i * z->stride[0]) + (j * z->stride[1])];
//         }
//         Complex* col_fft = fft_1d(col, rows, inverse);
//         for (int i = 0; i < rows; i++) {
//             z->array[(i * z->stride[0]) + (j * z->stride[1])] = col_fft[i];
//         }
//         free(col_fft);
//     }
//     free(col);

//     if (inverse) {
//         double scale = 1.0 / (rows * cols);
//         for (int i = 0; i < rows; i++) {
//             for (int j = 0; j < cols; j++) {
//                 z->array[(i * z->stride[0]) + (j * z->stride[1])] = *cScale(
//                     &z->array[(i * z->stride[0]) + (j * z->stride[1])], 
//                     scale
//                 );
//             }
//         }
//     }
// }

// ComplexMatrix* c_empty (int shape[2]) {
//     ComplexMatrix* new_matrix = (ComplexMatrix*)malloc_trace(sizeof(ComplexMatrix));
//     new_matrix->array = (Complex*)malloc_trace((size_t)(shape[0] * shape[1]) * sizeof(Complex));
//     new_matrix->shape[0] = shape[0];
//     new_matrix->shape[1] = shape[1];
//     new_matrix->stride[0] = shape[1];
//     new_matrix->stride[1] = 1;

//     return new_matrix;
// }

// ComplexMatrix* c_zero_matrix (int shape[2]) {
//     ComplexMatrix* z = c_empty(shape);
//     for (int i = 0; i < shape[0]; i++) {
//         for (int j = 0; j < shape[1]; j++) {
//             z->array[
//                 i * z->stride[0] +
//                 j * z->stride[1]
//             ] = *Init(0, 0, coordinate);
//         }
//     }

//     return z;
// }

// ComplexMatrix* c_ones_matrix (int shape[2]) {
//     ComplexMatrix* z = c_empty(shape);
//     for (int i = 0; i < shape[0]; i++) {
//         for (int j = 0; j < shape[1]; j++) {
//             z->array[
//                 i * z->stride[0] +
//                 j * z->stride[1]
//             ] = *Init(1, 0, coordinate);
//         }
//     }

//     return z;
// }

// ComplexMatrix* c_random_matrix (int shape[2]) {
//     ComplexMatrix* z = c_empty(shape);
//     for (int i = 0; i < shape[0]; i++) {
//         for (int j = 0; j < shape[1]; j++) {
//             double r1 = rand();
//             double r2 = rand();
//             z->array[
//                 i * z->stride[0] +
//                 j * z->stride[1]
//             ] = *Init((double)r1 / (double)RAND_MAX, 
//                       (double)r2 / (double)RAND_MAX, 
//                       coordinate
//                     );
//         }
//     }

//     return z;
// }

// ComplexMatrix* c_gaussian_matrix (int shape[2], double mean, double std) {
//     srand(time(NULL));
//     ComplexMatrix* z = c_empty(shape);
//     for (int i = 0; i < shape[0]; i++) {
//         for (int j = 0; j < shape[1]; j++) {
//             double r1 = rand();
//             double r2 = rand();
//             z->array[
//                 i * z->stride[0] +
//                 j * z->stride[1]
//             ] = *Init(normal_random(mean, std), 
//                       normal_random(mean, std), 
//                       coordinate
//                     );
//         }
//     }

//     return z;
// }

// ComplexMatrix* c_sin (ComplexMatrix* z1, FreeFlag free) {
//     ComplexMatrix* m3 = c_empty(z1->shape);
//     for (int i = 0; i < z1->shape[0]; i++) {
//         for (int j = 0; j < z1->shape[1]; j++) {
//             m3->array[
//                 (i * m3->stride[0]) + 
//                 (j * m3->stride[1])
//             ] = *cSin(&z1->array[(i * z1->stride[0]) + (j * z1->stride[1])]);
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(z1);
//             break;
        
//         default:
//             break;
            
//     }

//     return m3;
// }

// ComplexMatrix* c_cos (ComplexMatrix* z1, FreeFlag free) {
//     ComplexMatrix* m3 = c_empty(z1->shape);
//     for (int i = 0; i < z1->shape[0]; i++) {
//         for (int j = 0; j < z1->shape[1]; j++) {
//             m3->array[
//                 (i * m3->stride[0]) + 
//                 (j * m3->stride[1])
//             ] = *cCos(&z1->array[(i * z1->stride[0]) + (j * z1->stride[1])]);
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(z1);
//             break;
        
//         default:
//             break;
            
//     }

//     return m3;
// }

// ComplexMatrix* c_log (ComplexMatrix* z1, FreeFlag free) {
//     ComplexMatrix* m3 = c_empty(z1->shape);
//     for (int i = 0; i < z1->shape[0]; i++) {
//         for (int j = 0; j < z1->shape[1]; j++) {
//             m3->array[
//                 (i * m3->stride[0]) + 
//                 (j * m3->stride[1])
//             ] = *cLog(&z1->array[(i * z1->stride[0]) + (j * z1->stride[1])]);
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(z1);
//             break;
        
//         default:
//             break;
            
//     }

//     return m3;
// }

// ComplexMatrix* c_negate (ComplexMatrix* z1, FreeFlag free) {
//     ComplexMatrix* m3 = c_empty(z1->shape);
//     for (int i = 0; i < z1->shape[0]; i++) {
//         for (int j = 0; j < z1->shape[1]; j++) {
//             m3->array[
//                 (i * m3->stride[0]) + 
//                 (j * m3->stride[1])
//             ] = *negate(&z1->array[(i * z1->stride[0]) + (j * z1->stride[1])]);
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(z1);
//             break;
        
//         default:
//             break;
            
//     }

//     return m3;
// }

// ComplexMatrix* c_add_matrix (ComplexMatrix* m1, ComplexMatrix* m2, FreeFlag free) {

//     ComplexMatrix* m3 = c_empty(m1->shape);
//     for (int i = 0; i < m1->shape[0]; i++) {
//         for (int j = 0; j < m2->shape[1]; j++) {
//             m3->array[
//                 i * m3->stride[0] +
//                 j * m3->stride[1]
//             ] = *cAdd(&m1->array[
//                 i * m1->stride[0] +
//                 j * m1->stride[1]
//             ], &m2->array[
//                 i * m2->stride[0] +
//                 j * m2->stride[1]
//             ]);
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(m1);
//             break;
        
//         case FREE_2 :
//             free_complexmatrix(m2);
//             break;
        
//         case FREE_BOTH :
//             free_complexmatrix(m1);
//             free_complexmatrix(m2);
//             break;
        
//         default:
//             break;

//     }

//     return m3;
// }

// ComplexMatrix* c_sub_matrix (ComplexMatrix* m1, ComplexMatrix* m2, FreeFlag free) {

//     ComplexMatrix* m3 = c_empty(m1->shape);
//     for (int i = 0; i < m1->shape[0]; i++) {
//         for (int j = 0; j < m2->shape[1]; j++) {
//             m3->array[
//                 i * m3->stride[0] +
//                 j * m3->stride[1]
//             ] = *cSub(&m1->array[
//                 i * m1->stride[0] +
//                 j * m1->stride[1]
//             ], &m2->array[
//                 i * m2->stride[0] +
//                 j * m2->stride[1]
//             ]);
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(m1);
//             break;
        
//         case FREE_2 :
//             free_complexmatrix(m2);
//             break;
        
//         case FREE_BOTH :
//             free_complexmatrix(m1);
//             free_complexmatrix(m2);
//             break;
        
//         default:
//             break;

//     }

//     return m3;
// }

// ComplexMatrix* c_mult_matrix (ComplexMatrix* m1, ComplexMatrix* m2, FreeFlag free) {

//     ComplexMatrix* m3;

//     if (m1 && m2) {
//         m3 = c_empty(m1->shape);
//     } else {
//         exit(EXIT_FAILURE);
//     }
    
//     if (m1 && m2 && m3) {
//         for (int i = 0; i < m1->shape[0]; i++) {
//             for (int j = 0; j < m2->shape[1]; j++) {
//                 m3->array[
//                     i * m3->stride[0] +
//                     j * m3->stride[1]
//                 ] = *cMult(&m1->array[
//                     i * m1->stride[0] +
//                     j * m1->stride[1]
//                 ], &m2->array[
//                     i * m2->stride[0] +
//                     j * m2->stride[1]
//                 ]);
//             }
//         }
//     } else {
//         exit(EXIT_FAILURE);
//     }

//     if (m1 && m2) {
//         switch (free) {
        
//             case FREE_1 :
//                 free_complexmatrix(m1);
//                 break;
        
//             case FREE_2 :
//                 free_complexmatrix(m2);
//                 break;
        
//             case FREE_BOTH :
//                 free_complexmatrix(m1);
//                 free_complexmatrix(m2);
//                 break;
        
//             default:
//                 break;

//         }
//     } else {
//         exit(EXIT_FAILURE);
//     }

//     return m3;
// }

// ComplexMatrix* c_matmul_matrix (ComplexMatrix* z1, ComplexMatrix* z2, FreeFlag free) {

//     int new_shape[2] = {z1->shape[0], z2->shape[1]};
//     ComplexMatrix* m3 = c_empty(new_shape);

//     for (int i = 0; i < m3->shape[0]; i++) {
//         for (int j=0; j < m3->shape[1]; j++) {
//             m3->array[i * m3->stride[0] + j * m3->stride[1]] = *Init(0, 0, coordinate);
//             for (int k = 0; k < z1->shape[1]; k++) {
//                 m3->array[i * m3->stride[0] + j * m3->stride[1]] = *cAdd( 
//                     &m3->array[i * m3->stride[0] + j * m3->stride[1]], 
//                     cMult(&z1->array[i * z1->stride[0] + k * z1->stride[1]],
//                           &z2->array[k * z2->stride[0] + j * z2->stride[1]]
//                         )
//                     );
//             }
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(z1);
//             break;
        
//         case FREE_2 :
//             free_complexmatrix(z2);
//             break;
        
//         case FREE_BOTH :
//             free_complexmatrix(z1);
//             free_complexmatrix(z2);
//             break;
        
//         default:
//             break;
            
//     }

//     return m3;
// }

// ComplexMatrix* c_scalar (ComplexMatrix* m1, Complex* value, FreeFlag free) {

//     ComplexMatrix* new_matrix = c_empty(m1->shape);

//     for (int i = 0; i < m1->shape[0]; i++) {
//         for (int j = 0; j < m1->shape[1]; j++) {
//             new_matrix->array[
//                 i * new_matrix->stride[0] +
//                 j * new_matrix->stride[1]
//             ] = *cMult(&m1->array[i * m1->stride[0] + j * m1->stride[1]], value); 
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(m1);
//             break;
        
//         default:
//             break;
            
//     }

//     return new_matrix;
// }

// ComplexMatrix* c_transpose_matrix (ComplexMatrix* m1, FreeFlag free) {
    
//     int shape[2] = {m1->shape[1], m1->shape[0]};

//     ComplexMatrix* new_matrix = c_empty(shape);

//     for (int i = 0; i < new_matrix->shape[0]; i++) {
//         for (int j = 0; j < new_matrix->shape[1]; j++) {
//             new_matrix->array[
//                 (i * new_matrix->stride[0]) + 
//                 (j * new_matrix->stride[1])
//             ] = m1->array[(i * m1->stride[1]) + (j * m1->stride[0])];
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(m1);
//             break;
        
//         default:
//             break;
            
//     }

//     return new_matrix;
// };

// ComplexMatrix* c_hermitian (ComplexMatrix* m1, FreeFlag free) {

//     int shape[2] = {m1->shape[1], m1->shape[0]};

//     ComplexMatrix* new_matrix = c_empty(shape);

//     for (int i = 0; i < new_matrix->shape[0]; i++) {
//         for (int j = 0; j < new_matrix->shape[1]; j++) {
//             new_matrix->array[
//                 (i * new_matrix->stride[0]) + 
//                 (j * new_matrix->stride[1])
//             ] = *conjugate(&m1->array[(i * m1->stride[1]) + (j * m1->stride[0])]);
//         }
//     }

//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(m1);
//             break;
        
//         default:
//             break;
            
//     }

//     return new_matrix;
// }

// ComplexMatrix* c_ele_power (ComplexMatrix* m1, double n, FreeFlag free) {

//     ComplexMatrix* m2 = c_empty(m1->shape);
//     for (int i = 0; i < m1->shape[0]; i++) {
//         for (int j = 0; j < m1->shape[1]; j++) {
//             m2->array[
//                 (i * m2->stride[0]) + 
//                 (j * m2->stride[1])
//             ] = *cPow(&m1->array[(i * m1->stride[0]) + (j * m1->stride[1])], n);
//         }
//     }
//     switch (free) {
        
//         case FREE_1 :
//             free_complexmatrix(m1);
//             break;
        
//         default:
//             break;
            
//     }
//     return m2;
// }

// ComplexTensor* createTensor (int shape[2]) {
//     assert (shape[0] > 0 && shape[1] > 0);

//     ComplexTensor* matrix = (ComplexTensor*)malloc_trace(sizeof(ComplexTensor));
//     matrix->gradient = c_zero_matrix(shape);
//     matrix->tensor_matrix = c_empty(shape);

//     matrix->power = 1;
//     matrix->scalar = 1;
//     matrix->parents[0] = NULL;
//     matrix->parents[1] = NULL;
//     matrix->requires_grad = FALSE;
//     matrix->creation_operation = OP_NULL;

//     return matrix;
// }

// ComplexTensor* Zeros (int shape[2], Bool requires_grad) {

//     ComplexTensor* zeros = createTensor(shape);
//     free_complexmatrix(zeros->tensor_matrix);
//     zeros->tensor_matrix = c_zero_matrix(shape);
//     zeros->requires_grad = requires_grad;

//     return zeros;
// }

// ComplexTensor* Ones (int shape[2], Bool requires_grad) {

//     ComplexTensor* ones = createTensor(shape);
//     free_complexmatrix(ones->tensor_matrix);
//     ones->tensor_matrix = c_ones_matrix(shape);
//     ones->requires_grad = requires_grad;

//     return ones;
// }

// ComplexTensor* Random (int shape[2], Bool requires_grad) {

//     ComplexTensor* ones = createTensor(shape);
//     ones->tensor_matrix = c_random_matrix(shape);
//     ones->requires_grad = requires_grad;

//     return ones;
// }

// ComplexTensor* Gaussian (int shape[2], double mean, double std, Bool requires_grad) {
    
//     srand(time(NULL));
//     ComplexTensor* normal = createTensor(shape);
//     normal->tensor_matrix = c_gaussian_matrix(shape, mean, std);
//     normal->requires_grad = requires_grad;

//     return normal;
// }

// ComplexTensor* Sin (ComplexTensor* z1) {

//     ComplexTensor* new = createTensor(z1->tensor_matrix->shape);

//     new->creation_operation = OP_SIN;
//     new->parents[0] = z1;
    
//     free_complexmatrix(new->tensor_matrix);
//     new->tensor_matrix = c_sin(z1->tensor_matrix, FREE_0);

//     if (z1->requires_grad == TRUE) {
//         new->requires_grad = TRUE;
//     } else {
//         new->requires_grad = FALSE;
//     }

//     return new;
// };

// ComplexTensor* Cos (ComplexTensor* z1) {

//     ComplexTensor* new = createTensor(z1->tensor_matrix->shape);

//     new->creation_operation = OP_COS;
//     new->parents[0] = z1;
    
//     free_complexmatrix(new->tensor_matrix);
//     new->tensor_matrix = c_cos(z1->tensor_matrix, FREE_0);

//     if (z1->requires_grad == TRUE) {
//         new->requires_grad = TRUE;
//     } else {
//         new->requires_grad = FALSE;
//     }

//     return new;
// };

// ComplexTensor* Log (ComplexTensor* z1) {

//     ComplexTensor* new = createTensor(z1->tensor_matrix->shape);

//     new->creation_operation = OP_LOG;
//     new->parents[0] = z1;
    
//     free_complexmatrix(new->tensor_matrix);
//     new->tensor_matrix = c_log(z1->tensor_matrix, FREE_0);

//     if (z1->requires_grad == TRUE) {
//         new->requires_grad = TRUE;
//     } else {
//         new->requires_grad = FALSE;
//     }

//     return new;
// };

// ComplexTensor* Transpose (ComplexTensor* z1) {

//     ComplexTensor* new = createTensor(z1->tensor_matrix->shape);

//     new->creation_operation = OP_TRANSPOSE;
//     new->parents[0] = z1;
    
//     free_complexmatrix(new->tensor_matrix);
//     new->tensor_matrix = c_transpose_matrix(z1->tensor_matrix, FREE_0);

//     if (z1->requires_grad == TRUE) {
//         new->requires_grad = TRUE;
//     } else {
//         new->requires_grad = FALSE;
//     }

//     return new;
// }

// ComplexTensor* Scalar (ComplexTensor* z1, Complex* z) {

//     ComplexTensor* new = createTensor(z1->tensor_matrix->shape);

//     new->creation_operation = OP_SCALAR;
//     new->parents[0] = z1;
    
//     free_complexmatrix(new->tensor_matrix);
//     new->tensor_matrix = c_scalar(z1->tensor_matrix, z, FREE_0);

//     if (z1->requires_grad == TRUE) {
//         new->requires_grad = TRUE;
//     } else {
//         new->requires_grad = FALSE;
//     }

//     return new;
// };

// ComplexTensor* Add (ComplexTensor* z1, ComplexTensor* z2) {
//     assert (z1->tensor_matrix->shape[0] == z2->tensor_matrix->shape[0] && 
//             z1->tensor_matrix->shape[1] == z2->tensor_matrix->shape[1]);

//     ComplexTensor* z3 = createTensor(z1->tensor_matrix->shape);
//     z3->parents[0] = z1;
//     z3->parents[1] = z2;
//     z3->creation_operation = OP_ADD;
//     free_complexmatrix(z3->tensor_matrix);

//     z3->tensor_matrix = c_add_matrix(z1->tensor_matrix, z2->tensor_matrix, FREE_0);

//     if (z1->requires_grad == TRUE || z2->requires_grad == TRUE) {
//         z3->requires_grad = TRUE;
//     } else {
//         z3->requires_grad = FALSE;
//     }

//     return z3;
// };

// ComplexTensor* Sub (ComplexTensor* z1, ComplexTensor* z2) {
//     assert (z1->tensor_matrix->shape[0] == z2->tensor_matrix->shape[0] && 
//             z1->tensor_matrix->shape[1] == z2->tensor_matrix->shape[1]);

//     ComplexTensor* z3 = createTensor(z1->tensor_matrix->shape);
//     z3->parents[0] = z1;
//     z3->parents[1] = z2;
//     z3->creation_operation = OP_SUBTRACT;
//     free_complexmatrix(z3->tensor_matrix);

//     z3->tensor_matrix = c_sub_matrix(z1->tensor_matrix, z2->tensor_matrix, FREE_0);

//     if (z1->requires_grad == TRUE || z2->requires_grad == TRUE) {
//         z3->requires_grad = TRUE;
//     } else {
//         z3->requires_grad = FALSE;
//     }

//     return z3;
// };

// ComplexTensor* Mult (ComplexTensor* z1, ComplexTensor* z2) {
//     assert (z1->tensor_matrix->shape[0] == z2->tensor_matrix->shape[0] && 
//             z1->tensor_matrix->shape[1] == z2->tensor_matrix->shape[1]);

//     ComplexTensor* z3 = createTensor(z1->tensor_matrix->shape);
//     z3->parents[0] = z1;
//     z3->parents[1] = z2;
//     z3->creation_operation = OP_MULTIPLY;
//     free_complexmatrix(z3->tensor_matrix);

//     z3->tensor_matrix = c_mult_matrix(z1->tensor_matrix, z2->tensor_matrix, FREE_0);

//     if (z1->requires_grad == TRUE || z2->requires_grad == TRUE) {
//         z3->requires_grad = TRUE;
//     } else {
//         z3->requires_grad = FALSE;
//     }

//     return z3;
// };

// ComplexTensor* MatMul (ComplexTensor* z1, ComplexTensor* z2) {
//     assert (z1->tensor_matrix->shape[1] == z2->tensor_matrix->shape[0]);

//     ComplexTensor* z3 = createTensor(z1->tensor_matrix->shape);
//     z3->parents[0] = z1;
//     z3->parents[1] = z2;
//     z3->creation_operation = OP_MATMUL;
//     free_complexmatrix(z3->tensor_matrix);

//     z3->tensor_matrix = c_matmul_matrix(z1->tensor_matrix, z2->tensor_matrix, FREE_0);

//     if (z1->requires_grad == TRUE || z2->requires_grad == TRUE) {
//         z3->requires_grad = TRUE;
//     } else {
//         z3->requires_grad = FALSE;
//     }

//     return z3;
// }

// ComplexTensor* Hermitian (ComplexTensor* z) {
//     ComplexTensor* new = createTensor(z->tensor_matrix->shape);

//     new->creation_operation = OP_HERMITIAN;
//     new->parents[0] = z;
    
//     free_complexmatrix(new->tensor_matrix);
//     new->tensor_matrix = c_transpose_matrix(z->tensor_matrix, FREE_0);

//     if (z->requires_grad == TRUE) {
//         new->requires_grad = TRUE;
//     } else {
//         new->requires_grad = FALSE;
//     }

//     return new;
// }

// ComplexTensor* Ele_Pow (ComplexTensor* z1, double n) {

//     ComplexTensor* z3 = createTensor(z1->tensor_matrix->shape);
//     z3->parents[0] = z1;
//     z3->creation_operation = OP_ELE_POW;

//     free_complexmatrix(z3->tensor_matrix);
//     z3->tensor_matrix = c_ele_power(z1->tensor_matrix, n, FREE_0);

//     if (z1->requires_grad == TRUE) {
//         z3->requires_grad = TRUE;
//     } else {
//         z3->requires_grad = FALSE;
//     }

//     return z3;
// }

// void FFT (ComplexTensor* z) {
//     fft_2d(z->tensor_matrix, 0);
// }

// void iFFT (ComplexTensor* z) {
//     fft_2d(z->tensor_matrix, 1);
// }

// void backward (ComplexTensor* Z, ComplexMatrix* backward_gradient) {

//     if (Z->requires_grad == FALSE) {
//         free_complexmatrix(backward_gradient);
//         return;
//     }

//     Z->gradient = c_add_matrix(Z->gradient, backward_gradient, FREE_BOTH);

//     if (Z->parents[0] == NULL && Z->parents[1] == NULL) {
//         return;
//     }

//     switch (Z->creation_operation) {

//         case OP_ADD:
//             backward(Z->parents[0], c_mult_matrix(Z->gradient, c_ones_matrix(Z->gradient->shape), FREE_2));
//             backward(Z->parents[1], c_mult_matrix(Z->gradient, c_ones_matrix(Z->gradient->shape), FREE_2));
//             break;

//         case OP_SUBTRACT:
//             backward(Z->parents[0], c_mult_matrix(Z->gradient, c_ones_matrix(Z->gradient->shape), FREE_2));
//             backward(Z->parents[1], c_mult_matrix(Z->gradient, c_negate(c_ones_matrix(Z->gradient->shape), FREE_1), FREE_2));
//             break;
    
//         case OP_MULTIPLY:
//             backward(Z->parents[0], c_mult_matrix(Z->gradient, Z->parents[1]->tensor_matrix, FREE_0));
//             backward(Z->parents[1], c_mult_matrix(Z->gradient, Z->parents[0]->tensor_matrix, FREE_0));
//             break;
        
//         case OP_SIN:
//             backward(Z->parents[0], c_mult_matrix(Z->gradient, c_cos(Z->parents[0]->tensor_matrix, 0), FREE_2));
//             break;
        
//         case OP_COS:
//             backward(Z->parents[0], c_mult_matrix(Z->gradient, c_negate(c_sin(Z->parents[0]->tensor_matrix, FREE_0), FREE_1), FREE_2));
//             break;
        
//         case OP_MATMUL:
//             backward(Z->parents[0], c_matmul_matrix(Z->gradient, c_hermitian(Z->parents[1]->tensor_matrix, FREE_0), FREE_2));
//             backward(Z->parents[1], c_hermitian(c_matmul_matrix(c_hermitian(Z->gradient, FREE_0), Z->parents[0]->tensor_matrix, FREE_1), FREE_1));
//             break;
        
//         case OP_LOG:
//             backward( 
//                 Z->parents[0], 
//                 c_mult_matrix(
//                     Z->gradient, 
//                     c_ele_power(Z->parents[0]->tensor_matrix, -1.0f, FREE_0),
//                     FREE_2
//                 )
//             );
//             break;


//         default:
//             break;
//     }
// }

// void Backward (ComplexTensor* z) {
//     ComplexMatrix* back_grad = c_ones_matrix(z->tensor_matrix->shape);
//     backward(z, back_grad);
// }


// gcc ctensor.c complex.c -o exec