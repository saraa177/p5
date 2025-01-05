#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "matrix.h"

struct matrix_struct {
    unsigned int rows;
    unsigned int cols;
    double** data;
};

struct vector_struct {
    unsigned int length;
    double* data;
};


typedef struct matrix_struct* matrix;
typedef struct vector_struct* vector;

/** create a (zero filled) NxM double matrix. */
matrix createMatrix(unsigned short n, unsigned short m) {
    unsigned int i;
    unsigned int j;
    matrix a;
    if (n == 0 || m == 0) {
        return NULL;
    }

    a = (matrix)malloc(sizeof(struct matrix_struct));
    if (a == NULL) return NULL;

    a->rows = n;
    a->cols = m;

    
    a->data = (double**)malloc(n * sizeof(double*));
    if (a->data == NULL) {
        free(a); 
        return NULL;
    }

    for (i = 0; i < n; i++) {
        a->data[i] = (double*)malloc(m * sizeof(double));
        if (a->data[i] == NULL) {
            
            for (j = 0; j < i; j++) {
                free(a->data[j]);
            }
            free(a->data);
            free(a);
            return NULL;
        }
    }

    // Fill the matrix with zeros
    for (i = 0; i < n; i++) {
        for (unsigned int j = 0; j < m; j++) {
            a->data[i][j] = 0.0;
        }
    }

    return a;
}

/** remove a matrix freeing allocated memory. */
void rmMatrix(matrix a) {
    if (a != NULL) {
        for (unsigned int i = 0; i < a->rows; i++) {
            free(a->data[i]);
        }
        free(a->data);
        free(a);
    }
}

/** return the number of matrix rows. */
unsigned short rows(matrix a) {
    if (a != NULL) {
        return (a->rows);
    }
    return 0;
}

/** return the number of matrix columns. */
unsigned short cols(matrix a) {
    if (a != NULL) {
        return (a->cols);
    }
    return 0;
}

/** access to the internal NxM data. */
double** data(matrix a) {
    if (a == NULL)
        return NULL;

    
    double** result = (double**)malloc(a->rows * sizeof(double*));
    if (result == NULL) {
        return NULL;
    }

    for (unsigned int i = 0; i < a->rows; i++) {
        result[i] = a->data[i];
    }

    return result;
}

/** get the matrix entry with index r,c */
double getEntry(matrix a, unsigned short r, unsigned short c) {
    if (a == NULL) {
        printf("Error: Matrix is NULL.\n");
        exit(-1);
    }

    if (r >= a->rows || c >= a->cols) {
        printf("Error: Invalid row (%hu) or column (%hu) index.\n", r, c);
        exit(-1);
    }

    return a->data[r][c];
}

/** set the matrix entry with index r,c */
void setEntry(matrix a, unsigned short r, unsigned short c, double v) {
    if (a == NULL) {
        printf("Error: Matrix is NULL.\n");
        exit(-1);
    }

    if (r >= a->rows || c >= a->cols) {
        printf("Error: Invalid row (%hu) or column (%hu) index.\n", r, c);
        exit(-1);
    }

    a->data[r][c] = v;
}

/** Calculate the matrix product C = A * B. */
matrix matrixDotMatrix(matrix a, matrix b) {
    matrix c;
    unsigned int i, j, k;
    if (a == NULL || b == NULL) {
        printf("Error: One or both matrices are NULL.\n");
        exit(-1);
    }

    if (a->cols != b->rows) {
        printf("Error: Matrix dimensions do not match for multiplication.\n");
        exit(-1);
    }

    c = (matrix)malloc(sizeof(struct matrix_struct));
    if (c == NULL) {
        printf("Error: Memory allocation failed for result matrix.\n");
        exit(-1);
    }

    c->rows = a->rows;
    c->cols = b->cols;

    
    c->data = (double**)malloc(c->rows * sizeof(double*));
    for (i = 0; i < c->rows; i++) {
        c->data[i] = (double*)malloc(c->cols * sizeof(double));
    }

    
    for (i = 0; i < c->rows; i++) {
        for (j = 0; j < c->cols; j++) {
            c->data[i][j] = 0;
            for (k = 0; k < a->cols; k++) {
                c->data[i][j] += a->data[i][k] * b->data[k][j];
            }
        }
    }

    return c;
}

/** Calculate the matrix addition C = A + B. */
matrix matrixPlusMatrix(matrix a, matrix b) {
    matrix c;
    unsigned int i, j;
    if (a == NULL || b == NULL) {
        printf("Error: One or both matrices are NULL.\n");
        exit(-1);
    }

    if (a->rows != b->rows || a->cols != b->cols) {
        printf("Error: Matrix dimensions do not match for addition.\n");
        exit(-1);
    }

    c = (matrix)malloc(sizeof(struct matrix_struct));
    if (c == NULL) {
        printf("Error: Memory allocation failed for result matrix.\n");
        exit(-1);
    }

    c->rows = a->rows;
    c->cols = a->cols;

    
    c->data = (double**)malloc(c->rows * sizeof(double*));
    for (i = 0; i < c->rows; i++) {
        c->data[i] = (double*)malloc(c->cols * sizeof(double));
    }

    for (i = 0; i < c->rows; i++) {
        for (j = 0; j < c->cols; j++) {
            c->data[i][j] = a->data[i][j] + b->data[i][j];
        }
    }

    return c;
}

/** Calculate the matrix vector multiplication y = A * x. */
vector matrixDotVector(matrix a, vector x) {
    unsigned int i, j;
    vector y;
    if (a == NULL || x == NULL) {
        printf("Error: One or both matrices are NULL.\n");
        exit(-1);
    }

    if (a->cols != x->length) {
        printf("Error: Matrix dimensions do not match for multiplication.\n");
        exit(-1);
    }

    y = (vector)malloc(sizeof(struct vector_struct));
    if (y == NULL) {
        printf("Error: Memory allocation failed for result vector.\n");
        exit(-1);
    }

    y->length = a->rows;
    y->data = (double*)malloc(y->length * sizeof(double));
    if (y->data == NULL) {
        printf("Error: Memory allocation failed for result vector data.\n");
        exit(-1);
    }

    for (i = 0; i < a->rows; i++) {
        y->data[i] = 0;
        for (j = 0; j < a->cols; j++) {
            y->data[i] += a->data[i][j] * x->data[j];
        }
    }

    return y;
}
