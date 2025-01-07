#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

struct matrix_struct {
    unsigned int rows;
    unsigned int cols;
    double** data;
};

/** create a (zero filled) NxM double matrix. */
matrix createMatrix(unsigned short n, unsigned short m) {
    matrix mat = calloc(1, sizeof(struct matrix_struct) + n * sizeof(double *) + m * n * sizeof(double));
    unsigned short i;
    mat->rows = n;
    mat->cols = m;
    mat->data = (double **)(mat + 1);

    for (i = 0; i < n; i++) {
        mat->data[i] = (double *)mat->data + n + i * m;
    }
    return mat;
}

/** remove a matrix freeing allocated memory. */
void rmMatrix(matrix a) {
    free(a);
}

/** return the number of matrix rows. */
unsigned short rows(matrix a) {
    return a->rows;
}

/** return the number of matrix columns. */
unsigned short cols(matrix a) {
    return a->cols;
}

/** access to the internal NxM data. */
double** data(matrix a) {
    return a->data;
}

/** get the matrix entry with index r,c */
double getEntry(matrix a, unsigned short r, unsigned short c) {

    if (r >= a->rows || c >= a->cols) {
        printf("Error: Wrong dimensions for matrix");
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

    if (a->cols != b->rows || a->rows != b->cols) {
        printf("Error: Matrix dimensions do not match for multiplication.\n");
        exit(-1);
    }

    c = createMatrix(a->rows, b->cols);
    if (c == NULL) {
        printf("Error: Memory allocation failed for result matrix.\n");
        exit(-1);
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

    c = createMatrix(a->rows, a->cols);
    if (c == NULL) {
        printf("Error: Memory allocation failed for result matrix.\n");
        exit(-1);
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

    if (a->cols != size(x)) {
        printf("Error: Matrix dimensions do not match for multiplication.\n");
        exit(-1);
    }

    y = createVector(a->rows);
    if (y == NULL) {
        printf("Error: Memory allocation failed for result vector.\n");
        exit(-1);
    }

    for (i = 0; i < a->rows; i++) {
        setValue(y, i, 0);
        for (j = 0; j < a->cols; j++) {
            setValue(y, i, getValue(y, i) + a->data[i][j] * getValue(x, j));
        }
    }

    return y;
}
