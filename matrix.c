#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "matrix.h"

/** create a (zero filled) NxM double matrix.
 *  @param n the number of rows, with 0<n
 *  @param m the number of columns with 0<m
 *  @return a NxM double matrix reference
 */
matrix createMatrix(unsigned short n, unsigned short m){
    unsigned int i;
    if (n == 0 || m == 0){
        return NULL;
    }

    matrix mx = (matrix)malloc(sizeof(struct matrix_struct));
    if (mx == NULL)
    return NULL;

    mx->rows = n;
    mx->cols = m;

    mx->data = (double *)malloc(n * m * sizeof(double));

    if (mx->data == NULL) {
        free(mx); 
        return NULL;
    }

    for (i = 0; i < n * m; i++) {
    mx->data[i] = 0.0;
}

}
