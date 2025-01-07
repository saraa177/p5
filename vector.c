#include <stdio.h>
#include <stdlib.h>
#include "vector.h"
#include "matrix.h"

struct vector_struct {
    unsigned short length;  
    double* data;          
};

/** create a (zero filled) N-dimensional double vector.
 *  @param n the length of the vector 0 < n
 *  @return a N-dimensional double vector reference
 */
vector createVector(unsigned short n){
    vector v;
    if (n == 0){
        exit(-1);
    }

    v = (vector)calloc(1, sizeof(struct vector_struct) + n * sizeof(double));
    if (v == NULL) {
        exit(-1);
    }

    v->length = n;
    v->data = (double *)v + 1;
    if (v->data == NULL) {
        printf("Error: Memory allocation failed for vector data.\n");
        free(v);
        exit(-1);
    }

    return v;
}

/**
 * remove a vector freeing allocated memory.
 *  @param v valid vector reference
 */
void rmVector(vector v){
    free(v);
}

/** @return length N of the vector */
unsigned short size(vector v){
    return v->length;
}

/** access to the internal data array. */
double* values(vector v){
    if (v != NULL && v->data != NULL) {
        return v->data;
    }
    return NULL;
}

/** get the vector entry with index j
 * @param v vector reference
 * @param j the index number  0 <= j < N
 * @return the entry
 * @print error and exits if dimensions wrong
 */
double getValue(vector v, unsigned short j){
    if (v == NULL || v->data == NULL) {
        printf("Error: Vector data is NULL.\n");
        exit(1);
    }

    if (j >= v->length) {
        printf("Error: Index %u is out of bounds. Valid range: 0 to %u.\n", j, v->length - 1);
        exit(1);
    }

    return v->data[j];
}

/** set the vector entry with index j
 * @param v vector reference
 * @param j the index number  0 <= j < N
 * @param value the entry to set
 * @print error and exits if dimensions wrong
 */
void setValue(vector v, unsigned short j, double value){
    if (v == NULL || v->data == NULL) {
        printf("Error: Vector data is NULL. Cannot set value.\n");
        exit(-1);
    }

    if (j >= v->length) {
        printf("Error: Index %u is out of bounds. Valid range: 0 to %u.\n", j, v->length - 1);
        exit(-1);
    }

    v->data[j] = value;
}

/** scalar product between two vectors.
 * @param a vector reference
 * @param b vector reference
 * @return scalar product a * b
 * @print error and exits if dimensions wrong
 */
double vectorDotVector(vector a, vector b){
    unsigned int i;
    double dot_product = 0.0;
    
    if (a == NULL || a->data == NULL || b == NULL || b->data == NULL) {
        printf("Error: One or both vectors have NULL data.\n");
        exit(1);  
    }
    
    if (a->length != b->length) {
        printf("Error: Vectors have different dimensions. a.length = %u, b.length = %u.\n", a->length, b->length);
        exit(1);  
    }

    for (i = 0; i < a->length; i++) {
        dot_product += a->data[i] * b->data[i];
    }
    
    return dot_product;
}

/** addition of two vectors.
 * @param a vector reference
 * @param b vector reference
 * @return c = a + b vector
 * @print error and exits if dimensions wrong
 */
vector vectorPlusVector(vector a, vector b){
    unsigned int i;
    vector c;

    if (a == NULL || a->data == NULL || b == NULL || b->data == NULL) {
        printf("Error: One or both vectors have NULL data. Cannot add vectors.\n");
        exit(1);  
    }

    if (a->length != b->length) {
        printf("Error: Vectors have different dimensions. a.length = %u, b.length = %u.\n", a->length, b->length);
        exit(1);  
    }

    c = createVector(a->length);
    for (i = 0; i < a->length; i++) {
        c->data[i] = a->data[i] + b->data[i];
    }

    return c;
}
