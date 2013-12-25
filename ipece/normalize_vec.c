#include <math.h>

void normalize_vec(float *input_vec, float *normalized_vec) {
    
    int i;
    float norm = 0;
    
    for (i=0; i<3; i++) {
        norm += input_vec[i]*input_vec[i];
    }
    norm = sqrt(norm);
    
    if (norm > 10e-8) {
        for (i=0; i<3; i++) {
            normalized_vec[i] = input_vec[i]/norm;
        }
    }
    else {
        for (i=0; i<3; i++) {
            normalized_vec[i] = 0;
        }
    }
}

