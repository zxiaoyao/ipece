#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ipece.h"

extern BOX box;

void rotate(ATOMS atoms, float *rotation) 
{
    
    float rotation_matrix[3][3];
    float r_old[3];
    float cos_alpha = cos(rotation[3]);
    float sin_alpha = sin(rotation[3]);
    int ia;
    int i;
    int j;
    //printf("rotation component, x:%.2f, y:%.2f, z:%.2f, angle:%.2f\n", rotation[0], rotation[1], rotation[2], rotation[3]);
    rotation_matrix[0][0] = rotation[0]*rotation[0]*(1-cos_alpha) + cos_alpha;
    rotation_matrix[0][1] = rotation[0]*rotation[1]*(1-cos_alpha) - rotation[2]*sin_alpha;
    rotation_matrix[0][2] = rotation[0]*rotation[2]*(1-cos_alpha) + rotation[1]*sin_alpha;
    
    rotation_matrix[1][0] = rotation[0]*rotation[1]*(1-cos_alpha) + rotation[2]*sin_alpha;
    rotation_matrix[1][1] = rotation[1]*rotation[1]*(1-cos_alpha) + cos_alpha;
    rotation_matrix[1][2] = rotation[1]*rotation[2]*(1-cos_alpha) - rotation[0]*sin_alpha;
    
    rotation_matrix[2][0] = rotation[0]*rotation[2]*(1-cos_alpha) - rotation[1]*sin_alpha;
    rotation_matrix[2][1] = rotation[1]*rotation[2]*(1-cos_alpha) + rotation[0]*sin_alpha;
    rotation_matrix[2][2] = rotation[2]*rotation[2]*(1-cos_alpha) + cos_alpha;
    
    for (ia = 0; ia<atoms.n; ia++) {
        for(i=0;i<3;i++) r_old[i] = atoms.array[ia].r[i];
        for(i=0;i<3;i++) atoms.array[ia].r[i] = 0;
        
        for(i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                atoms.array[ia].r[i] += rotation_matrix[i][j] * r_old[j];
            }
        }
    }
}

void rotate_axis(float *rotation) 
{
    float rotation_matrix[3][3];
    float old_axis[3];
    float cos_alpha = cos(rotation[3]);
    float sin_alpha = sin(rotation[3]);
    int i;
    int j;
    
    rotation_matrix[0][0] = rotation[0]*rotation[0]*(1-cos_alpha) + cos_alpha;
    rotation_matrix[0][1] = rotation[0]*rotation[1]*(1-cos_alpha) - rotation[2]*sin_alpha;
    rotation_matrix[0][2] = rotation[0]*rotation[2]*(1-cos_alpha) + rotation[1]*sin_alpha;
    
    rotation_matrix[1][0] = rotation[0]*rotation[1]*(1-cos_alpha) + rotation[2]*sin_alpha;
    rotation_matrix[1][1] = rotation[1]*rotation[1]*(1-cos_alpha) + cos_alpha;
    rotation_matrix[1][2] = rotation[1]*rotation[2]*(1-cos_alpha) - rotation[0]*sin_alpha;
    
    rotation_matrix[2][0] = rotation[0]*rotation[2]*(1-cos_alpha) - rotation[1]*sin_alpha;
    rotation_matrix[2][1] = rotation[1]*rotation[2]*(1-cos_alpha) + rotation[0]*sin_alpha;
    rotation_matrix[2][2] = rotation[2]*rotation[2]*(1-cos_alpha) + cos_alpha;
    
    for(i=0;i<3;i++) old_axis[i] = box.buried_cylinder_axis[i];
    for(i=0;i<3;i++) box.buried_cylinder_axis[i] = 0;
    
    for(i=0;i<3;i++) {
        for(j=0;j<3;j++) {
            box.buried_cylinder_axis[i] += rotation_matrix[i][j] * old_axis[j];
        }
    }
}

