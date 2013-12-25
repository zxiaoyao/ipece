#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ipece.h"

extern IPECE_PRM ipece_prm;
extern BOX box;

void coor2latt(float *coordi, int *lattice_point) {
    int i;
    
    for (i=0;i<3;i++)
        lattice_point[i] = (int)((coordi[i] - ((float)box.lower_lattice[i] * ipece_prm.lattice_scale)) / ipece_prm.lattice_scale + 0.5);
}

void latt2coor(int *lattice_point, float *coordi) {
    int i;
    
    for (i=0;i<3;i++)
        coordi[i]  = ((float)(lattice_point[i] + box.lower_lattice[i])) * ipece_prm.lattice_scale;
}

long index_lattice(int *lattice_point) {
    return 
    lattice_point[0] * box.n_lattice[1] * box.n_lattice[2]
    + lattice_point[1] * box.n_lattice[2]
    + lattice_point[2] ;
}

void idx2latt(long index, int *lattice_point) {
    long remainder;
    lattice_point[0] = index / (box.n_lattice[1] * box.n_lattice[2]);
    remainder = fmod(index, (box.n_lattice[1] * box.n_lattice[2]));
    lattice_point[1] = remainder / box.n_lattice[2];
    lattice_point[2] = fmod(remainder, box.n_lattice[2]);
}

float distsq(float *r1, float *r2) {
    float dsq = 0;
    int i;
    
    for(i=0;i<3;i++) dsq += (r1[i]-r2[i])*(r1[i]-r2[i]);
    return dsq;
}
