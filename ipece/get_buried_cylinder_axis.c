#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ipece.h"

extern IPECE_PRM ipece_prm;
extern BOX box;

void get_buried_cylinder_axis() {
    int lattice_point[3];
    int mem_latt_z1;
    int mem_latt_z2;
    float coordi[3] = {0,0,0};
    float sum_cylinder_surface[2];
    float sum_cylinder_surface_sq[2];
    float sum_cylinder_radii;
    float sum_cylinder_radii_sq;
    float n_surface_points;
    float cylinder_axis[3] = {0,0,0};
    float sum_cylinder_axis[3] = {0,0,0};
    float sum_axis_multi_z[3] = {0,0,0};
    float avg_cylinder_axis[3];
    float dev_axis_multi_z[3];
    float fitted_axis[3];
    int i;
    
    coordi[2] = -ipece_prm.half_mem_thickness;
    coor2latt(coordi, lattice_point);  /* Get lattice points of membrane boundary */
    mem_latt_z1 = lattice_point[2]-ipece_prm.axis_extention;
    
    coordi[2] = ipece_prm.half_mem_thickness;
    coor2latt(coordi, lattice_point);
    mem_latt_z2 = lattice_point[2]+ipece_prm.axis_extention;
    
    for (lattice_point[2] = mem_latt_z1; lattice_point[2] <= mem_latt_z2; lattice_point[2]++ ) {
        /* Looping on Z direction, along membrane normal */
        for (i=0;i<2;i++)
            sum_cylinder_surface[i] = 0;
        for (i=0;i<2;i++)
            sum_cylinder_surface_sq[i] = 0;
        sum_cylinder_radii = 0;
        sum_cylinder_radii_sq = 0;
        n_surface_points = 0;
        for (lattice_point[0] = 0; lattice_point[0] < box.n_lattice[0]; lattice_point[0]++ ) {
            for (lattice_point[1] = 0; lattice_point[1] < box.n_lattice[1]; lattice_point[1]++ ) {
                if (box.lattice_prop[index_lattice(lattice_point)] == 's') {
                    n_surface_points++;
                    for (i=0;i<2;i++)
                        sum_cylinder_surface[i] += lattice_point[i];
                    /*
                    for (i=0;i<2;i++)
                        sum_cylinder_surface_sq[i] += lattice_point[i]*lattice_point[i];
                    radii_sq = lattice_point[0]*lattice_point[0] + lattice_point[1]*lattice_point[1]
                    sum_cylinder_radii_sq += radii_sq;
                    sum_cylinder_radii += sqrt(radii_sq);
                    */
                }
            }
        }
        
        for (i=0;i<2;i++)
            cylinder_axis[i] = sum_cylinder_surface[i] / n_surface_points;
        cylinder_axis[2] = lattice_point[2];
        /*
        avg_cylinder_radii = sum_cylinder_radii / n_surface_points;
        avg_cylinder_radii_sq = sum_cylinder_radii_sq / n_surface_points;
        */
        
        for (i=0;i<3;i++)
            sum_cylinder_axis[i] += cylinder_axis[i];
        for (i=0;i<3;i++)
            sum_axis_multi_z[i] += cylinder_axis[i] * cylinder_axis[2];
    }
    
    for (i=0;i<3;i++)
        avg_cylinder_axis[i] = sum_cylinder_axis[i]/(mem_latt_z2-mem_latt_z1+1);
    for (i=0;i<3;i++)
        dev_axis_multi_z[i] = sum_axis_multi_z[i]/(mem_latt_z2-mem_latt_z1+1) - avg_cylinder_axis[i]*avg_cylinder_axis[2];
    for (i=0;i<3;i++)
        fitted_axis[i] = dev_axis_multi_z[i] / dev_axis_multi_z[2];
    normalize_vec(fitted_axis,box.buried_cylinder_axis);
}
