#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ipece.h"

extern IPECE_PRM ipece_prm;
extern BOX box;

void add_mem_atoms(ATOMS *atoms_p) 
{
    int   lattice_point[3];
    float coordi[3] = {0,0,0};
    int   ngh_lattice_point[3];
    float ngh_coordi[3];
    int mem_latt1[3] = {-1,-1,-1};
    int mem_latt2[3] = {-1,-1,-1};
    int added;
    int n_mem_atoms = 0;
    int i;
    char sbuffer[160];

    probe(*atoms_p, ipece_prm.mem_radius);

    coordi[2] = -ipece_prm.half_mem_thickness;
    coor2latt(coordi, lattice_point);  /* Get lattice points of membrane boundary */
    mem_latt1[2] = lattice_point[2];
    
    coordi[2] = ipece_prm.half_mem_thickness;
    coor2latt(coordi, lattice_point);
    mem_latt2[2] = lattice_point[2];

    for (lattice_point[0] = 0; lattice_point[0] < box.n_lattice[0]; lattice_point[0]++) {
        for (lattice_point[1] = 0; lattice_point[1] < box.n_lattice[1]; lattice_point[1]++) {
            for (lattice_point[2] = mem_latt1[2]; lattice_point[2] <= mem_latt2[2]; lattice_point[2]++) {
                if (box.lattice_prop[index_lattice(lattice_point)] == 'p') {
                    for (i=0;i<2;i++) {
                        if (mem_latt1[i] == -1) 
                            mem_latt1[i] = lattice_point[i];
                        if (mem_latt2[i] == -1) 
                            mem_latt2[i] = lattice_point[i];
                        
                        if (lattice_point[i] < mem_latt1[i])
                            mem_latt1[i] = lattice_point[i];
                        if (lattice_point[i] > mem_latt2[i])
                            mem_latt2[i] = lattice_point[i];
                    }
                }
            }
        }
    }
    
    for (i=0;i<2;i++) {
        mem_latt1[i] -= ipece_prm.boundary_extention[i]/ipece_prm.lattice_scale;
        mem_latt2[i] += ipece_prm.boundary_extention[i]/ipece_prm.lattice_scale;
    }
    
    for (lattice_point[0] = mem_latt1[0]; lattice_point[0] < mem_latt2[0]; lattice_point[0]++) {
        for (lattice_point[1] = mem_latt1[1]; lattice_point[1] < mem_latt2[1]; lattice_point[1]++) {
            for (lattice_point[2] = mem_latt1[2]; lattice_point[2] <= mem_latt2[2]; lattice_point[2]++) {
                if (box.lattice_prop[index_lattice(lattice_point)] == 'o') {
                    latt2coor(lattice_point,coordi);
                    added = 0; /* added = 1 if another atom within certain distance has been added (to make sure added atoms are seperated by given distance. */
                    for (ngh_lattice_point[0] = lattice_point[0]-ipece_prm.mem_place_sep/ipece_prm.lattice_scale; ngh_lattice_point[0]<=lattice_point[0]+ipece_prm.mem_place_sep/ipece_prm.lattice_scale; ngh_lattice_point[0]++) {
                        for (ngh_lattice_point[1] = lattice_point[1]-ipece_prm.mem_place_sep/ipece_prm.lattice_scale; ngh_lattice_point[1]<=lattice_point[1]+ipece_prm.mem_place_sep/ipece_prm.lattice_scale; ngh_lattice_point[1]++) {
                            for (ngh_lattice_point[2] = lattice_point[2]-ipece_prm.mem_place_sep/ipece_prm.lattice_scale; ngh_lattice_point[2]<=lattice_point[2]+ipece_prm.mem_place_sep/ipece_prm.lattice_scale; ngh_lattice_point[2]++) {
                                latt2coor(ngh_lattice_point,ngh_coordi);
                                if (distsq(coordi, ngh_coordi) < ipece_prm.mem_place_sep) {
                                    
					if (ngh_lattice_point[0]>0 && ngh_lattice_point[1]>0 && ngh_lattice_point[2]>0) {
                                    if (box.lattice_prop[index_lattice(ngh_lattice_point)] == 'a'
                                    || box.lattice_prop[index_lattice(ngh_lattice_point)] == 'p') {
                                        added = 1;
                                        break;
                                    }
					}
                                }
                            }
                            if (added) break;
                        }
                        if (added) break;
                    }
                    if (!added) {
                        box.lattice_prop[index_lattice(lattice_point)] = 'a';
                        
                        atoms_p->n++;
                        n_mem_atoms++;
                        atoms_p->array = realloc(atoms_p->array, atoms_p->n * sizeof(ATOM));
                        memset(&atoms_p->array[atoms_p->n-1],0,sizeof(ATOM));
                        sprintf(atoms_p->array[atoms_p->n-1].pdb_line,"HETATM      ");
                        strcat(atoms_p->array[atoms_p->n-1].pdb_line,ipece_prm.mem_name.txt);
                        sprintf(atoms_p->array[atoms_p->n-1].pdb_line+26,"    %8.3f%8.3f%8.3f\n",coordi[0],coordi[1],coordi[2]);
                        char sbuffer[160];
                        sprintf(sbuffer,"%06i",n_mem_atoms+999);
                        strncpy(atoms_p->array[atoms_p->n-1].pdb_line+23,sbuffer,3);
                        atoms_p->array[atoms_p->n-1].pdb_line[12]=sbuffer[3];
                        strncpy(atoms_p->array[atoms_p->n-1].pdb_line+14,sbuffer+4,2);
                        for(i=0;i<3;i++)
                            atoms_p->array[atoms_p->n-1].r[i] = coordi[i];
                        
                        /*
                        atoms_p->n++;
                        n_mem_atoms++;
                        atoms_p->array = realloc(atoms_p->array, atoms_p->n * sizeof(ATOM));
                        memset(&atoms_p->array[atoms_p->n-1],0,sizeof(ATOM));
                        sprintf(atoms_p->array[atoms_p->n-1].pdb_line,"HETATM      ");
                        strcat(atoms_p->array[atoms_p->n-1].pdb_line,ipece_prm.mem_name.txt);
                        sprintf(atoms_p->array[atoms_p->n-1].pdb_line+26,"    %8.3f%8.3f%8.3f\n",coordi[0],coordi[1],coordi[2]);
                        sprintf(sbuffer,"%06i",n_mem_atoms+999);
                        strncpy(atoms_p->array[atoms_p->n-1].pdb_line+23,sbuffer,3);
                        atoms_p->array[atoms_p->n-1].pdb_line[12]=sbuffer[3];
                        strncpy(atoms_p->array[atoms_p->n-1].pdb_line+14,sbuffer+4,2);

                        for(i=0;i<3;i++)
                            atoms_p->array[atoms_p->n-1].r[i] = coordi[i];
                        */
                    }
                }
            }
        }
    }
    printf("%i membrane atoms added\n",n_mem_atoms);
}

                    

