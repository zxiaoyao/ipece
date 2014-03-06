#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

extern IPECE_PRM ipece_prm;
extern BOX box;
#define RADIUS_MAX 2.5

void probe(ATOMS atoms, float probe_radius) 
{    
    int ia;
    BOUNDARY boundary;
    int i;
    int lattice_point[3];
    int ngh_lattice_point[3];
    int latt_z1;
    int latt_z2;
    float coordi[3];
    float ngh_coordi[3];
    float radius;
    float radsq;
    float sum_rad, sum_rad2;
    float sum_rad_sq, sum_rad2_sq;
    LONGS probes;
    int i_index;
    int iprobe;
    int new_delete = 1;
    int exposed_surface;
    
    /* find box boundary. */
    for (i=0;i<3;i++) boundary.min[i] = atoms.array[0].r[i];
    for (i=0;i<3;i++) boundary.max[i] = atoms.array[0].r[i];
    
    for (ia = 1; ia < atoms.n; ia++) {
        if (!strncmp(atoms.array[ia].pdb_line, "ATOM  ", 6)) {
            for (i=0;i<3;i++) {
                if (atoms.array[ia].r[i] < boundary.min[i] ) {
                    boundary.min[i] = atoms.array[ia].r[i];
                }
                if (atoms.array[ia].r[i] > boundary.max[i] ) {
                    boundary.max[i] = atoms.array[ia].r[i];
                }
            }
        }
    }
    
    for (i=0;i<3;i++) boundary.min[i] -= ipece_prm.boundary_extention[i]+RADIUS_MAX;
    for (i=0;i<3;i++) boundary.max[i] += ipece_prm.boundary_extention[i]+RADIUS_MAX;
    
    for (i=0;i<3;i++) {
        box.lower_lattice[i] = boundary.min[i] / ipece_prm.lattice_scale - 1;
        box.upper_lattice[i] = boundary.max[i] / ipece_prm.lattice_scale + 1;
        box.n_lattice[i] = box.upper_lattice[i] - box.lower_lattice[i] + 1;
    }
    
    // labels of all lattice points
    box.lattice_prop = realloc(box.lattice_prop, box.n_lattice[0]*box.n_lattice[1]*box.n_lattice[2]*sizeof(char));
    memset(box.lattice_prop,0,box.n_lattice[0]*box.n_lattice[1]*box.n_lattice[2]*sizeof(char));
    
    for (ia=0; ia<atoms.n; ia++) {
    	// only "ATOM" is taken into account, not "HEATOM"
        if (strncmp(atoms.array[ia].pdb_line, "ATOM  ", 6)) continue;
        coor2latt(atoms.array[ia].r,lattice_point);
        
        radius = get_rad(atoms.array[ia]);
        sum_rad = radius+probe_radius;
        sum_rad2 = radius+2.*probe_radius;
        radsq  = radius*radius;
        sum_rad_sq = sum_rad*sum_rad;
        sum_rad2_sq = sum_rad2*sum_rad2;
        for (ngh_lattice_point[0] = lattice_point[0]-sum_rad2/ipece_prm.lattice_scale; ngh_lattice_point[0]<=lattice_point[0]+sum_rad2/ipece_prm.lattice_scale; ngh_lattice_point[0]++) {
            for (ngh_lattice_point[1] = lattice_point[1]-sum_rad2/ipece_prm.lattice_scale; ngh_lattice_point[1]<=lattice_point[1]+sum_rad2/ipece_prm.lattice_scale; ngh_lattice_point[1]++) {
                for (ngh_lattice_point[2] = lattice_point[2]-sum_rad2/ipece_prm.lattice_scale; ngh_lattice_point[2]<=lattice_point[2]+sum_rad2/ipece_prm.lattice_scale; ngh_lattice_point[2]++) {
                    float distsq_2ngh;
                    latt2coor(ngh_lattice_point,ngh_coordi);
                    distsq_2ngh = distsq(atoms.array[ia].r, ngh_coordi);
                    
                    //printf("%d, %d\n",index_lattice(ngh_lattice_point), box.n_lattice[0]*box.n_lattice[1]*box.n_lattice[2]);
                    //if (distsq(atoms.array[ia].r, ngh_coordi) < radsq) {
                    // 'p' refers to being inside the atom.
                    if (distsq_2ngh <= radsq) {
                        box.lattice_prop[index_lattice(ngh_lattice_point)] = 'p';
                    }
                    else if (distsq_2ngh <= sum_rad_sq) {
                    	// 's' means "surface"?
                        if ( box.lattice_prop[index_lattice(ngh_lattice_point)] != 'p' ) {
                            box.lattice_prop[index_lattice(ngh_lattice_point)] = 's';
                        }
                    }
                    else if (distsq_2ngh <= sum_rad2_sq) {
                        if ( box.lattice_prop[index_lattice(ngh_lattice_point)] == 0 ) {
                            box.lattice_prop[index_lattice(ngh_lattice_point)] = 'h';
                        }
                    }
                }
            }
        }
    }
    //write_probe("testprot.pdb");
    // documentation is really important, now it's really hard for me to understand this piece of code without it.
    if (!ipece_prm.inner_mem) {
        int latt_mem_z1, latt_mem_z2;
        coordi[2] = -6.;
        coor2latt(coordi, lattice_point);  /* Get lattice points of mid part of membrane */
        latt_z1 = lattice_point[2];
        
        coordi[2] = -coordi[2];
        coor2latt(coordi, lattice_point);
        latt_z2 = lattice_point[2];
        
        coordi[2] = -ipece_prm.half_mem_thickness-ipece_prm.axis_extention;  // make sure it's still in the box
        coor2latt(coordi, lattice_point);  /* Get lattice points of membrane boundary */
        latt_mem_z1 = lattice_point[2];
        
        coordi[2] = -coordi[2];
        coor2latt(coordi, lattice_point);
        latt_mem_z2 = lattice_point[2];
        
        for (lattice_point[0] = box.n_lattice[0]-1; lattice_point[0]>=0; lattice_point[0]--) {
            for (lattice_point[1] = box.n_lattice[1]-1; lattice_point[1]>=0; lattice_point[1]--) {
                float temp_counter = 0.;
                for (lattice_point[2] = latt_z2; lattice_point[2]>=latt_z1; lattice_point[2]--) {
                    if ( box.lattice_prop[index_lattice(lattice_point)] == 'p' ) temp_counter += 1.0;
                    else if ( box.lattice_prop[index_lattice(lattice_point)] == 's' ) temp_counter += 1.0;
                }
                
                /* if mostly protein in the middle [-6, 6], fill the rest of grid point on the same z coordinates as 'p' */
				// because anyway there is no inner membrane placed inside
                if (temp_counter / (float)(latt_z2-latt_z1+1) >= 0.9) {
                    for (lattice_point[2] = latt_mem_z2; lattice_point[2]>=latt_mem_z1; lattice_point[2]--) {
                        if (box.lattice_prop[index_lattice(lattice_point)] != 'p') box.lattice_prop[index_lattice(lattice_point)] = 's';
                    }
                }
            }
        }
    }
    
    /* start from each edge to find solvent exposed grids, without probing */
    if (ipece_prm.inner_mem) {
		// only test the upper part of the protein, 'o' means outside
        for (lattice_point[0] = box.n_lattice[0]-1; lattice_point[0]>=0; lattice_point[0]--) {
            for (lattice_point[1] = box.n_lattice[1]-1; lattice_point[1]>=0; lattice_point[1]--) {
                for (lattice_point[2] = box.n_lattice[2]-1; lattice_point[2]>=0; lattice_point[2]--) {
                    if ( box.lattice_prop[index_lattice(lattice_point)] == 0 ) 
                        box.lattice_prop[index_lattice(lattice_point)] = 'o';
                    else
                        break;
                }
            }
        }
        for (lattice_point[1] = box.n_lattice[1]-1; lattice_point[1]>=0; lattice_point[1]--) {
            for (lattice_point[2] = box.n_lattice[2]-1; lattice_point[2]>=0; lattice_point[2]--) {
                for (lattice_point[0] = box.n_lattice[0]-1; lattice_point[0]>=0; lattice_point[0]--) {
                    if ( box.lattice_prop[index_lattice(lattice_point)] == 0 ) 
                        box.lattice_prop[index_lattice(lattice_point)] = 'o';
                    else
                        break;
                }
            }
        }
        for (lattice_point[2] = box.n_lattice[2]-1; lattice_point[2]>=0; lattice_point[2]--) {
            for (lattice_point[0] = box.n_lattice[0]-1; lattice_point[0]>=0; lattice_point[0]--) {
                for (lattice_point[1] = box.n_lattice[1]-1; lattice_point[1]>=0; lattice_point[1]--) {
                    if ( box.lattice_prop[index_lattice(lattice_point)] == 0 ) 
                        box.lattice_prop[index_lattice(lattice_point)] = 'o';
                    else
                        break;
                }
            }
        }
    }
    else {
        coordi[2] = -ipece_prm.half_mem_thickness-ipece_prm.axis_extention;
        coor2latt(coordi, lattice_point);  /* Get lattice points of membrane boundary */
        latt_z1 = lattice_point[2];
        
        coordi[2] = ipece_prm.half_mem_thickness+ipece_prm.axis_extention;
        coor2latt(coordi, lattice_point);
        latt_z2 = lattice_point[2];
        
        for (lattice_point[1] = box.n_lattice[1]-1; lattice_point[1]>=0; lattice_point[1]--) {
            for (lattice_point[2] = latt_z2; lattice_point[2]>=latt_z1; lattice_point[2]--) {
                for (lattice_point[0] = box.n_lattice[0]-1; lattice_point[0]>=0; lattice_point[0]--) {
                    if ( box.lattice_prop[index_lattice(lattice_point)] == 0 ) 
                        box.lattice_prop[index_lattice(lattice_point)] = 'o';
                    else
                        break;
                }
            }
        }
        for (lattice_point[2] = latt_z2; lattice_point[2]>=latt_z1; lattice_point[2]--) {
            for (lattice_point[0] = box.n_lattice[0]-1; lattice_point[0]>=0; lattice_point[0]--) {
                for (lattice_point[1] = box.n_lattice[1]-1; lattice_point[1]>=0; lattice_point[1]--) {
                    if ( box.lattice_prop[index_lattice(lattice_point)] == 0 ) 
                        box.lattice_prop[index_lattice(lattice_point)] = 'o';
                    else
                        break;
                }
            }
        }
    }
	
	
    memset(&probes,0,sizeof(LONGS));
    for (i_index = 0; i_index < box.n_lattice[0]*box.n_lattice[1]*box.n_lattice[2]; i_index++) {
        if ( box.lattice_prop[i_index] == 0 || (box.lattice_prop[i_index] == 'h' && ipece_prm.add_mem)) {
            if (!ipece_prm.inner_mem) {
                idx2latt(i_index,lattice_point);
                if (lattice_point[2]<latt_z1 || lattice_point[2]>latt_z2) {
                    if (box.lattice_prop[i_index] == 0) box.lattice_prop[i_index] = 'u'; /* solution */
                    continue;
                }
            }
            probes.n++;
            probes.array = realloc(probes.array, probes.n*sizeof(long));
            probes.array[probes.n-1] = i_index;
        }
    }
    
    while (new_delete) {
        /* 
        START! Debug/Checking 
        printf("Number of probes left: %i\n",probes.n);
        END!   Debug/Checking
        */
        new_delete = 0;
        for (iprobe=probes.n-1;iprobe>=0;iprobe--) {
            /* If one of the neighbor lattice points belongs to outside, the probe is outside. */
            for (i=0;i<3;i++)
                ngh_lattice_point[i] = -1;
            while (box.lattice_prop[probes.array[iprobe]+index_lattice(ngh_lattice_point)] != 'o') {
                ngh_lattice_point[2]++;
                
                if (ngh_lattice_point[2]>1) {
                    ngh_lattice_point[2] -= 3;
                    ngh_lattice_point[1]++;
                    
                    if (ngh_lattice_point[1]>1) {
                        ngh_lattice_point[1] -= 3;
                        ngh_lattice_point[0]++;
                        
                        if (ngh_lattice_point[0]>1) {
                            ngh_lattice_point[0] -= 3;
                            break;
                        }
                    }
                }
            }
            if ( box.lattice_prop[probes.array[iprobe]+index_lattice(ngh_lattice_point)] == 'o' ) {
                box.lattice_prop[probes.array[iprobe]] = 'o';
                memmove(&probes.array[iprobe], &probes.array[iprobe+1], (probes.n-iprobe-1)*sizeof(long));
                probes.n--;
                new_delete = 1;
            }
        }
    }
    
    for (iprobe = 0; iprobe < probes.n; iprobe++) {
        if (box.lattice_prop[probes.array[iprobe]] == 0) box.lattice_prop[probes.array[iprobe]] = 'c';
    }
    
    for (lattice_point[0] = 0; lattice_point[0] < box.n_lattice[0]; lattice_point[0]++) {
        for (lattice_point[1] = 0; lattice_point[1] < box.n_lattice[1]; lattice_point[1]++) {
            for (lattice_point[2] = 0; lattice_point[2] < box.n_lattice[2]; lattice_point[2]++) {
                if (!ipece_prm.inner_mem) {
                    if (lattice_point[2] < latt_z1) continue;
                    if (lattice_point[2] > latt_z2) continue;
                }
                
                if (box.lattice_prop[index_lattice(lattice_point)] == 'h') {
                    exposed_surface = 0;
                    /* 
                    START! Debug/Checking 
                    printf("debug\n");
                    END!   Debug/Checking 
                    */
                    
                    for (ngh_lattice_point[0] = lattice_point[0]-1.1*probe_radius/ipece_prm.lattice_scale; ngh_lattice_point[0]<=lattice_point[0]+1.1*probe_radius/ipece_prm.lattice_scale; ngh_lattice_point[0]++) {
                        for (ngh_lattice_point[1] = lattice_point[1]-1.1*probe_radius/ipece_prm.lattice_scale; ngh_lattice_point[1]<=lattice_point[1]+1.1*probe_radius/ipece_prm.lattice_scale; ngh_lattice_point[1]++) {
                            for (ngh_lattice_point[2] = lattice_point[2]-1.1*probe_radius/ipece_prm.lattice_scale; ngh_lattice_point[2]<=lattice_point[2]+1.1*probe_radius/ipece_prm.lattice_scale; ngh_lattice_point[2]++) {
                                if (index_lattice(ngh_lattice_point) < 0) continue;
                                if (index_lattice(ngh_lattice_point) >= box.n_lattice[0]*box.n_lattice[1]*box.n_lattice[2]) continue;
                                //printf("%d %d %d %d\n",ngh_lattice_point[0], ngh_lattice_point[1], ngh_lattice_point[2], index_lattice(ngh_lattice_point));
                                if (box.lattice_prop[index_lattice(ngh_lattice_point)] == 'o') {
                                    exposed_surface = 1;
                                    break;
                                }
                            }
                            if (exposed_surface) break;
                        }
                        if (exposed_surface) break;
                    }
                    if (!exposed_surface) box.lattice_prop[index_lattice(lattice_point)] = 'i';
                }
            }
        }
    }
}

void write_probe(char *file_name) {
    float r[3];
    int counter = 0;
    int i_index, lattice_point[3];
    FILE *pdb_out_fp;
    char pdbline[160];
    pdb_out_fp = fopen(file_name, "w");
    for (i_index = 0; i_index < box.n_lattice[0]*box.n_lattice[1]*box.n_lattice[2]; i_index++) {
        counter++;
        idx2latt(i_index,lattice_point);
        latt2coor(lattice_point, r);
        
        //if (box.lattice_prop[i_index] == 0) box.lattice_prop[i_index] = '_';
        
        int resSeq = counter/1000.;
        while (resSeq >9999) resSeq -= 10000;
        
        sprintf(pdbline, "ATOM  00000  %c     %c Y%4d    %8.3f%8.3f%8.3f",
                            box.lattice_prop[i_index],
                            box.lattice_prop[i_index],
                            resSeq,
                            r[0], r[1], r[2]);
                            
        fprintf(pdb_out_fp,"%s\n",pdbline);

    }
    fclose(pdb_out_fp);

}
