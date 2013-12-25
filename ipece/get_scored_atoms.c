
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

extern IPECE_PRM ipece_prm;
extern BOX box;

ATOMS get_scored_atoms(ATOMS atoms) 
{
    ATOMS scored_atoms;
    int   ia;
    int   i_sa_type;
    int   scored;
    int   lattice_point[3];
    int   ngh_lattice_point[3];
    
    memset(&scored_atoms,0,sizeof(ATOMS));
    
    for (ia=0; ia<atoms.n; ia++) {  /* Loop over all atoms. */
        for (i_sa_type = 0; i_sa_type < ipece_prm.n_sa_type; i_sa_type++) {
            /* Loop over the scored atom type list. */
            /* Check if "ia" belongs to the scored atom list */
            if (!strncmp(atoms.array[ia].pdb_line+17, ipece_prm.sa_type[i_sa_type].res,3)
             && !strncmp(atoms.array[ia].pdb_line+12, ipece_prm.sa_type[i_sa_type].atom,4) ) {
                coor2latt(atoms.array[ia].r, lattice_point);
                scored = 0;
                /* only put the atom into the score list if it is within
                surface_exp_rad to the solvent. this condition is added for
                the outer membrane proteins, where there are many exposed ionizable
                residues in the transmembrane region facing the pore side */
                for (ngh_lattice_point[0] = lattice_point[0]-ipece_prm.surface_exp_rad/ipece_prm.lattice_scale; ngh_lattice_point[0]<=lattice_point[0]+ipece_prm.surface_exp_rad/ipece_prm.lattice_scale; ngh_lattice_point[0]++) {
                    for (ngh_lattice_point[1] = lattice_point[1]-ipece_prm.surface_exp_rad/ipece_prm.lattice_scale; ngh_lattice_point[1]<=lattice_point[1]+ipece_prm.surface_exp_rad/ipece_prm.lattice_scale; ngh_lattice_point[1]++) {
                        for (ngh_lattice_point[2] = lattice_point[2]-ipece_prm.surface_exp_rad/ipece_prm.lattice_scale; ngh_lattice_point[2]<=lattice_point[2]+ipece_prm.surface_exp_rad/ipece_prm.lattice_scale; ngh_lattice_point[2]++) {
                            if (box.lattice_prop[index_lattice(ngh_lattice_point)] == 'o'
                            || box.lattice_prop[index_lattice(ngh_lattice_point)] == 's') {
                                scored = 1;
                                scored_atoms.n++;
                                scored_atoms.array = realloc(scored_atoms.array, scored_atoms.n*sizeof(ATOM));
                                memcpy(&scored_atoms.array[scored_atoms.n-1], &atoms.array[ia], sizeof(ATOM));
                                scored_atoms.array[scored_atoms.n-1].score = ipece_prm.sa_type[i_sa_type].buried_score * scored_atoms.array[scored_atoms.n-1].acc;
                                break;
                            }
                        }
                        if (scored) break;
                    }
                    if (scored) break;
                }
            }
        }
    }
    return scored_atoms;
}
