#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ipece.h"

IPECE_PRM ipece_prm;
BOX box;

int main(int argc, char *argv[]) 
{
    
    ATOMS all_atoms;
    ATOMS scored_atoms;
    HISTORIES histories;
    float central_point[3];
    float neg_central[3];
    int i;
    float score;
    float new_score;
    float min_score;
    float deltaS;
    int intl_case;
    float translation[3];
    float rotation[4];
    float phi;
    int i_iter;
    char sbuffer[TXT_LEN_MAX];
    int ia;
    int i_hist;
    int i_ion;
    long dummy = 45031240;
    long *idum = &dummy;
    char *pdb_out_fname;
    FILE *pdb_out_fp;
    
    if (argc<4) {
        printf("Usage:\n    ipece parameter_file input_pdb output_pdb\nNotice:\n    You also need to have acc.atm file in the same directory, which is output from program SURFV\n");
        exit(0);
    }
    
    printf("\n     IPECE - Implement Protein Environment for Continumm Electrostatics\n");
    printf("                           Version 1.0.3\n\n");
    memset(&histories,0,sizeof(HISTORIES));
    memset(&box, 0, sizeof(BOX));
    ipece_prm = read_ipece_prm(argv[1]);
    all_atoms = read_atoms(argv[2]);
    pdb_out_fname = argv[3];
    if (!all_atoms.n) {
        printf("No atoms read in from %s\n", argv[2]);
        exit(-1);
    }
    printf("Number of atoms %i\n",all_atoms.n);
    
    /* 
    for (i=0;i<all_atoms.n;i++) {
        printf("%s\n",all_atoms.array[i].pdb_line);
    }
    */
    
    if (ipece_prm.add_mem && ipece_prm.spe_mem == 0) {
        read_acc_atom(all_atoms, "acc.atm");
        find_center_nonhet(all_atoms, central_point);
        for(i=0;i<3;i++) neg_central[i] = -central_point[i];
        
        translate(all_atoms, neg_central);
        save_hist('t', neg_central, ipece_prm.axis_score_weight, &histories);
        probe(all_atoms,ipece_prm.mem_radius);
        scored_atoms = get_scored_atoms(all_atoms);  /* make another list and copy of atoms. */
        
        min_score = get_score(scored_atoms);
        intl_case = 1;
        printf("Score 1 = %f\n", min_score);
        
        rotation[0] = 0;
        rotation[1] = -1;
        rotation[2] = 0;
        rotation[3] = ipece_prm.PI/2;
        
        if (ipece_prm.inner_mem) {
            backup(scored_atoms);
            rotate(scored_atoms, rotation);
            score = get_score(scored_atoms);
            recover(scored_atoms);
        }
        else {
            backup(all_atoms);
            rotate(all_atoms, rotation);
            probe(all_atoms,ipece_prm.mem_radius);
            scored_atoms = get_scored_atoms(all_atoms);
            score = get_score(scored_atoms);
            recover(all_atoms);
        }
        
        printf("Score 2 = %f\n", score);
        if (score < min_score) {
            min_score = score;
            intl_case = 2;
        }
        
        rotation[0] = 1;
        rotation[1] = 0;
        rotation[2] = 0;
        rotation[3] = ipece_prm.PI/2;
        
        if (ipece_prm.inner_mem) {
            backup(scored_atoms);
            rotate(scored_atoms, rotation);
            score = get_score(scored_atoms);
            recover(scored_atoms);
        }
        else {
            backup(all_atoms);
            rotate(all_atoms, rotation);
            probe(all_atoms,ipece_prm.mem_radius);
            scored_atoms = get_scored_atoms(all_atoms);
            score = get_score(scored_atoms);
            recover(all_atoms);
        }
        
        printf("Score 3 = %f\n", score);
        if (score < min_score) {
            min_score = score;
            intl_case = 3;
        }
        
        printf("Initial orientation case = %i\n", intl_case);
        if (intl_case == 2) {
            rotation[0] = 0;
            rotation[1] = -1;
            rotation[2] = 0;
            rotation[3] = ipece_prm.PI/2;
            rotate(all_atoms, rotation);
            rotate(scored_atoms, rotation);
            save_hist('r', rotation, score, &histories);
        }
        else if (intl_case == 3) {
            rotation[0] = 1;
            rotation[1] = 0;
            rotation[2] = 0;
            rotation[3] = ipece_prm.PI/2;
            rotate(all_atoms, rotation);
            rotate(scored_atoms, rotation);
            save_hist('r', rotation, score, &histories);
        }
        
        probe(all_atoms,ipece_prm.mem_radius);
        get_buried_cylinder_axis();
        if (!ipece_prm.inner_mem) {
            scored_atoms = get_scored_atoms(all_atoms);
        }
        printf("List of atoms being scored:\n");
        for(ia = 0; ia<scored_atoms.n; ia++) {
            strncpy(sbuffer,scored_atoms.array[ia].pdb_line+13,14);
            sbuffer[14] = '\0';
            printf("%s acc=%8.3f\n",sbuffer,scored_atoms.array[ia].acc);
        }
        printf("%f,%f,%f\n", box.buried_cylinder_axis[0],box.buried_cylinder_axis[1],box.buried_cylinder_axis[2]);
        
        score = get_score(scored_atoms);
        
        translation[0]=0;
        translation[1]=0;
        for (i_iter=0;i_iter<ipece_prm.n_iteration;i_iter++) {
            if (ran2(idum) < 0.5) {
                translation[2] = 2.*(ran2(idum) - 0.5)*ipece_prm.translation_max;
                backup(scored_atoms);
                translate(scored_atoms, translation);
                new_score = get_score(scored_atoms);
                deltaS = exp(-ipece_prm.beta * (new_score - score));
                if (ran2(idum) < deltaS) {
                    translate(all_atoms, translation);
                    score = new_score;
                    save_hist('t', translation, score, &histories);
                }
                else {
                    recover(scored_atoms);
                }
            }
            else {
                phi = ran2(idum)*2.*ipece_prm.PI;
                rotation[0] = cos(phi);
                rotation[1] = sin(phi);
                rotation[2] = 0.;
                rotation[3] = 2.*(ran2(idum)-0.5)*ipece_prm.rotation_max;
                backup(scored_atoms);
                backup_axis();
                
                rotate(scored_atoms,rotation);
                rotate_axis(rotation);
                
                new_score = get_score(scored_atoms);
                deltaS = exp(-ipece_prm.beta * (new_score - score));
                if (ran2(idum) < deltaS) {
                    rotate(all_atoms,rotation);
                    score = new_score;
                    save_hist('r', rotation, score, &histories);
                }
                else {
                    recover(scored_atoms);
                    recover_axis();
                }
            }
            
            if (!fmod(i_iter,ipece_prm.axis_orientation_update)) {
                probe(all_atoms,ipece_prm.mem_radius);
                get_buried_cylinder_axis();
                if (!ipece_prm.inner_mem) {
                    scored_atoms = get_scored_atoms(all_atoms);
                }
                printf("iter=%i,score=%f\n",i_iter,score);
            }
        }
        /*
        START! Debug/Checking 
        
        for (i_hist=0;i_hist<histories.n;i_hist++) {
            if (histories.array[i_hist].motion_type == 't')
                printf("i_hist translate %6.3f %6.3f %6.3f\n", histories.array[i_hist].move[0],histories.array[i_hist].move[1],histories.array[i_hist].move[2]);
            else if (histories.array[i_hist].motion_type == 'r')
                printf("i_hist rotate    %6.3f %6.3f %6.3f %6.3f\n", histories.array[i_hist].move[0],histories.array[i_hist].move[1],histories.array[i_hist].move[2],histories.array[i_hist].move[3]);
            else
                printf("Error\n");
        }
        END!   Debug/Checking 
        */
        minimum_score_position(all_atoms,&histories,ipece_prm.n_iteration);
        add_mem_atoms(&all_atoms);
        write_probe("probe.pdb");
    }
	
	// Use the membrane specified in ipece.prm 
    else if (ipece_prm.add_mem && ipece_prm.spe_mem == 1) { 
		find_center_nonhet(all_atoms, central_point);
        for(i=0;i<3;i++) neg_central[i] = -central_point[i];
        
        translate(all_atoms, neg_central);
        save_hist('t', neg_central, ipece_prm.axis_score_weight, &histories);
		
		float lvector;
		lvector = sqrt(ipece_prm.mem_dirc[0] * ipece_prm.mem_dirc[0] 
					 + ipece_prm.mem_dirc[1] * ipece_prm.mem_dirc[1]
					 + ipece_prm.mem_dirc[2] * ipece_prm.mem_dirc[2]);
		if (lvector < 0.00001) printf("Error: length of specified membrane vector too small\n");
		rotation[0] = ipece_prm.mem_dirc[1] / lvector;
		rotation[1] = - ipece_prm.mem_dirc[0] / lvector;
		rotation[2] = 0.0;
		rotation[3] = acos(ipece_prm.mem_dirc[2]/lvector);
		
		rotate(all_atoms, rotation);
		save_hist('r', rotation, score, &histories);
		add_mem_atoms(&all_atoms);
	}
		
    if (ipece_prm.add_ion) {
        for (i_ion = 0; i_ion<ipece_prm.n_ion; i_ion++) {
            add_ion_atoms(&all_atoms,i_ion);
        }
    }

    if (ipece_prm.add_mem) {
        if (ipece_prm.move_back) {
            for (i_hist=histories.n-1; i_hist>=0; i_hist--) {
                back_move(all_atoms,histories.array[i_hist]);
            }
            histories.n = 0;
        }
    }
    
    pdb_out_fp = fopen(pdb_out_fname, "w");
    for (ia = 0; ia < all_atoms.n; ia++) {
        strncpy(sbuffer,all_atoms.array[ia].pdb_line,54);
        sprintf(sbuffer+30,"%8.3f%8.3f%8.3f",all_atoms.array[ia].r[0],all_atoms.array[ia].r[1],all_atoms.array[ia].r[2]);
        fprintf(pdb_out_fp, "%s%s",sbuffer,all_atoms.array[ia].pdb_line+54);
    }
    fclose(pdb_out_fp);
    
    /*     pdb_out_fp = fopen("debug.pdb", "w");
    for (ia = 0; ia < scored_atoms.n; ia++) {
        strncpy(sbuffer,scored_atoms.array[ia].pdb_line,54);
        sprintf(sbuffer+30,"%8.3f%8.3f%8.3f",scored_atoms.array[ia].r[0],scored_atoms.array[ia].r[1],scored_atoms.array[ia].r[2]);
        fprintf(pdb_out_fp, "%s%s",sbuffer,scored_atoms.array[ia].pdb_line+54);
    }
    fclose(pdb_out_fp);
    */
    return 0;
}


