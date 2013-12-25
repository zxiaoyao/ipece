#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ipece.h"

IPECE_PRM read_ipece_prm() {
    IPECE_PRM ipece_prm;
    FILE *ipece_param_fp;
    STRING buffer_str;
    char buffer2[TXT_LEN_MAX];
    
    memset(&ipece_prm, 0, sizeof(IPECE_PRM));
    ipece_param_fp = fopen("ipece.prm","r");
    if (!ipece_param_fp) {
        printf ("Error: Input file \"ipece.prm\" not found.\n");
        exit(-1);
    }
    
    while (fgets(buffer_str.txt, TXT_LEN_MAX, ipece_param_fp)) {
        buffer_str = strip_comment(buffer_str.txt);
        
        if (!strncmp(buffer_str.txt, "LATTICE_SCALE", 13)) {
            ipece_prm.lattice_scale = atof(&buffer_str.txt[13]);
        }
        
        if (!strncmp(buffer_str.txt, "ADD_MEM", 7)) {
            buffer_str = strip_spc(&buffer_str.txt[7]);
            if (buffer_str.txt[0] == 't') ipece_prm.add_mem = 1;
        }
        if (!strncmp(buffer_str.txt, "INNER_MEM", 9)) {
            buffer_str = strip_spc(&buffer_str.txt[9]);
            if (buffer_str.txt[0] == 't') ipece_prm.inner_mem = 1;
        }
        if (!strncmp(buffer_str.txt, "MEM_THICKNESS", 13)) {
            ipece_prm.mem_thickness = atof(&buffer_str.txt[13]);
        }
        if (!strncmp(buffer_str.txt, "N_ITERATION", 11)) {
            ipece_prm.n_iteration = atoi(&buffer_str.txt[11]);
        }
        
        if (!strncmp(buffer_str.txt, "TRANSLATION_MAX", 15)) {
            ipece_prm.translation_max = atof(&buffer_str.txt[15]);
        }
        if (!strncmp(buffer_str.txt, "ROTATION_MAX", 12)) {
            ipece_prm.rotation_max = atof(&buffer_str.txt[12]);
        }
        if (!strncmp(buffer_str.txt, "BOUNDARY_EXTENTION_X", 20)) {
            ipece_prm.boundary_extention[0] = atof(&buffer_str.txt[20]);
        }
        if (!strncmp(buffer_str.txt, "BOUNDARY_EXTENTION_Y", 20)) {
            ipece_prm.boundary_extention[1] = atof(&buffer_str.txt[20]);
        }
        if (!strncmp(buffer_str.txt, "BOUNDARY_EXTENTION_Z", 20)) {
            ipece_prm.boundary_extention[2] = atof(&buffer_str.txt[20]);
        }
        if (!strncmp(buffer_str.txt, "BETA", 4)) {
            ipece_prm.beta = atof(&buffer_str.txt[4]);
        }
        if (!strncmp(buffer_str.txt, "SURFACE_EXP_RAD", 15)) {
            ipece_prm.surface_exp_rad = atof(&buffer_str.txt[15]);
        }
        if (!strncmp(buffer_str.txt, "MEM_NAME", 8)) {
            strncpy(ipece_prm.mem_name.txt, buffer_str.txt+32,14);
            ipece_prm.mem_name.txt[14]='\0';
        }
        if (!strncmp(buffer_str.txt, "MEM_RADIUS", 10)) {
            ipece_prm.mem_radius = atof(&buffer_str.txt[16]);
        }
        if (!strncmp(buffer_str.txt, "MEM_PLACE_SEP", 13)) {
            ipece_prm.mem_place_sep = atof(&buffer_str.txt[16]);
        }
        if (!strncmp(buffer_str.txt, "AXIS_SCORE_WEIGHT", 17)) {
            ipece_prm.axis_score_weight = atof(&buffer_str.txt[17]);
        }
        if (!strncmp(buffer_str.txt, "AXIS_ORIENTATION_UPDATE", 23)) {
            ipece_prm.axis_orientation_update = atoi(&buffer_str.txt[23]);
        }
        if (!strncmp(buffer_str.txt, "AXIS_EXTENTION", 14)) {
            ipece_prm.axis_extention = atof(&buffer_str.txt[14]);
        }
        if (!strncmp(buffer_str.txt, "MOVE_BACK", 9)) {
            buffer_str = strip_spc(&buffer_str.txt[9]);
            if (buffer_str.txt[0] == 't') ipece_prm.move_back = 1;
        }
        if (!strncmp(buffer_str.txt, "SCORED_ATOM", 11)) {
            ipece_prm.n_sa_type++;
            ipece_prm.sa_type = realloc(ipece_prm.sa_type, ipece_prm.n_sa_type*sizeof(SA_TYPE));
            memset(&ipece_prm.sa_type[ipece_prm.n_sa_type-1],0,sizeof(SA_TYPE));
            strncpy(ipece_prm.sa_type[ipece_prm.n_sa_type-1].res, buffer_str.txt+32, 3);
            strncpy(ipece_prm.sa_type[ipece_prm.n_sa_type-1].atom, buffer_str.txt+40, 4);
            strncpy(buffer2, buffer_str.txt+48, 8);
            buffer2[8]='\0';
            ipece_prm.sa_type[ipece_prm.n_sa_type-1].buried_score = atof(buffer2);
            
            strncpy(buffer2, buffer_str.txt+56, 8);
            buffer2[8]='\0';
            ipece_prm.sa_type[ipece_prm.n_sa_type-1].edge_score = atof(buffer2);
        }
        
        if (!strncmp(buffer_str.txt, "ADD_ION", 7)) {
            buffer_str = strip_spc(&buffer_str.txt[7]);
            if (buffer_str.txt[0] == 't') ipece_prm.add_ion = 1;
        }

        if (!strncmp(buffer_str.txt, "CAV_POS", 7)) {
            sscanf(buffer_str.txt,"CAV_POS %f %f %f",&ipece_prm.cav[0],&ipece_prm.cav[1],&ipece_prm.cav[2]);
        }

        if (!strncmp(buffer_str.txt, "CAV_THR", 7)) {
			ipece_prm.cav_thr = atof(&buffer_str.txt[7]);
        }

        /*
        ION_NAME                         O   HOH Z      1.2     0.4
        */
        if (!strncmp(buffer_str.txt, "ION_NAME", 8)) {
            ipece_prm.n_ion++;
            ipece_prm.ion_name = realloc(ipece_prm.ion_name,ipece_prm.n_ion*sizeof(STRING));
            ipece_prm.ion_radius = realloc(ipece_prm.ion_radius,ipece_prm.n_ion*sizeof(float));
            ipece_prm.ion_place_sep = realloc(ipece_prm.ion_place_sep,ipece_prm.n_ion*sizeof(float));
            
            strncpy(ipece_prm.ion_name[ipece_prm.n_ion-1].txt, buffer_str.txt+32,10);
            ipece_prm.ion_name[ipece_prm.n_ion-1].txt[10]='\0';
            
            strncpy(buffer2, buffer_str.txt+48, 8);
            buffer2[8]='\0';
            ipece_prm.ion_radius[ipece_prm.n_ion-1]=atof(buffer2);

            strncpy(buffer2, buffer_str.txt+56, 8);
            buffer2[8]='\0';
            ipece_prm.ion_place_sep[ipece_prm.n_ion-1]=atof(buffer2);
        }
		
		if (!strncmp(buffer_str.txt, "SPECIFY_MEM", 11)) {
			if (buffer_str.txt[14] == 't') {
				ipece_prm.spe_mem = 1;
				strncpy(buffer2, buffer_str.txt + 19, 10); buffer2[10] = '\0';
				//printf("buffer2: %s\n", buffer2);

				ipece_prm.mem_dirc[0] = atof(buffer2);
				strncpy(buffer2, buffer_str.txt + 29, 10); buffer2[10] = '\0';
				//printf("buffer2: %s\n", buffer2);

				ipece_prm.mem_dirc[1] = atof(buffer2);
				strncpy(buffer2, buffer_str.txt + 39, 10); buffer2[10] = '\0';
				//printf("buffer2: %s\n", buffer2);

				ipece_prm.mem_dirc[2] = atof(buffer2);
			}
			else {
				ipece_prm.spe_mem = 0;
				ipece_prm.mem_dirc[0] = 0.0;
				ipece_prm.mem_dirc[1] = 0.0;
				ipece_prm.mem_dirc[2] = 0.0;
			}
			printf("spe_mem: %d\n", ipece_prm.spe_mem);
			printf("mem_dirc: %.2f, %.2f, %.2f\n", ipece_prm.mem_dirc[0], ipece_prm.mem_dirc[1], ipece_prm.mem_dirc[2]);
		}
    }
    
    ipece_prm.half_mem_thickness = ipece_prm.mem_thickness / 2;
    ipece_prm.PI = atan(1.)*4.;
    return ipece_prm;
}
