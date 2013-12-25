#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

void read_acc_atom(ATOMS atoms, char *acc_fname) 
{
    FILE *acc_fp;
    char  line[TXT_LEN_MAX];
    int ka;
    int ka_start;
    
    if ( (acc_fp = fopen(acc_fname, "r")) == NULL ) {
        printf("Can't open file %s\n", acc_fname);
        exit(-1);
    }
    
    /* 
    .0123456789 123456789 123456789 12345
    . atm   res no.           area
    . N     THR A   5          7.674
    . CA    THR A   5         22.905
    . C     THR A   5          2.004
    . O     THR A   5         14.921
    . NH1   ARG C 151  43.413
    */
    ka = 0;
    //fgets(line, TXT_LEN_MAX, acc_fp); /* Skip the first line of acc.atm. */
    memset(line, 0, TXT_LEN_MAX*sizeof(char));
    while (fgets(line, TXT_LEN_MAX, acc_fp)) {
        ka_start = ka;
		char buff_line[TXT_LEN_MAX], buff_atom[TXT_LEN_MAX];
		strncpy(buff_line, line+17, 4); buff_line[4] = '\0';
		strncpy(buff_atom, atoms.array[ka].pdb_line+22, 4); buff_atom[4] = '\0';
		
        while ( strncmp(line+7, atoms.array[ka].pdb_line+12, 4)   // atom name 
			 || strncmp(line+12, atoms.array[ka].pdb_line+17, 3)   // residue name
			 || strncmp(line+16, atoms.array[ka].pdb_line+21, 1)  // chain ID
			 || (atoi(buff_line) != atoi(buff_line))) { /* Keep within the loop if not same residue name or not same atom name. */
            ka++;
            if (ka == atoms.n) ka = 0; /* Search hit bottom, continue from top */
            if (ka == ka_start) { /* Loop back to beginning, no match found */
                ka_start = -1; 
                break;
            }
        }
        if (ka_start != -1) {
            atoms.array[ka].acc = atof(line+21);
            ka++;
            if (ka == atoms.n) ka = 0; /* Search hit bottom, continue from top */
        }
        else {
            printf("Warning! one line in acc.atm not recognized:\n%s\nFile acc.atm needs to be calculated from same input PDB file.\n",line);
        }
    }
    /*
    for (ka=0;ka<atoms.n;ka++) {
        printf("%8.3f,%s\n",atoms.array[ka].acc,atoms.array[ka].pdb_line);
    }
    */
}
            
