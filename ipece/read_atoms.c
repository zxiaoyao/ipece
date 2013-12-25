#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

ATOMS read_atoms(char *pdb_fname) {
    ATOMS atoms;
    FILE *pdb_fp;
    char  line[TXT_LEN_MAX];
    char  sbuffer[TXT_LEN_MAX];
    int i;
    
    memset(&atoms, 0, sizeof(ATOMS));
    
    if ( (pdb_fp = fopen(pdb_fname, "r")) == NULL ) {
        printf("Can't open file %s\n", pdb_fname);
        exit(-1);
    }
    
    memset(line,0,TXT_LEN_MAX*sizeof(char));
    while (fgets(line, TXT_LEN_MAX, pdb_fp)) {
        if (strncmp("ATOM  ", line, 6) && strncmp("HETATM", line, 6)) continue;
        
        atoms.n++;
        atoms.array = realloc(atoms.array, atoms.n * sizeof(ATOM));
        if (!atoms.array) printf("Out of memory\n"), exit(-1);
        memset(&atoms.array[atoms.n-1],0,sizeof(ATOM));
        
        strcpy(atoms.array[atoms.n-1].pdb_line, line);
        for (i=0; i<3; i++) {
            strncpy(sbuffer, (line+30+8*i), 8); 
            sbuffer[8] = '\0'; 
            atoms.array[atoms.n-1].r[i] = atof(sbuffer);
        }
        memset(line,0,TXT_LEN_MAX*sizeof(char));
    }
    fclose(pdb_fp);
    return atoms;
}
