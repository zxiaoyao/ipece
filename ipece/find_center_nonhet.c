#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

int find_center_nonhet(ATOMS atoms, float *central_point) {
    int iatom;
    int i;
    
    for (i=0;i<3;i++) central_point[i] = 0;
    
    for (iatom = 0; iatom<atoms.n; iatom++) {
        if ( !strncmp("ATOM  ", atoms.array[iatom].pdb_line, 6) ) {
            for (i=0;i<3;i++) central_point[i] += atoms.array[iatom].r[i];
        }
    }
    
    for (i=0;i<3;i++) central_point[i] = central_point[i]/(float)atoms.n;
    
    return 0;
}

