#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

void translate(ATOMS atoms, float *translation) {
    int iatom;
    int i;
    
    for (iatom = 0; iatom<atoms.n; iatom++) {
        for(i=0;i<3;i++) atoms.array[iatom].r[i] += translation[i];
    }
}

