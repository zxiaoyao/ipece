#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

float get_rad(ATOM atom) {
    char elem[3];
    strncpy(elem, atom.pdb_line+12, 2);
    elem[2] = '\0';
    
    if (!strcmp(elem, " C")) return 1.6;
    else if (!strcmp(elem, " N")) return 1.4;
    else if (elem[1]=='H') return 1.0;
    else if (!strcmp(elem, " O")) return 1.2;
    else return 1.6;
}
