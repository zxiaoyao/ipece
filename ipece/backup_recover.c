#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

extern BOX box;

void backup(ATOMS atoms) {
    int ia;
    int i;
    
    for (ia=0; ia<atoms.n; ia++) {
        for(i=0;i<3;i++)
            atoms.array[ia].backup_r[i] = atoms.array[ia].r[i];
    }
}

void recover(ATOMS atoms) {
    int ia;
    int i;
    
    for (ia=0; ia<atoms.n; ia++) {
        for(i=0;i<3;i++)
            atoms.array[ia].r[i] = atoms.array[ia].backup_r[i];
    }
}

void backup_axis() {
    int i;
    
    for(i=0;i<3;i++)
        box.backup_axis[i] = box.buried_cylinder_axis[i];
}

void recover_axis() {
    int i;

    for(i=0;i<3;i++)
        box.buried_cylinder_axis[i] = box.backup_axis[i];
}

