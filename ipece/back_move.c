#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

void back_move(ATOMS atoms, HISTORY history) 
{
    float move[4];
    int i;
    
    memcpy(move, history.move, 4*sizeof(float));
    if (history.motion_type == 't') {
        for(i=0;i<3;i++)
            move[i] = -move[i];
        translate(atoms,move);
    }
    if (history.motion_type == 'r') {
        move[3] = -move[3];
        rotate(atoms,move);
    }
}
