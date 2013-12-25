#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

void save_hist(char motion_type, float *move, float score, HISTORIES *histories_p) {
    
    int i;
    
    histories_p->n++;
    histories_p->array = realloc(histories_p->array, histories_p->n * sizeof(HISTORY));
    memset(&histories_p->array[histories_p->n-1], 0, sizeof(HISTORY));
    
    histories_p->array[histories_p->n-1].score = score;
    histories_p->array[histories_p->n-1].motion_type = motion_type;
    if (motion_type == 't') {
        for(i=0;i<3;i++) histories_p->array[histories_p->n-1].move[i] = move[i];
    }
    if (motion_type == 'r') {
        for(i=0;i<4;i++) histories_p->array[histories_p->n-1].move[i] = move[i];
    }
}

