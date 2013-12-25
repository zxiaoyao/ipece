#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

void minimum_score_position(ATOMS atoms,HISTORIES *histories_p, int n_iteration) 
{
    int i_hist;
    int min_i_hist;
    float min_score;
    
    printf("Total number of steps %d\n",histories_p->n);
    if (histories_p->n > 1) {
        min_score = histories_p->array[1].score;
        min_i_hist = 1;
        for (i_hist=2; i_hist<histories_p->n; i_hist++) {
            if (histories_p->array[i_hist].score < min_score) {
                min_score = histories_p->array[i_hist].score;
                min_i_hist = i_hist;
            }
        }
    }
    else {
        min_score = histories_p->array[0].score;
        min_i_hist = 0;
    }
    if (n_iteration == 0)
        min_i_hist = -1;

    printf("Minimun score position is at step %d, and the score is %f\n",min_i_hist,min_score);
	
    for (i_hist=histories_p->n-1; i_hist>min_i_hist; i_hist--) {
        back_move(atoms,histories_p->array[i_hist]);
    }
    histories_p->n = min_i_hist+1;
}

