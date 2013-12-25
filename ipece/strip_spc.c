/* 
NAME
        strip_spc - strip space and tab off from both sides of a string

SYNOPSIS
        #include <mcce.h>
        STRING  strip_spc(char *str);

DESCRIPTION
        strip space and tab off from both sides of input string.
        return the output string in a STRING structure.

SEE ALSO
        cut_comment

EXAMPLE
        char str[];
        STRING new_str;
        ...
        new_str = strip_spc(str);

AUTHOR
        Yifan Song, 05/27/2003 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

STRING  strip_spc(char *str) {
    
    STRING new_str;
    int len = strlen(str);
    int i;
    int j;
    char spc[] = " \t";
    
    memset(&new_str, 0, sizeof(STRING));
    
    for (i = 0; i < len; i++) {
        if ( !strchr(spc,str[i]) ) {
            break;
        }
    }
    
    for (j = strlen(str) - 1; j>=0; j--) {
        if ( !strchr(spc,str[j]) ) {
            break;
        }
    }
    strncpy(new_str.txt, str+i, (len-i));
    new_str.txt[j-i+1] = '\0';
    
    return new_str;
}

