/* 
NAME
        strip_comment - strip comment off from string

SYNOPSIS
        #include <mcce.h>
        STRING strip_comment(char *str);

DESCRIPTION
        strip comment off from input string.
        (comments begin with "REMARK" or "#" or "!")
        return the output string in a STRING structure.

SEE ALSO
        strip_spc

EXAMPLE
        char str[];
        STRING new_str;
        ...
        new_str = strip_comment(str);

AUTHOR
        Yifan Song, 05/27/2003 
*/      

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ipece.h"

STRING strip_comment(char *str) 
{
    
    int i;
    STRING      str_new;
    char        char_remark[] = "#!\n";
    char        str_remark[] = "REMARK";
    
    for (i = 0; i< strlen(str); i++) {
        if (i+6<strlen(str)) {
            if (!strncmp(str+i,str_remark,6)) {
                break;
            }
        }
        if (strchr(char_remark,str[i])) {
            break;
        }
    }
    strncpy(str_new.txt,str,i);
    str_new.txt[i] = '\0';
    
    return str_new;
}

