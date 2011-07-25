#ifndef CONFUSE_H
#define CONFUSE_H

#include <confuse.h>
#include "match.h"

int parse_confuse(char const *,Parameters *);
int write_confuse(char const *filename, 
        Parameters *what,
        SimInspiralTable default_siminsp,
        PPNParamStruc default_ppn);

#endif


