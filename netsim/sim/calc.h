#ifndef CALC_H
#define CALC_H

#include <R.h>
#include <Rdefines.h>

#include "ext_parser.h"

SEXP calc(SEXP str, SEXP args);
SEXP calcv(SEXP str, SEXP args);

#endif
