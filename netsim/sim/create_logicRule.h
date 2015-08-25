#ifndef CREATE_LOGICRULE_H
#define CREATE_LOGICRULE_H

#include "r_aux.h"
#include "globali.h"
#include "ext_parser.h"
#include "sample.h"

//typedef VETTOREd *(*FPRAND)(VETTOREd *norm_lev);

//VETTOREd *f(int n, VETTOREd *norm_lev);

VETTOREi *create_logicRule(VETTOREi *ris, const VETTOREi *v, int f_pr_and, muParserHandle_t hParser);

#endif
