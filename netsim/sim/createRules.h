#ifndef CREATE_RULES_H
#define CREATE_RULES_H

#include "r_aux.h"
#include "globali.h"
#include "ext_parser.h"
#include "create_logicRule.h"

LISTA *createRules1(LISTA *ris, const MATRICEd *m, const GString *f_pr_and, int *len, muParserHandle_t hParser);

SEXP createRules(SEXP m, SEXP f_pr_and);

#endif
