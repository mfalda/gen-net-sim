#ifndef SIMULATEPROFILES_H
#define SIMULATEPROFILES_H

#ifdef WIN32
	#include <io.h>
#else
	#include <unistd.h>
#endif
#include <dirent.h>  // per stat()
#include <sys/types.h>  // per stat()
#include <sys/stat.h>   // per stat()

#include "r_aux.h"
#include "globali.h"
#include "createRules.h"
#include "write_table.h"
#include "read_table.h"
#include "dinamica.h"
#include "lsoda.h"

enum TParam { ESTERNO, UNIF, NORM, LOG_NORM, SN, ESTERNO_UNIF, ESTERNO_NORM, ESTERNO_LOG_NORM, SN_UNIF, SN_NORM, SN_LOG_NORM };

enum NParam { LAMBDA, ALPHA, BETA, XMIN, XMAX, XZERO };

#define MAXPAR 6

MATRICEd *simulateprofiles1(MATRICEd *ris, MATRICEd *weights, LISTA *regole, GString *f_pr_and, VETTOREd *X0, VETTOREd *Xmin, VETTOREd *Xmax, VETTOREd *lambda, GString *act_fun, VETTOREd *alpha, VETTOREd *beta, double sd_noise, VETTOREd *times, GString *method, MATRICEd *ext_in, LISTA *ext_fun, int save, int rete, int ind_itera, int ko, int pad_reti, int pad_exp, int pad_ko, VETTOREd *ext1, VETTOREd *ext2, muParserHandle_t *hParsers, double stat_thr, double stat_width, int *stat);

SEXP simulateprofiles(SEXP N, SEXP weights, SEXP Rules, SEXP f_pr_and, SEXP act_fun, SEXP alpha, SEXP beta, SEXP lambda, SEXP Xmin, SEXP Xmax, SEXP X0, SEXP param, SEXP ko_experim, SEXP sd_noise, SEXP times, SEXP stat_thr, SEXP stat_width, SEXP method, SEXP ext_in, SEXP ext_fun, SEXP num_exp, SEXP save);

#endif
