#ifndef WRITETABLE_H
#define WRITETABLE_H

#include "r_aux.h"

// weights: write.table_m_i(M*Mneg, file = etich,quote=FALSE,sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = as.character(seq(1,N,1)),col.names = NA, qmethod = c("escape", "double"))
// SIMdata: write.table_m_d(D, file = etich,quote=FALSE,sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = as.character(seq(1,N,1)),col.names = NA, qmethod = c("escape", "double"))
// Rules: write.table_m_i(REG, file = etich,quote=FALSE,sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = as.character(seq(1,N,1)),col.names = NA, qmethod = c("escape", "double"))
// parameters: write.table_v3_d(cbind(lambda,alpha,beta), file = etich,quote=FALSE,sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = as.character(seq(1,N,1)),col.names = NA, qmethod = c("escape", "double"))
// Rules_simulateprofiles: write.table_m_i(REG, file = etich,quote=FALSE,sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = as.character(seq(1,N,1)),col.names = NA, qmethod = c("escape", "double"))

void write_m_i(const char *nome, MATRICEi *m);
SEXP write_m_int(SEXP nome, SEXP m);

void write_m_d(const char *nome, MATRICEd *m);
SEXP write_m_double(SEXP nome, SEXP x);

//~ #ifdef MDEBUG
	//~ #define write_mn_d(A, B, C) _write_m_d(A, B, C, #B, __FILE__, __LINE__)
//~ #else
	//~ #define write_mn_d(A, B, C) _write_m_d(A, B, C)
//~ #endif
//~ #ifdef MDEBUG
	//~ void _write_mn_d(const char *nome, MATRICEd *m, const LISTA *nomi);
//~ #else
	//~ void _write_mn_d(const char *nome, MATRICEd *m, const LISTA *nomi);
//~ #endif

//~ SEXP write_param_double(SEXP nome, SEXP m, SEXP nomi);

void write_vn_d(const char *nomefile, VETTOREd *v, const char *nome);
SEXP write_vparam_double(SEXP nomefile, SEXP m, SEXP nome);

void write_mn_d(const char *nomefile, MATRICEd *m, const char *nome);
SEXP write_mparam_double(SEXP nomefile, SEXP m, SEXP nome);

#endif
