#ifndef READTABLE_H
#define READTABLE_H

#include "r_aux.h"

#define MAX_BUF 256 // lunghezza massima del nome del parametro

LISTA *read_lst_i(LISTA *ris, int n, const char *nomefile);

MATRICEi *read_m_i(MATRICEi *ris, const char *nomefile);

MATRICEd *read_m_d(MATRICEd *ris, int n, const char *nomefile);

VETTOREd *read_vn_d(VETTOREd *ris, const char *nomefile);

#endif
