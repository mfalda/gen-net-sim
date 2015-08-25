#ifndef MOD3_H
#define MOD3_H

#include "r_aux.h"
#include "globali.h"
#include "sampleB.h"
#include "scoremodular.h"
#include "hist.h"
#include "triangola.h"

// Cg = 0.0
MATRICEi *mod3(MATRICEi *ris, int Ng, int Ng_in, const VETTOREd *STin, const VETTOREd *STout, int max_con, double Cg);

#endif
