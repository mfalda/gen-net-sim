#ifndef PROBMOD_UND_H
#define PROBMOD_UND_H

#include <math.h>

#include "r_aux.h"
#include "scoremodular.h"

struct RisProbM {
	double score;
	MATRICEi *conn_matr;
	MATRICEd *score_matr;
	char label[3];
	VETTOREi *indices;
};

MATRICEd *probmod2_und(MATRICEd *ris, MATRICEi *M, const VETTOREi *h, const VETTOREi *Sout, const VETTOREd *STout, const VETTOREd *Freq_out, const VETTOREd *toll, struct RisProbM *ris1);

#endif
