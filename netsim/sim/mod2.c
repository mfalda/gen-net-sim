#include "mod2.h"

MATRICEi *mod2(MATRICEi *ris, int Ng, double Cg, int max_con)
{
	int r;

	_Intestazione("\n***mod2***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tNg =  %d\n", Ng);
	fprintf(fp_det, "\tCg =  %.16g\n", Cg);
	fprintf(fp_det, "\tmax_con =  %d\n", max_con);
#endif

	CREAm_i(ris, Ng, Ng);
	InitMatr_i(ris, 0);
	// m[2:Ng,1]<-1
	for (r = 2; r <= Ng; r++)
		ASSEGNAm_i(ris, r, 1, 1);

	if (Cg > 0.0 && Ng > 2)
		//m <- triangola(M=m,Cg=Cg,max_con)
		ris = triangola(ris, Cg, max_con);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "mod2 output:\n");
	fprintf(fp_det, "\tm = ");
	_StampaRawMatr_i(ris);
#endif

	return ris;
}

