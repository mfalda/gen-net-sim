#include "mod1.h"

#define g_ind globali.mod1.ind
#define g_tmp1_i globali.mod1.tmp1_i
#define g_x globali.mod1.x
#define g_s globali.mod1.s

MATRICEi *mod1(MATRICEi *ris, int Ng, double Cg, int max_con)
{
	_Intestazione("\n***mod1***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tNg =  %d\n", Ng);
	fprintf(fp_det, "\tCg =  %.16g\n", Cg);
	fprintf(fp_det, "\tmax_con =  %d\n", max_con);
#endif

	CREAm_i(ris, Ng, Ng);
	InitMatr_i(ris, 0);
	if (Ng > 1) {
		// g_x <- seq(1,Ng-1,1)
		g_x = seq_i(g_x, 1, Ng - 1, 1);
		// g_ind<-cbind(g_x,g_x+1)
		g_tmp1_i = somma_vs_i(g_tmp1_i, g_x, 1);
		g_ind = cbind2v_i(g_ind, g_x, g_tmp1_i);
		// m[g_ind] <- 1
		assegna1_ms_indx2_i(ris, g_ind, 1);
		// m[Ng][1] = 1;
		ASSEGNAm_i(ris, Ng, 1, 1);
		g_tmp1_i = vettore2s_i(g_tmp1_i, 1, -1);
		g_s = sampleB(g_s, g_tmp1_i, 1, 0);
		if (ACCEDIv_i(g_s, 1) == 1) {
			// m[Ng][1] = 0;
			ASSEGNAm_i(ris, Ng, 1, 0);
			// m[1][Ng] = 1;
			ASSEGNAm_i(ris, 1, Ng, 1);
		}
	}
	else
		// m[1][1] = 1;
		ASSEGNAm_i(ris, 1, 1, 1);
	if (Cg > 0.0 && Ng > 2)
		ris = triangola(ris, Cg, max_con);

	//~ CANCELLAm_i(g_ind);
	//~ CANCELLAv_i(g_tmp1_i);
	//~ CANCELLAv_i(g_x);
	//~ CANCELLAv_i(g_s);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "mod1 output:\n");
	fprintf(fp_det, "\tm = ");
	_StampaRawMatr_i(ris);
#endif

	return ris;
}
