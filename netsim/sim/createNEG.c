#include "createNEG.h"

#define g_tmp_i globali.createNEG.tmp_i
#define g_segno globali.createNEG.segno
#define g_ind globali.createNEG.ind

MATRICEi *createNEG1(MATRICEi *ris, const MATRICEi *m)
{
	int i, l;

	_Intestazione("\n***createNEG***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tM = ");
	_StampaRawMatr_i(m);
#endif

	CREAm_i(ris, LENGTHm1_i(m), LENGTHm2_i(m));
	InitMatr_i(ris, 0);
	// g_ind<-which(M==1,arr.g_ind=TRUE)
	g_ind = which_m_indxeq_i(g_ind, m, 1);
	l = LENGTHv_i(g_ind);
	// g_segno<-sample(c(-1,1),L,replace=TRUE)
	g_tmp_i = vettore2s_i(g_tmp_i, -1, 1);
	g_segno = sample(g_segno, g_tmp_i, l, 1, "createNEG");
	// M[g_ind]<-segno
	for (i = 1; i <= l; i++)
		ASSEGNAmv_i(ris, ACCEDIv_i(g_ind, i), ACCEDIv_i(g_segno, i));
	//~ CANCELLAv_i(g_tmp_i);
	//~ CANCELLAv_i(g_segno);
	//~ CANCELLAv_i(g_ind);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "createNEG output:\n\tM = ");
	_StampaRawMatr_i(ris);
#endif

	return ris;
}

SEXP createNEG(SEXP m)
{
	int nProtected = 0;
	MATRICEi *m1, *ris1 = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** createNEG ***\n");

	m1 = inMATRICE_i(m, &nProtected);

	ris1 = createNEG1(ris1, m1);

	ris = daMATRICE_i(ris1, &nProtected);

	CANCELLAm_i(m1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
