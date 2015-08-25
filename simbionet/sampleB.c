#include "sampleB.h"

#ifdef MDEBUG
	VETTOREi *_sampleB_p(VETTOREi *ris, const VETTOREi *x, int k, int replace, VETTOREd *p, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi *_sampleB_p(VETTOREi *ris, const VETTOREi *x, int k, int replace, VETTOREd *p)
#endif
{
	_Intestazione("\n***sampleB_p***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tx = ");
	_StampaRawVett_i(x);
	fprintf(fp_det, "\tsize =  %d\n", k);
	if (replace)
		fprintf(fp_det, "\treplace =  TRUE\n");
	else
		fprintf(fp_det, "\treplace =  FALSE\n");
	fprintf(fp_det, "\tprob = ");
	_StampaRawVett_d(p);
#endif

	if (x->dim == 1)
		ris = rep_v_i(ris, x, k);
	else
		ris = sample_p(ris, x, k, replace, p, "sampleB");

#ifdef DET
	fprintf(fp_det, "sampleB_p output:\n");
	fprintf(fp_det, "\tres = ");
	_StampaRawVett_i(ris);
#endif

	return ris;
}

#ifdef MDEBUG
	VETTOREi *_sampleB(VETTOREi *ris, const VETTOREi *x, int k, int replace, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi *_sampleB(VETTOREi *ris, const VETTOREi *x, int k, int replace)
#endif
{
	_Intestazione("\n***sampleB***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tx = ");
	_StampaRawVett_i(x);
	fprintf(fp_det, "\tsize =  %d\n", k);
	if (replace)
		fprintf(fp_det, "\treplace =  TRUE\n");
	else
		fprintf(fp_det, "\treplace =  FALSE\n");
#endif

	if (x->dim == 1)
		ris = rep_v_i(ris, x, k);
	else
		ris = sample(ris, x, k, replace, "sampleB");

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "sampleB output:\n");
	fprintf(fp_det, "\tres = ");
	_StampaRawVett_i(ris);
#endif

	return ris;
}
