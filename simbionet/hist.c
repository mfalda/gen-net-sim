#include "hist.h"

// Logica di R (binning.c)

#ifdef MDEBUG
	VETTOREi *_hist1(VETTOREi *ris, VETTOREi *x, VETTOREi *breaks, int right, int include_border, int naok, const char *nome, const char *nomefile, int linea)
#else
	VETTOREi *_hist1(VETTOREi *ris, VETTOREi *x, VETTOREi *breaks, int right, int include_border, int naok)
#endif
{
	int i, lo, hi;
	int n, nb1, new1, lft;

	_Intestazione("\n*** hist1 ***\n");

	assert(x != NULL && breaks != NULL);
#ifdef FDEBUG
	_StampaVett_i(x);
	_StampaVett_i(breaks);
	fprintf(fp_fdbg, "right = %d\n", right);
	fprintf(fp_fdbg, "include_border = %d\n", include_border);
	fprintf(fp_fdbg, "naok = %d\n", naok);
#endif
	n = LENGTHv_i(x);
	nb1 = LENGTHv_i(breaks);
	lft = !right;

	CREAv_i(ris, nb1 - 1);
	InitVett_i(ris, 0);

	for (i = 1 ; i <= n ; i++) {
		if (R_FINITE(ACCEDIv_i(x, i))) {
			lo = 1;
			hi = nb1;
			if (ACCEDIv_i(breaks, lo) <= ACCEDIv_i(x, i) &&
			        (ACCEDIv_i(x, i) < ACCEDIv_i(breaks, hi) ||
			         (ACCEDIv_i(x, i) == ACCEDIv_i(breaks, hi) && include_border))) {
				while (hi - lo >= 2) {
					new1 = (int) (hi + lo) / 2;
					if (ACCEDIv_i(x, i) > ACCEDIv_i(breaks, new1) || (lft && ACCEDIv_i(x, i) == ACCEDIv_i(breaks, new1)))
						lo = new1;
					else
						hi = new1;
				}
				ASSEGNAv_i(ris, lo, ACCEDIv_i(ris, lo) + 1);
			}
		}
		else if (!naok)
			error("NA's in \"hist\", ... NAOK=FALSE)");
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	_StampaVett_i(ris);
#endif

	StrBilanciam();

	return ris;
}

SEXP hist(SEXP x, SEXP breaks, SEXP right, SEXP include_border, SEXP naok)
{
	int nProtected = 0;
	int right1, include_border1, naok1;
	VETTOREi *x1, *breaks1;
	VETTOREi *ris1 = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** hist ***\n");

	x1 = inVETTORE_i(x, &nProtected);
	breaks1 = inVETTORE_i(breaks, &nProtected);
	right1 = INTEGER_VALUE(right);
	include_border1 = INTEGER_VALUE(include_border);
	naok1 = INTEGER_VALUE(naok);

	ris1 = hist1(ris1, x1, breaks1, right1, include_border1, naok1);
	ris = daVETTORE_i(ris1, &nProtected);

	CANCELLAv_i(x1);
	CANCELLAv_i(breaks1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
