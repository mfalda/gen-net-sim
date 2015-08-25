#include "trpol2_R.h"

VETTOREd *trpol2(VETTOREd *ris, int n, double x)
{
	double mu = 10.0;
	double s;
	double pu = 0.0;
	int i, j;
	VETTOREd *pol = NULL;

	CREAv_d(pol, 100);
	for (i = 0; i < n; i++) {
		for (j = 1; j <= 100; j++) {
			mu = (mu + 2.0) / 2.0;
			ASSEGNAv_d(pol, j, mu);
			//~ pol[j] = mu = (mu + 2.0) / 2.0;
		}
		s = 0.0;
		for (j = 1; j <= 100; j++) {
			s = x * s + ACCEDIv_d(pol, j);
		}
		pu += s;
	}
	CANCELLAv_d(pol);
	CREAv_d(ris, 1);
	ASSEGNAv_d(ris, 1, pu);
	printf("%.16g\n", pu);
	return ris;
}

SEXP ret_trpol2(SEXP n, SEXP x)
{
	int nProtected = 0;
	int n1;
	VETTOREd *ris1 = NULL;
	double x1;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** ret_trpol2 ***\n");

	n1 = INTEGER_VALUE(n);
	x1 = NUMERIC_VALUE(x);
	ris1 = trpol2(ris1, n1, x1);

	ris = daVETTORE_d(ris1, &nProtected);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
