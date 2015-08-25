#include "mandel1.h"

// param1 = xmin, xmax, ymin, ymax, passo, juliax, juliay
MATRICEi *mandel1(MATRICEi *ris, VETTOREd *param, double ris_x, double ris_y, int n)
{
	_Complex double c, z;
	double x, y;
	int i, j = 1, k;

	_Intestazione("\n*** mandel1 ***\n");

	for (y = ACCEDIv_d(param, 3); y <= ACCEDIv_d(param, 4); y += ris_y) {
		k = 1;
		for (x = ACCEDIv_d(param, 1); x <= ACCEDIv_d(param, 2); x +=ris_x) {
			z = (x + ACCEDIv_d(param, 6)) + (y + ACCEDIv_d(param, 7)) * 1i;
			c = z;
			for (i = 0; i < n; i++) {
				z = z * z + c;
				if(creal(z) * creal(z) + cimag(z) * cimag(z) > 4.0) {
					ASSEGNAm_i(ris, k, j, i);
					break;
				}
			}
			k++;
		}
		j++;
	}

	StrBilanciam();

	return ris;
}

SEXP mandel(SEXP param, SEXP n)
{
	VETTOREd *param1 = NULL;
	MATRICEi *ris1 = NULL;
	int nProtected = 0;
	double ris_x, ris_y;
	int n1, n_elem;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** mandel ***\n");

	param1 = inVETTORE_d(param, &nProtected);
	n1 = INTEGER_VALUE(n);
	ris_x = (ACCEDIv_d(param1, 2) - ACCEDIv_d(param1, 1)) / ACCEDIv_d(param1, 5);
	ris_y = (ACCEDIv_d(param1, 4) - ACCEDIv_d(param1, 3)) / ACCEDIv_d(param1, 5);
	n_elem = (int) ACCEDIv_d(param1, 5) + 1;
	CREAm_i(ris1, n_elem, n_elem);
	InitMatr_i(ris1, 0);
	ris1 = mandel1(ris1, param1, ris_x, ris_y, n1);
	ris = daMATRICE_i(ris1, &nProtected);

	CANCELLAv_d(param1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
