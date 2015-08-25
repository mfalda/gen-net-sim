#include "mandel_calc.h"

// Callback function for parser errors
void OnError(muParserHandle_t hParser)
{
	char err[256];

	snprintf(err, 256, "Error:\n");
	snprintf(err, 256, "%s------\n", err);
	snprintf(err, 256, "%sMessage:  \"%s\"\n", err, mupGetErrorMsg(hParser));
	snprintf(err, 256, "%sToken:    \"%s\"\n", err, mupGetErrorToken(hParser));
	snprintf(err, 256, "%sPosition: %d\n", err, mupGetErrorPos(hParser));
	snprintf(err, 256, "%sErrc:     %d\n", err, mupGetErrorCode(hParser));

	error(err);
}

// param1 = xmin, xmax, ymin, ymax, passo, juliax, juliay
MATRICEi *mandel1(muParserHandle_t hParser_real, muParserHandle_t hParser_imag, MATRICEi *ris, VETTOREd *param, GString *reale, GString *img, double ris_x, double ris_y, int n)
{
	double c_imag, c_real, z_imag, z_real;
	double x, y, fVal;
	double args[3];
	int i, j = 1, k;
	GString *errore = NULL;

	_Intestazione("\n*** mandel1 ***\n");

	mupDefineVar(hParser_real, "x", &args[0]);
	mupDefineVar(hParser_real, "y", &args[1]);
	mupDefineVar(hParser_real, "z", &args[2]);

	mupDefineVar(hParser_imag, "x", &args[0]);
	mupDefineVar(hParser_imag, "y", &args[1]);
	mupDefineVar(hParser_imag, "z", &args[2]);

	for (y = ACCEDIv_d(param, 3); y <= ACCEDIv_d(param, 4); y += ris_y) {
		k = 1;
		for (x = ACCEDIv_d(param, 1); x <= ACCEDIv_d(param, 2); x +=ris_x) {
			z_imag = y + ACCEDIv_d(param, 6);
			z_real = x + ACCEDIv_d(param, 7);
			c_imag = z_imag;
			c_real = z_real;
			for (i = 0; i < n; i++) {
				args[0] = z_real;
				args[1] = z_imag;
				args[2] = c_real;
#ifdef FDEBUG
				fprintf(fp_fdbg, "Calcolo la funzione 'reale' in (%3.3f, %3.3f, %3.3f): ", args[0], args[1], args[2]);
#endif
				fVal = mupEval(hParser_real);
				if (!mupError(hParser_real))
					z_real = fVal;
#ifdef FDEBUG
				fprintf(fp_fdbg, "%3.3f\n", z_real);
				fprintf(fp_fdbg, "Calcolo la funzione 'imag' in (%3.3f, %3.3f, %3.3f): ", args[0], args[1], args[2]);
#endif
				args[2] = c_imag;
				fVal = mupEval(hParser_imag);
				if (!mupError(hParser_imag))
					z_imag = fVal;
#ifdef FDEBUG
				fprintf(fp_fdbg, "%3.3f\n", z_imag);
#endif
				if((z_real * z_real + z_imag * z_imag) > 4.0) {
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

SEXP mandel(SEXP param, SEXP real, SEXP img, SEXP n)
{
	VETTOREd *param1 = NULL;
	MATRICEi *ris1 = NULL;
	int nProtected = 0;
	double ris_x, ris_y;
	GString *real1, *img1;
	int n1, n_elem;
	muParserHandle_t hParser_real, hParser_imag;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** mandel ***\n");

	param1 = inVETTORE_d(param, &nProtected);
	real1 = inSTRINGA(real, &nProtected, "real");
	img1 = inSTRINGA(img, &nProtected, "imag");
	n1 = INTEGER_VALUE(n);
	ris_x = (ACCEDIv_d(param1, 2) - ACCEDIv_d(param1, 1)) / ACCEDIv_d(param1, 5);
	ris_y = (ACCEDIv_d(param1, 4) - ACCEDIv_d(param1, 3)) / ACCEDIv_d(param1, 5);
	n_elem = (int) ACCEDIv_d(param1, 5) + 1;
	CREAm_i(ris1, n_elem, n_elem);
	InitMatr_i(ris1, 0);
	hParser_real = mupCreate(muBASETYPE_FLOAT);
	hParser_imag = mupCreate(muBASETYPE_FLOAT);
	mupSetErrorHandler(hParser_real, OnError);
	mupSetErrorHandler(hParser_imag, OnError);
	mupSetExpr(hParser_real, real1->str);
	mupSetExpr(hParser_imag, img1->str);
	ris1 = mandel1(hParser_real, hParser_imag, ris1, param1, real1, img1, ris_x, ris_y, n1);
	mupRelease(hParser_real);
	mupRelease(hParser_imag);
	ris = daMATRICE_i(ris1, &nProtected);

	CANCELLAv_d(param1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
