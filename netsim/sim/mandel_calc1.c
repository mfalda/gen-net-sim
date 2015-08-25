#include "mandel_calc1.h"

// param1 = xmin, xmax, ymin, ymax, passo, juliax, juliay
MATRICEi *mandel1(MATRICEi *ris, VETTOREd *param, GString *reale, GString *img, double ris_x, double ris_y, int n)
{
	double c_imag, c_real, z_imag, z_real;
	double x, y;
	double args[3];
	int i, j = 1, k;
	GString *errore = NULL;

	_Intestazione("\n*** mandel1 ***\n");

#ifdef FDEBUG
	fprintf(fp_fdbg, "Controllo la funzione 'reale' definita come '%s' ... ", reale->str);
#endif
	errore = DefF("reale", reale->str, 3, 0);
	if (errore != NULL) {
		error(errore->str);
		g_string_free(errore, TRUE);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "ok\n");
	fprintf(fp_fdbg, "Controllo la funzione 'imag' definita come '%s' ... ", img->str);
#endif
	errore = DefF("imag", img->str, 3, 0);
	if (errore != NULL) {
		error(errore->str);
		g_string_free(errore, TRUE);
	}
#ifdef FDEBUG
	fprintf(fp_fdbg, "ok\n");
#endif
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
				errore = Calcola("reale", args, &z_real);
				if (errore != NULL) {
					error(errore->str);
					g_string_free(errore, TRUE);
				}
#ifdef FDEBUG
				fprintf(fp_fdbg, "%3.3f\n", z_real);
				fprintf(fp_fdbg, "Calcolo la funzione 'imag' in (%3.3f, %3.3f, %3.3f): ", args[0], args[1], args[2]);
#endif
				args[2] = c_imag;
				errore = Calcola("imag", args, &z_imag);
				if (errore != NULL) {
					error(errore->str);
					g_string_free(errore, TRUE);
				}
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
	Ripristina();
	ris1 = mandel1(ris1, param1, real1, img1, ris_x, ris_y, n1);
	Cancella();
	ris = daMATRICE_i(ris1, &nProtected);

	CANCELLAv_d(param1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
