#include "lsoda_oscillatore.h"

// oscillatore: m * x'' + b * x' + k * x = F cos(o * t)
int funzione(double t, const double y[], double f[],
      void *params)
{
	double m = ((Params *) params)->m;
	double b = ((Params *) params)->b;
	double k = ((Params *) params)->k;
	double F = ((Params *) params)->F;
	double o = ((Params *) params)->o;

	f[0] = y[1];
	f[1] = -(k / m) * y[0] - (b / m) * y[1] + F * cos(o * t);
	return GSL_SUCCESS;
}

int jacobiano (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
	double m = ((Params *) params)->m;
	double b = ((Params *) params)->b;
	double k = ((Params *) params)->k;
	double F = ((Params *) params)->F;
	double o = ((Params *) params)->o;

	gsl_matrix_view dfdy_mat
	= gsl_matrix_view_array (dfdy, 2, 2);
	gsl_matrix * matr = &dfdy_mat.matrix;
	gsl_matrix_set (matr, 0, 0, 0.0);
	gsl_matrix_set (matr, 0, 1, 1.0);
	gsl_matrix_set (matr, 1, 0, -k / m);
	gsl_matrix_set (matr, 1, 1, -b / m);
	dfdt[0] = 0.0;
	dfdt[1] = -o * F * sin(o * t);
	return GSL_SUCCESS;
}

MATRICEd *lsoda_oscillatore1(MATRICEd *ris, LISTA *params, VETTOREd *X0, VETTOREd *times, const char *metodo, double atol, double rtol, double stat_thr, double stat_width)
{
	gsl_odeiv_step_type *T =NULL;
	int i, j;
	MATRICEd *tmp = NULL;
	Params p;

	_Intestazione("\n*** lsoda_oscillatore1 ***\n");

	// M<-parms[[1]]
	CtrlLlst(params, 1);
	p.m = ACCEDIv_d(ACCEDIlst(params, 1, vd), 1);
	// R<-parms[[2]]
	// da parametro
	// N<-parms[[3]]
	CtrlLlst(params, 2);
	p.b = ACCEDIv_d(ACCEDIlst(params, 2, vd), 1);
	// k<-parms[[4]]
	CtrlLlst(params, 3);
	p.k = ACCEDIv_d(ACCEDIlst(params, 3, vd), 1);
	// alpha<-parms[[5]]
	CtrlLlst(params, 4);
	p.F = ACCEDIv_d(ACCEDIlst(params, 4, vd), 1);
	// theta<-parms[[6]]
	CtrlLlst(params, 5);
	p.o = ACCEDIv_d(ACCEDIlst(params, 5, vd), 1);

	if (!strncmp(metodo, "rkf45", 5))
		T = gsl_odeiv_step_rkf45;
	else if (!strncmp(metodo, "rkck", 4))
		T = gsl_odeiv_step_rkck;
	else if (!strncmp(metodo, "rk8pd", 5))
		T = gsl_odeiv_step_rk8pd;
	else if (!strncmp(metodo, "rk2imp", 6))
		T = gsl_odeiv_step_rk2imp;
	else if (!strncmp(metodo, "rk4imp", 6))
		T = gsl_odeiv_step_rk4imp;
	else if (!strncmp(metodo, "gear1", 5))
		T = gsl_odeiv_step_gear1;
	else if (!strncmp(metodo, "gear2", 5))
		T = gsl_odeiv_step_gear2;
	else
		error("The requested stepping function does not exist!\n");

	gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 2);
	gsl_odeiv_control *c = gsl_odeiv_control_y_new (atol, rtol);
	gsl_odeiv_evolve *e	= gsl_odeiv_evolve_alloc (2);

	gsl_odeiv_system sys = {funzione, NULL, 2, &p};

	double t = 0.0;
	double h = 1e-6;
	double y[2] = { ACCEDIv_d(X0, 1), ACCEDIv_d(X0, 2) };

	CREAm_d(tmp, LENGTHv_d(times), 3);
	int st = 0;
	double y1 = 0.0;
	double max_stat;
	if (stat_width == 0.0)
		max_stat = LENGTHv_d(times);
	else
		max_stat = (double) stat_width * LENGTHv_d(times);
	for (i = 1; (!stat_width || st < max_stat) && i <= LENGTHv_d(times); i++) {
		double ti = ACCEDIv_d(times, i);

		while (t < ti) {
			gsl_odeiv_evolve_apply (e, c, s,
			                        &sys,
			                        &t, ti, &h,
			                        y);
		}

		if (fabs(y[0] - y1) < stat_thr)
			st++;
	  ASSEGNAm_d(tmp, i, 1, t);
	  ASSEGNAm_d(tmp, i, 2, y[0]);
	  ASSEGNAm_d(tmp, i, 3, y[1]);
		y1 = y[0];
	}
	CREAm_d(ris, i - 1, 3);
	for (j = 1; j < i; j++)
		copia1_m_riga_d(ris, j, tmp, j);
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	StrBilanciam();

   return ris;
}

SEXP lsoda_oscillatore(SEXP parms, SEXP X0, SEXP times, SEXP metodo, SEXP atol, SEXP rtol, SEXP stat_thr, SEXP stat_width)
{
	int i, ll, nProtected = 0;
	VETTOREd *X01 = NULL, *times1 = NULL;
	MATRICEd *ris1 = NULL;
	LISTA *parms1 = NULL;
	GString *metodo1 = NULL;
	double atol1, rtol1, stat_thr1, stat_width1;
#ifdef MDEBUG
	GString **nomi;
	char tmp[15];
#endif
	enum TIPO *tipi;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** lsoda_oscillatore ***\n");

	ll = length(parms);
	tipi = mia_alloc(ll, enum TIPO);
	if (ll > 0 && tipi == NULL) {
		Rprintf("Not enough memory (lsoda_oscillatore # %d, tipi)", __LINE__ - 2);
		error("");
	}
#ifdef MDEBUG
	nomi = mia_alloc(ll, GString *);
	if (ll > 0 && nomi == NULL) {
		Rprintf("Not enough memory (lsoda_oscillatore # %d, nomi)", __LINE__ - 2);
		error("");
	}
#endif
	X01 = inVETTORE_d(X0, &nProtected);
	times1 = inVETTORE_d(times, &nProtected);
	metodo1 = inSTRINGA(metodo, &nProtected, "metodo");
	atol1 = NUMERIC_VALUE(atol);
	rtol1 = NUMERIC_VALUE(rtol);
	stat_thr1 = NUMERIC_VALUE(stat_thr);
	stat_width1 = NUMERIC_VALUE(stat_width);
	tipi[0] = VETTd;
	tipi[1] = VETTd;
	tipi[2] = VETTd;
	tipi[3] = VETTd;
	tipi[4] = VETTd;
	for (i = 0; i < ll; i++) {
#ifdef MDEBUG
		if (i < ll)
			snprintf(tmp, 15, "Lista parms %d", i + 1);
		else
			tmp[0] = '\0';
		CREAstr(nomi[i], tmp);
#endif
	}
	parms1 = inLISTA(parms, &nProtected, ll, tipi, nomi);
	libera(tipi);
#ifdef MDEBUG
	for (i = 0; i < ll; i++)
		CANCELLAstr(nomi[i]);
	libera(nomi);
#endif

	ris1 = lsoda_oscillatore1(ris1, parms1, X01, times1, metodo1->str, atol1, rtol1, stat_thr1, stat_width1);
	ris = daMATRICE_d(ris1, &nProtected);

	CancellaLISTA(parms1, true);

	CANCELLAv_d(X01);
	CANCELLAv_d(times1);
	CANCELLAstr(metodo1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}

