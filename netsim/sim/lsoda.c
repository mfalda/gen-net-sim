#include "lsoda.h"

#define g_tipi globali.tipi
#define g_ntipi globali.n_tipi

#define g_y globali.lsoda.y
#define g_y1 globali.lsoda.y1
#define g_y_prec globali.lsoda.y_prec
#define g_targetT globali.lsoda.targetT
#define g_targ globali.lsoda.targ
#define g_tmp_ris globali.lsoda.tmp_ris

int funzione(double t, const double *y, double *f,
             void *params)
{
	int i, dim = ((Params *) params)->dim;

	VETTOREd *alpha = ((Params *) params)->alpha;
	VETTOREd *beta = ((Params *) params)->theta;
	VETTOREd *k = ((Params *) params)->k;
	LISTA *R = ((Params *) params)->r;
	GString *af = ((Params *) params)->af;
	MATRICEd* M = ((Params *) params)->M;
	MATRICEi *N = ((Params *) params)->N;
	MATRICEd* ext_in = ((Params *) params)->ei;
	LISTA *ext_fun = ((Params *) params)->ef;
	VETTOREd *xmin = ((Params *) params)->xmin;
	double sd_noise = ((Params *) params)->sd;
	muParserHandle_t *hParsers = ((Params *) params)->hp;

	CREAv_d(g_y1, dim);
	for (i = 1; i <= dim; i++)
		ASSEGNAv_d(g_y1, i, y[i - 1]);
	// g_targetT<-target(y,R,M,N,EXT.IN,EXT.FUN,t)
	g_targ = target1(g_targ, g_y1, R, M, N, ext_in, ext_fun, t, hParsers, sd_noise);
	if (!strncmp(af->str, "linear", 6))
		// g_targetT <-  g_targetT * (1 - xmin) + xmin
		g_targetT = f_aux8_d(g_targetT, g_targ, xmin);
	else
		// g_targetT <- (1 / (1 + exp( -alpha*(g_targ - theta)))) * (1 - xmin) + xmin
		g_targetT = f_aux5_d(g_targetT, alpha, g_targ, beta, xmin);

	//func < -list(k*(g_targetT - y))
	for (i = 1; i <= dim; i++) {
		f[i - 1] = ACCEDIv_d(k, i) * (ACCEDIv_d(g_targetT, i) - y[i - 1]);
	}
	//~ CANCELLAv_d(g_y1);
	//~ CANCELLAv_d(g_targetT);
	//~ CANCELLAv_d(g_targ);

	return GSL_SUCCESS;
}

MATRICEd *lsoda1(MATRICEd *ris, LISTA *parms, LISTA *R, LISTA *ext_fun, VETTOREd *X0, VETTOREd *times,const char *metodo, double atol, double rtol, muParserHandle_t *hParsers, double stat_thr, double stat_width, int *stat, double sd_noise)
{
	gsl_odeiv_step_type *T = NULL;
	int i, j;
	double ti;
	double t = 0.0;
	double h = 1e-6;
	int st = 0, st_all;
	double y1 = 0.0;
	double max_stat;
	Params p;

	_Intestazione("\n***lsoda***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tX0 = ");
	_StampaRawVett_d(X0);
	fprintf(fp_det, "\ttimes = ");
	_StampaRawVett_d(times);
#endif

	// M<-parms[[1]]
	CtrlLlst(parms, 1);
	p.M = ACCEDIlst(parms, 1, md);
	// R<-parms[[2]]
	// da parametro
	// N<-parms[[3]]
	CtrlLlst(parms, 2);
	p.N = ACCEDIlst(parms, 2, mi);
	// k<-parms[[4]]
	CtrlLlst(parms, 3);
	p.k = ACCEDIlst(parms, 3, vd);
	// alpha<-parms[[5]]
	CtrlLlst(parms, 4);
	p.alpha = ACCEDIlst(parms, 4, vd);
	// theta<-parms[[6]]
	CtrlLlst(parms, 5);
	p.theta = ACCEDIlst(parms, 5, vd);
	// act.fun<-parms[[7]]
	CtrlLlst(parms, 6);
	p.af = ACCEDIlst(parms, 6, str);
	// EXT.IN<-parms[[8]]
	CtrlLlst(parms, 7);
	p.ei = ACCEDIlst(parms, 7, md);
	// EXT.FUN<-parms[[9]]
	p.ef = ext_fun;
	p.r = R;
	p.dim = LENGTHv_d(X0);
	CtrlLlst(parms, 8);
	p.xmin = ACCEDIlst(parms, 8, vd);
	p.hp = hParsers;
	p.sd = sd_noise;

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

	gsl_odeiv_step *s = gsl_odeiv_step_alloc(T, p.dim);
	gsl_odeiv_control *c = gsl_odeiv_control_y_new(atol, rtol);
	gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(p.dim);

	gsl_odeiv_system sys = { funzione, NULL, LENGTHv_d(X0), &p };

	CREAv_d(g_y, p.dim);
	CREAv_d(g_y_prec, p.dim);
	g_y = copia_v_d(g_y, X0, 1, LENGTHv_d(X0));
	InitVett_d(g_y_prec, 0.0);

	CREAm_d(g_tmp_ris, p.dim, LENGTHv_d(times));

	if (stat_width == 0.0)
		max_stat = LENGTHv_d(times);
	else
		max_stat = (double) stat_width * LENGTHv_d(times);
	for (i = 1; (!stat_width || st < max_stat) && i <= LENGTHv_d(times); i++) {
		ti = ACCEDIv_d(times, i);
		while (t < ti) {
			gsl_odeiv_evolve_apply (e, c, s,
												&sys,
												&t, ti, &h,
												g_y->dati);
		}
#ifdef DET
		fprintf(fp_det, "lsoda step (t = %.16g): ", t);
		for (j = 1; j <= LENGTHv_d(X0); j++) {
			fprintf(fp_det, " %.16g", ACCEDIv_d(g_y, j));
		}
		fprintf(fp_det, "\n");

#endif
		st_all = 1;
		for (j = 1; j <= LENGTHv_d(X0); j++) {
			ASSEGNAm_d(g_tmp_ris, j, i, ACCEDIv_d(g_y, j));
			st_all &= (fabs(ACCEDIv_d(g_y, j) - ACCEDIv_d(g_y_prec, j)) < stat_thr);
		}
		if (st_all)
			st++;
		g_y_prec = copia_v_d(g_y_prec, g_y, 1, LENGTHv_d(g_y));
	}
	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	CREAm_d(ris, LENGTHv_d(X0), i - 1);
	for (j = 1; j < i; j++)
		copia1_m_colonna_d(ris, j, g_tmp_ris, j);
	if (i < LENGTHv_d(times)) {
		*stat = 1;
		Rprintf("Si ferma a %d su %d!\n", i, LENGTHv_d(times));
	}
	else
		*stat = 0;

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "lsoda output:\n");
	fprintf(fp_det, "\tt(out) = ");
	_StampaRawMatr_d(ris);
#endif

	return ris;
}

SEXP lsoda(SEXP parms, SEXP regole, SEXP ext_fun, SEXP X0, SEXP times, SEXP metodo, SEXP atol, SEXP rtol)
{
	int i, ll, ll1, ll2, ll3, nProtected = 0;
	VETTOREd *X01 = NULL, *times1 = NULL;
	MATRICEd *ris1 = NULL;
	LISTA *parms1 = NULL, *regole1 = NULL, *ext_fun1 = NULL;
	GString *metodo1 = NULL;
	double atol1, rtol1;
#ifdef MDEBUG

	GString **nomi;
	char tmp[10];
#endif

	enum TIPO tipi[7];
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** lsoda ***\n");

	ll1 = length(parms);
	ll2 = length(regole);
	ll3 = length(ext_fun);
	ll = max_s_i(ll1, max_s_i(ll2, ll3));
#ifdef MDEBUG
	nomi = (GString **) mia_alloc(ll, GString*);
	if (nomi == NULL) {
		Rprintf("Not enough memory (lsoda # %d, nomi)", __LINE__ - 2);
		error("");
	}
#endif
	X01 = inVETTORE_d(X0, &nProtected);
	times1 = inVETTORE_d(times, &nProtected);
	metodo1 = inSTRINGA(times, &nProtected, "metodo");
	atol1 = NUMERIC_VALUE(atol);
	rtol1 = NUMERIC_VALUE(rtol);
	tipi[0] = MATRd;
	tipi[1] = MATRi;
	tipi[2] = VETTd;
	tipi[3] = VETTd;
	tipi[4] = VETTd;
	tipi[5] = STRINGA;
	tipi[6] = MATRd;
	for (i = 0; i < ll; i++) {
#ifdef MDEBUG
		if (i < ll1)
			snprintf(tmp, 14, "Lista parms %d", i + 1);
		else
			tmp[0] = '\0';
		CREAstr(nomi[i], tmp);
#endif

	}
	parms1 = inLISTA(parms, &nProtected, ll1, tipi, nomi);
	for (i = 0; i < ll2; i++) {
#ifdef MDEBUG

		g_string_printf(nomi[i], "Lista regole %d", i + 1);
#endif

	}
	regole1 = inLISTA(regole, &nProtected, ll2, NULL, nomi);
	if (g_tipi == NULL) {
		g_ntipi = ll3;
		g_tipi = mia_alloc(ll3, enum TIPO);
		if (ll3 > 0 && g_tipi == NULL) {
			Rprintf("Not enough memory (lsoda # %d, tipi)", __LINE__ - 2);
		error("");
		}
	}
	assert(g_ntipi == ll3);
	for (i = 0; i < ll3; i++) {
		g_tipi[i] = STRINGA;
#ifdef MDEBUG
		g_string_printf(nomi[i], "Lista ext_fun %d", i + 1);
#endif
	}
	ext_fun1 = inLISTA(ext_fun, &nProtected, ll3, g_tipi, nomi);
#ifdef MDEBUG

	for (i = 0; i < ll; i++)
		CANCELLAstr(nomi[i]);
	libera(nomi);
#endif

	ris1 = lsoda1(ris1, parms1, regole1, ext_fun1, X01, times1, metodo1->str, atol1, rtol1, NULL, 0, 0, NULL, 0);
	ris = daMATRICE_d(ris1, &nProtected);

	CancellaLISTA(parms1, true);
	CancellaLISTA(regole1, true);
	CancellaLISTA(ext_fun1, true);

	CANCELLAv_d(X01);
	CANCELLAv_d(times1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);
	return ris;
}

