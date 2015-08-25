#include "simulateprofiles.h"

#define g_tipi globali.tipi
#define g_ntipi globali.n_tipi

#define g_ind globali.simulateprofiles.ind
#define g_tmpm1_i globali.simulateprofiles.tmpm1_i
#define g_tmpm1_d globali.simulateprofiles.tmpm1_d
#define g_reg globali.simulateprofiles.reg
#define g_M globali.simulateprofiles.M
#define g_Mdiscr globali.simulateprofiles.Mdiscr
#define g_Mneg globali.simulateprofiles.Mneg
#define g_R1 globali.simulateprofiles.R1
#define g_genenet globali.simulateprofiles.genenet

static const char *nomi_param[MAXPAR] = { "lambda", "alpha", "beta", "Xmin", "Xmax", "X0" };

int Esiste(const char *strPath, const char *nome)
{
	char tmp[256];

	strncpy(tmp, strPath, 256);
	strncat(tmp, nome, 256);
	if (g_file_test(nome, G_FILE_TEST_EXISTS)) {
		if(g_file_test(nome, G_FILE_TEST_IS_DIR))
			return 1;
		else
			return 2;
	}
	return 0;
}

int Conta(const char *dir, const char *pref)
{
	struct dirent *dp;
	int tot = 0;

	DIR *dfd = opendir(dir);
	if(dfd != NULL) {
		while((dp = readdir(dfd)) != NULL) {
			if (!strncmp(pref, dp->d_name, strlen(pref))) {
				tot++;
			}
		}
		closedir(dfd);
	}
	return tot;
}

MATRICEd *simulateprofiles1(MATRICEd *D, MATRICEd *weights, LISTA *R, GString *f_pr_and, VETTOREd *X0, VETTOREd *Xmin, VETTOREd *Xmax, VETTOREd *lambda, GString *act_fun, VETTOREd *alpha, VETTOREd *beta, double sd_noise, VETTOREd *times, GString *method, MATRICEd *ext_in, LISTA *ext_fun, int save, int rete, int ind_itera, int ko, int pad_reti, int pad_exp, int pad_ko, VETTOREd *param_xmin, VETTOREd *param_xzero, muParserHandle_t *hParsers, double stat_thr, double stat_width, int *stat)
{
	int i;
	VETTOREd *param[MAXPAR];
	enum TIPO tipi[9];
	char etich[50];

	_Intestazione("\n*** simulateprofiles1 ***\n");

	/*
	if(is.null(Rules)) {
		g_Mdiscr<-weights
	        g_ind<-which(g_Mdiscr!=0,arr.g_ind=TRUE)
	        g_Mdiscr[g_ind]<-1
	        ausR<-createRules(g_Mdiscr,f.pr.and)
		Rules<-ausR[[1]]
		max.lengthR<-ausR[[2]]
	        REG<-matrix(0,ncol=max.lengthR,nrow=N)
		max.L<-0
		for (i in (1:N)) {
			aus<-Rules[[i]]
			L<-length(aus)
			if (L>0) REG[i,1:L]<-aus
	   } ...
	*/

	tipi[0] = MATRd;
	tipi[1] = MATRi;
	tipi[2] = VETTd;
	tipi[3] = VETTd;
	tipi[4] = VETTd;
	tipi[5] = STRINGA;
	tipi[6] = MATRd;
	tipi[7] = VETTd;
	if (g_genenet == NULL) {
		CreaLISTA(g_genenet, tipi, 8);
		CtrlLlst(g_genenet, 6);
		CREAstr(ACCEDIlst(g_genenet, 6, str), "");
	}

	/*
	g_M<-weights
	g_Mdiscr<-M
	g_ind<-which(g_M!=0,arr.g_ind=TRUE)
	g_Mdiscr[g_ind]<-1
	g_Mneg<-g_Mdiscr*sign(g_M)
	g_M<-abs(g_M)
	R<-Rules
	*/

	g_ind = which_m_indxne_d(g_ind, weights, 0.0);
	g_Mdiscr = arrotonda_m_d(g_Mdiscr, weights);
	assegna1_ms_indx_i(g_Mdiscr, g_ind, 1);
	g_tmpm1_i = segno_m_d(g_tmpm1_i, weights);
	g_Mneg = moltiplica_mm_i(g_Mneg, g_Mdiscr, g_tmpm1_i);
	g_M = abs_m_d(g_M, weights);

	// genenet<-list(g_M,R,g_Mneg,lambda,alpha,beta,act.fun,EXT.IN,EXT.FUN)
	CtrlSlst(g_genenet, 1);
	ASSEGNAlst(g_genenet, 1, md, g_M);
	CtrlSlst(g_genenet, 2);
	ASSEGNAlst(g_genenet, 2, mi, g_Mneg);
	CtrlSlst(g_genenet, 3);
	ASSEGNAlst(g_genenet, 3, vd, lambda);
	CtrlSlst(g_genenet, 4);
	ASSEGNAlst(g_genenet, 4, vd, alpha);
	CtrlSlst(g_genenet, 5);
	ASSEGNAlst(g_genenet, 5, vd, beta);
	CtrlLlst(g_genenet, 6);
	g_string_assign(ACCEDIlst(g_genenet, 6, str), act_fun->str);
	CtrlSlst(g_genenet, 7);
	ASSEGNAlst(g_genenet, 7, md, ext_in);
	CtrlLlst(g_genenet, 8);
	ASSEGNAlst(g_genenet, 8, vd, Xmin);

	if (!strncmp(method->str, "Euler", 5))
		// D<-dinamica(g_genenet,x0,times)*Xmax
		D = dinamica1(D, g_genenet, R, ext_fun, X0, times, hParsers, stat_thr, stat_width, sd_noise);
	else
		/*
		out<-lsoda(y=x0, times=times, func=dinamica.lsoda, parms=g_genenet, rtol=1e-3, atol=1e-4)
		out.aus<-out[,2:(N+1)]
		D<-t(out.aus)*Xmax
		*/
		D = lsoda1(D, g_genenet, R, ext_fun, X0, times, method->str, 1e-4, 1e-3, hParsers, stat_thr, stat_width, stat, sd_noise);

	moltiplica1_mv_d(D, Xmax);

	if (save == 1) {
		if (ko >= 0)
			snprintf(etich, 50, "SIMdata_simulateprofiles%0*d_%0*d_ko%0*d.txt", pad_reti, rete, pad_exp, ind_itera, pad_ko, ko);
		else
			snprintf(etich, 50, "SIMdata_simulateprofiles%0*d_%0*d.txt", pad_reti, rete, pad_exp, ind_itera);
		write_m_d(etich, D);
	}

	param[LAMBDA] = lambda;
	param[ALPHA] = alpha;
	param[BETA] = beta;
	// fa differenza se non sono esterni, perche' in quel caso sono lunghi 2; salvo Xmin e X0 non normalizzati
	param[XMIN] = param_xmin;
	param[XMAX] = Xmax;
	param[XZERO] = param_xzero;

	if (save == 1) {
		for (i = 0; i < MAXPAR; i++) {
			if (ko >= 0)
				snprintf(etich, 50, "parameters_%s_simulateprofiles%0*d_%0*d_ko%0*d.txt", nomi_param[i], pad_reti, rete, pad_exp, ind_itera, pad_ko, ko);
			else
				snprintf(etich, 50, "parameters_%s_simulateprofiles%0*d_%0*d.txt", nomi_param[i], pad_reti, rete, pad_exp, ind_itera);
			write_vn_d(etich, param[i], nomi_param[i]);
		}
	}

	//~ if (LENGTHlst(R) == 0)
	//~ CancellaLISTA(R1, true);
	//~ CancellaLISTA(genenet, false);
	//~ CANCELLAv_i(g_ind);
	//~ CANCELLAm_d(g_tmpm1_d);
	//~ CANCELLAm_i(g_reg);
	//~ CANCELLAm_d(g_M);
	//~ CANCELLAm_d(g_Mdiscr);
	//~ CANCELLAm_d(g_Mneg);

	StrBilanciam();

	return D;
}

SEXP simulateprofiles(SEXP N, SEXP weights, SEXP Rules, SEXP f_pr_and, SEXP act_fun, SEXP alpha, SEXP beta, SEXP lambda, SEXP Xmin, SEXP Xmax, SEXP X0, SEXP param, SEXP ko_experim, SEXP sd_noise, SEXP times, SEXP stat_thr, SEXP stat_width, SEXP method, SEXP ext_in, SEXP ext_fun, SEXP num_exp, SEXP save)
{
	int i, j, k, k1, m, f, nf, ll, llr, llf, nProtected = 0, N1;
	double min_xzero, max_xzero;
	double stat_thr1, stat_width1;
	int stat;
	int save1, num_exp1, ko = 0, max_lengthR;
	double p1[MAXPAR], p2[MAXPAR], tmp;
	VETTOREd *parm[MAXPAR], *param_orig[2];
	double sd_noise1;
	VETTOREd *lambda1 = NULL, *X01 = NULL, *X02 = NULL, *Xmin1 = NULL, *Xmax1 = NULL, *alpha1 = NULL, *beta1 = NULL, *times1 = NULL;
	VETTOREi *param1 = NULL, *ko_experim1 = NULL, *ko_experim2 = NULL;
	MATRICEd *weights1 = NULL, *weights2 = NULL, *ext_in1 = NULL, *ris1 = NULL;
	LISTA *r1 = NULL, *r2 = NULL, *ext_fun1 = NULL;
	GString *f_pr_and1 = NULL, *act_fun1 = NULL, *method1 = NULL;
	char buf[25];
	char etich[50];
	char err[256];
	int pad_reti, pad_exp, pad_ko, num_ko;
	bool regole_nulle = false;
	SEXP ris;
	muParserHandle_t hParser, *hParsers = NULL;
#ifdef MDEBUG
	GString **nomi;
	char tmp1[10];
#endif
	VETTOREd *parm0[MAXPAR] = { NULL, NULL, NULL, NULL, NULL, NULL };

	_InitDbg(false, false, false);

	_Intestazione("\n*** simulateprofiles ***\n");

	N1 = INTEGER_VALUE(N);
	weights1 = inMATRICE_d(weights, &nProtected);
	llr = length(Rules);
	llf = length(ext_fun);
	ll = max_s_i(llr, llf);
#ifdef MDEBUG

	nomi = mia_alloc(ll, GString *);
	if (ll > 0 && nomi == NULL) {
		Rprintf("Not enough memory (simulateprofiles # %d, nomi)", __LINE__ - 2);
		error("");
	}
#endif
	for (i = 0; i < ll; i++) {
#ifdef MDEBUG
		if (i < llr)
			snprintf(tmp1, 10, "Lista r %d", i + 1);
		else
			tmp1[0] = '\0';
		CREAstr(nomi[i], tmp1);
#endif

	}
	r2 = inLISTA(Rules, &nProtected, llr, NULL, nomi);
	if (g_tipi == NULL) {
		g_ntipi = llf;
		g_tipi = mia_alloc(llf, enum TIPO);
		if (llf > 0 && g_tipi == NULL) {
			Rprintf("Not enough memory (simulateprofiles # %d, tipi)", __LINE__ - 2);
			error("");
		}
	}
	assert(g_ntipi == llf);
	for (i = 0; i < llf; i++) {
		g_tipi[i] = STRINGA;
#ifdef MDEBUG
		g_string_printf(nomi[i], "Lista ext_fun %d", i + 1);
#endif

	}
	ext_fun1 = inLISTA(ext_fun, &nProtected, llf, g_tipi, nomi);
#ifdef MDEBUG

	for (i = 0; i < ll; i++)
		CANCELLAstr(nomi[i]);
	libera(nomi);
#endif

	f_pr_and1 = inSTRINGA(f_pr_and, &nProtected, "formula");
	X01 = inVETTORE_d(X0, &nProtected);
	Xmin1 = inVETTORE_d(Xmin, &nProtected);
	Xmax1 = inVETTORE_d(Xmax, &nProtected);
	lambda1 = inVETTORE_d(lambda, &nProtected);
	act_fun1 = inSTRINGA(act_fun, &nProtected, "act_fun");
	alpha1 = inVETTORE_d(alpha, &nProtected);
	beta1 = inVETTORE_d(beta, &nProtected);
	sd_noise1 = NUMERIC_VALUE(sd_noise);
	times1 = inVETTORE_d(times, &nProtected);
	method1 = inSTRINGA(method, &nProtected, "method");
	ext_in1 = inMATRICE_d(ext_in, &nProtected);
	param1 = inVETTORE_i(param, &nProtected);
	if (LENGTHv_i(param1) != MAXPAR)
		error("'params' must be a vector of exactly 6 integers\n");
	ko_experim1 = inVETTORE_i(ko_experim, &nProtected);
	save1 = INTEGER_VALUE(save);
	num_exp1 = INTEGER_VALUE(num_exp);
	stat_thr1 = NUMERIC_VALUE(stat_thr);
	stat_width1 = NUMERIC_VALUE(stat_width);

	parm[LAMBDA] = lambda1;
	parm0[LAMBDA] = copia_v_d(parm0[LAMBDA], lambda1, 1, LENGTHv_d(lambda1)); // non N1 perche´ devo ancora controllare se la lunghezza va bene
	parm[ALPHA] = alpha1;
	parm0[ALPHA] = copia_v_d(parm0[ALPHA], alpha1, 1, LENGTHv_d(alpha1));
	parm[BETA] = beta1;
	parm0[BETA] = copia_v_d(parm0[BETA], beta1, 1, LENGTHv_d(beta1));
	parm[XMIN] = Xmin1;
	parm0[XMIN] = copia_v_d(parm0[XMIN], Xmin1, 1, LENGTHv_d(Xmin1));
	parm[XMAX] = Xmax1;
	parm0[XMAX] = copia_v_d(parm0[XMAX], Xmax1, 1, LENGTHv_d(Xmax1));
	parm[XZERO] = X01;
	parm0[XZERO] = copia_v_d(parm0[XZERO], X01, 1, LENGTHv_d(X01));

	param_orig[0] = NULL;
	param_orig[1] = NULL;

	for (j = 0; j < MAXPAR; j++) {
		p1[j] = ACCEDIv_d((parm[j]), 1);
		p2[j] = ACCEDIv_d((parm[j]), 2);
		if (ACCEDIv_i(param1, j + 1) == ESTERNO
				|| ACCEDIv_i(param1, j + 1) == ESTERNO_UNIF
				|| ACCEDIv_i(param1, j + 1) == ESTERNO_NORM
				|| ACCEDIv_i(param1, j + 1) == ESTERNO_LOG_NORM) {
			if (ACCEDIv_i(param1, j + 1) == ESTERNO && LENGTHv_d(parm0[j]) != N1) {
				snprintf(err, 128, "the vector '%s' has %d elements: %d expected", nomi_param[j], LENGTHv_d(parm0[j]), N1);
				error(err);
			}
			else if (ACCEDIv_i(param1, j + 1) != ESTERNO && LENGTHv_d(parm0[j]) != N1 + 2) {
				snprintf(err, 128, "the vector '%s' has %d elements: %d expected (external parameters plus two additional distribution parameters at the end)", nomi_param[j], LENGTHv_d(parm0[j]), N1 + 2);
				error(err);
			}
			// raccolgo i parametri delle distribuzioni
			if (ACCEDIv_i(param1, j + 1) == ESTERNO_UNIF
					|| ACCEDIv_i(param1, j + 1) == ESTERNO_NORM
					|| ACCEDIv_i(param1, j + 1) == ESTERNO_LOG_NORM) {
				p1[j] = ACCEDIv_d(parm[j], N1 - 2);
				p2[j] = ACCEDIv_d(parm[j], N1 - 1);
			}
			if (j == XMAX) {
				// ignoro gli ultimi due elementi (infatti N1 e` la lunghezza giusta)
				for (i = 1; i <= N1; i++) {
					if (ACCEDIv_d(parm0[XMAX], i) <= 0) {
						snprintf(err, 128, "the element %d of vector 'Xmax' (%.16g) is not greater than zero", i, ACCEDIv_d(parm0[XMAX], i));
						warning(err);
					}
					if (ACCEDIv_d(parm0[XMAX], i) <= ACCEDIv_d(parm0[XMIN], i)) {
						snprintf(err, 128, "the element %d of vector 'Xmin' (%.16g) is not lower than the corresponding 'Xmax' (%.16g)", i, ACCEDIv_d(parm0[XMIN], i), ACCEDIv_d(parm0[XMAX], i));
						error(err);
					}
				}
			}
			else if (j == XZERO) {
				for (i = 1; i <= N1; i++) {
					if ((ACCEDIv_d(parm0[XZERO], i) <= ACCEDIv_d(parm0[XMIN], i)) ||
							(ACCEDIv_d(parm0[XZERO], i) >= ACCEDIv_d(parm0[XMAX], i))) {
						snprintf(err, 128, "the element %d of vector 'X0' (%.16g) is not included in the corresponding ['Xmin', 'Xmax'] range, which is [%.16g, %.16g]", i, ACCEDIv_d(parm[XZERO], i), ACCEDIv_d(parm[XMIN], i), ACCEDIv_d(parm[XMAX], i));
						warning(err);
					}
				}
			}
		}
		else
			CREAv_d(parm0[j], N1);
	}
	if (!isNull(ko_experim)) {
		ko = 1;
		if (ACCEDIv_i(ko_experim1, 1) == 0)
			ko_experim2 = seq_i(ko_experim2, 1, N1, 1);
		else
			ko_experim2 = copia_v_i(ko_experim2, ko_experim1, 1, LENGTHv_i(ko_experim1));
	}
	else
		ko = -1;
	nf = 1;
	if (weights1 == NULL) {
		//~ Rprintf("i pesi sono nulli: considero le ");
		nf = Conta("./", "weights");
		//~ Rprintf("%d reti generate da SN\n", nf);
	}
	if (r2 == NULL) {
		regole_nulle = true;
		//~ Rprintf("le regole sono nulle: le leggero` da file\n");
	}
	InitGlobali();
	GetRNGstate();
	pad_reti = (int) ceil(log10(nf + 1));
	pad_exp = (int) ceil(log10(num_exp1 + 1));
	num_ko = ((ko >= 0) ? LENGTHv_i(ko_experim2) : 1);
	pad_ko = (int) ceil(log10(num_ko + 1));
	// inizializzo i parser
	hParser = InitCalc();
	hParsers = mia_alloc(LENGTHm2_d(ext_in1), muParserHandle_t);
	for (i = 0; i < LENGTHm2_d(ext_in1); i++)
		hParsers[i] = InitCalc();

	for (f = 0; f < nf; f++) {
		Rprintf("Simulateprofiles");
		if (!isNull(weights))
			Rprintf(":\n");
		else {
			Rprintf(" (network %d of %d):\n", f + 1, nf);
			snprintf(buf, 25, "weights%0*d.txt", pad_reti, f + 1);
			weights1 = read_m_d(weights1, N1, buf);
		}
		if (regole_nulle) {
			if (!isNull(weights))
				error("To read rules from file 'weights' must be NULL!\n");
			snprintf(buf, 25, "Rules%0*d.txt", pad_reti, f + 1);
			r1 = read_lst_i(r1, N1, buf);
		}
		else if (LENGTHlst(r2) == 0) {
			//~ Rprintf("\nla lista delle regole ha lunghezza zero: la genero casualmente\n");
			g_Mdiscr = arrotonda_m_d(g_Mdiscr, weights1);
			g_ind = which_m_indxne_i(g_ind, g_Mdiscr, 0);
			// questa copia server perche´ createRules1 vuole double, mentre Mdiscr e` di interi
			g_tmpm1_d = copia_m_d(g_tmpm1_d, weights1);
			assegna1_ms_indx_d(g_tmpm1_d, g_ind, 1.0);

			r1 = createRules1(r1, g_tmpm1_d, f_pr_and1, &max_lengthR, hParser);
			CREAm_i(g_reg, LENGTHm1_d(weights1), max_lengthR);
			InitMatr_i(g_reg, 0);
			for (i = 1; i <= LENGTHm1_d(weights1); i++) {
				CtrlLlst(r1, i);
				if (ACCEDIlst(r1, i, vi) != NULL)
					assegna1_mv_riga_i(g_reg, i, ACCEDIlst(r1, i, vi));
			}
			if (save1 == 1) {
				snprintf(etich, 50, "Rules_simulateprofiles%0*d.txt", pad_reti, f + 1);
				write_m_i(etich, g_reg);
			}
		}
		else {
			r1 = r2;
		}
		for (j = 0; j < MAXPAR; j++) {
			if (ACCEDIv_i(param1, j + 1) == SN
					|| ACCEDIv_i(param1, j + 1) == SN_UNIF
					|| ACCEDIv_i(param1, j + 1) == SN_NORM
					|| ACCEDIv_i(param1, j + 1) == SN_LOG_NORM) {
				if (!isNull(weights))
					error("To read parameters from SN files 'weights' must be NULL!\n");
				snprintf(buf, 25, "parameters_%s%0*d.txt", nomi_param[j], pad_reti, f + 1);
				// riutilizzo i vettori param0, ma con dimensioni maggiori
				CREAv_d(parm0[j], N1);
				parm0[j] = read_vn_d(parm0[j], buf);
				parm0[XMAX]->dim = N1;
			}
		}
		for (i = 0; i < num_exp1; i++) {
			if (ko >= 0)
				Rprintf("%d) ", i + 1);
			for (k = 0; k <= num_ko; k++) {
				for (j = 0; j < MAXPAR; j++) {
					switch (ACCEDIv_i(param1, j + 1)) {
						case ESTERNO:
							break;
						case UNIF:
						case ESTERNO_UNIF:
						case SN_UNIF:
							// nel caso di XZERO ho gia' XMIN e XMAX
							if (j == XZERO) {
								for (m = 1; m <= N1; m++) {
									if (ACCEDIv_i(param1, j + 1) == UNIF) {
										min_xzero = max_s_d(p1[XZERO], ACCEDIv_d(parm0[XMIN], m));
										max_xzero = min_s_d(p2[XZERO], ACCEDIv_d(parm0[XMAX], m));
									}
									else {
										min_xzero = max_s_d(ACCEDIv_d(parm0[XZERO], m) + p1[XZERO], ACCEDIv_d(parm0[XMIN], m));
										max_xzero = min_s_d(ACCEDIv_d(parm0[XZERO], m) + p2[XZERO], ACCEDIv_d(parm0[XMAX], m));
									}
									ASSEGNAv_d(parm0[XZERO], m, runif(min_xzero, max_xzero));
								}
							}
							// nel caso di XMAX potrei dover invertire max e min
							else if (j == XMAX) {
								for (m = 1; m <= N1; m++) {
									if (ACCEDIv_i(param1, j + 1) == UNIF)
										tmp = runif(p1[XMAX], p2[XMAX]);
									else
										tmp = ACCEDIv_d(parm0[XMAX], m) + runif(p1[XMAX], p2[XMAX]);
									if (tmp >= ACCEDIv_d(parm0[XMIN], m))
										ASSEGNAv_d(parm0[XMAX], m, tmp);
									else {
										ASSEGNAv_d(parm0[XMAX], m, ACCEDIv_d(parm0[XMIN], m));
										ASSEGNAv_d(parm0[XMIN], m, tmp);
									}
								}
							}
							else {
								if (ACCEDIv_i(param1, j + 1) == UNIF) {
									for (m = 1; m <= N1; m++)
										ASSEGNAv_d((parm0[j]), m, runif(p1[j], p2[j]));
								}
								else {
									for (m = 1; m <= N1; m++)
										ASSEGNAv_d(parm0[j], m, ACCEDIv_d(parm0[j], m) + runif(p1[j], p2[j]));
								}
							}
							break;
						case NORM:
						case ESTERNO_NORM:
						case SN_NORM:
							if (j == XMAX) {
								for (m = 1; m <= N1; m++) {
									if (ACCEDIv_i(param1, j + 1) == NORM)
										tmp = rnorm(p1[XMAX], p2[XMAX]);
									else
										tmp = ACCEDIv_d(parm0[XMAX], m) + rnorm(p1[XMAX], p2[XMAX]);
									if (tmp >= ACCEDIv_d(parm0[XMIN], m))
										ASSEGNAv_d(parm0[XMAX], m, tmp);
									else {
										ASSEGNAv_d(parm0[XMAX], m, ACCEDIv_d(parm0[XMIN], m));
										ASSEGNAv_d(parm0[XMIN], m, tmp);
									}
								}
							}
							else {
								if (ACCEDIv_i(param1, j + 1) == NORM) {
									for (m = 1; m <= N1; m++)
										ASSEGNAv_d(parm0[j], m, rnorm(p1[j], p2[j]));
								}
								else {
									for (m = 1; m <= N1; m++)
										ASSEGNAv_d(parm0[j], m, ACCEDIv_d(parm0[j], m) + rnorm(p1[j], p2[j]));
								}
							}
							break;
						case LOG_NORM:
						case ESTERNO_LOG_NORM:
						case SN_LOG_NORM:
							if (j == XMAX) {
								for (m = 1; m <= N1; m++) {
									if (ACCEDIv_i(param1, j + 1) == LOG_NORM)
										tmp = rlnorm(p1[XMAX], p2[XMAX]);
									else
										tmp = ACCEDIv_d(parm0[XMAX], m) + rlnorm(p1[XMAX], p2[XMAX]);
									if (tmp >= ACCEDIv_d(parm0[XMIN], m))
										ASSEGNAv_d(parm0[XMAX], m, tmp);
									else {
										ASSEGNAv_d(parm0[XMAX], m, ACCEDIv_d(parm0[XMIN], m));
										ASSEGNAv_d(parm0[XMIN], m, tmp);
									}
								}
							}
							else {
								if (ACCEDIv_i(param1, j + 1) == LOG_NORM) {
									for (m = 1; m <= N1; m++)
										ASSEGNAv_d(parm0[j], m, rlnorm(p1[j], p2[j]));
								}
								else {
									for (m = 1; m <= N1; m++)
										ASSEGNAv_d(parm0[j], m, ACCEDIv_d(parm0[j], m) + rlnorm(p1[j], p2[j]));
								}
							}
							break;
						case SN:
							break;
						default:
							error("Type of parameter in 'params' not allowed:\n\
								\t0 = external, 1 = uniform, 2 = normal, 3 = log-normal, 4 = from simulatenet,\n\
								\t5 = external + uniform, 6 = external + normal, 7 = external + log-normal,\n\
								\t8 = from simulatenet + uniform, 9 = from simulatenet + normal, 10 = from simulatenet + log-normal)\n");
					}
				}

				// normalizzo solo alla fine
				param_orig[0] = copia_v_d(param_orig[0], parm0[XMIN], 1, N1);
				dividi1_vv_d(parm0[XMIN], parm0[XMAX]);

				param_orig[1] = copia_v_d(param_orig[1], parm0[XZERO], 1, N1);
				dividi1_vv_d(parm0[XZERO], parm0[XMAX]);

				for (m = 1; m <= N1; m++) {
					if (ACCEDIv_d(parm0[XZERO], m) > 1 || ACCEDIv_d(parm0[XZERO], m) < ACCEDIv_d(parm0[XMIN], m)) {
						if (ACCEDIv_i(param1, j + 1) > LOG_NORM) {
							if (ACCEDIv_d(parm0[XZERO], m) > 1) {
								snprintf(err, 256, "the normalized element %d of 'X0' (%.16g) is greater than 1 and it has been set to 1", m, ACCEDIv_d(parm0[XZERO], m));
								warning(err);
								ASSEGNAv_d(parm0[XZERO], m, 1);
							}
							else {
snprintf(err, 256, "the normalized element %d of 'X0' (%.16g) is lower than the corresponding normalized 'Xmin' (%.16g) and it has been set equal to it", m, ACCEDIv_d(parm0[XZERO], m), ACCEDIv_d(parm0[XMIN], m));
								warning(err);
								ASSEGNAv_d(parm0[XZERO], m, ACCEDIv_d(parm0[XMIN], m));
							}
						}
						else {
							snprintf(err, 256, "the normalized element %d of 'X0' (%.16g) is not between the corresponding normalized 'Xmin' and 'Xmax' (%.16g) and 1", m, ACCEDIv_d(parm0[XZERO], m), ACCEDIv_d(parm0[XMIN], m));
							warning(err);
						}
					}
				}

				if (ko == -1) {
					ris1 = simulateprofiles1(ris1, weights1, r1, f_pr_and1, parm0[XZERO], parm0[XMIN], parm0[XMAX], parm0[LAMBDA], act_fun1, parm0[ALPHA], parm0[BETA], sd_noise1, times1, method1, ext_in1, ext_fun1, save1, f + 1, i + 1, -1, pad_reti, pad_exp, pad_ko, param_orig[0], param_orig[1], hParsers, stat_thr1, stat_width1, &stat);
					if (stat)
						Rprintf("-");
					else
						Rprintf(".");
					break;
				}
				else {
					// w[,i] <- 0; x0[i] <- 0
					weights2 = copia_m_d(weights2, weights1);
					X02 = copia_v_d(X02, parm0[XZERO], 1, N1);
					if (k > 0) {
						assegna1_ms_riga_d(weights2, k, 0.0);
						ASSEGNAv_d(X02, k, 0.0);
					}
					k1 = (k > 0) ? ACCEDIv_i(ko_experim2, k) : 0;
					ris1 = simulateprofiles1(ris1, weights2, r1, f_pr_and1, X02, parm0[XMIN], parm0[XMAX], parm0[LAMBDA], act_fun1, parm0[ALPHA], parm0[BETA], sd_noise1, times1, method1, ext_in1, ext_fun1, save1, f + 1, i + 1, k1, pad_reti, pad_exp, pad_ko, param_orig[0], param_orig[1], hParsers, stat_thr1, stat_width1, &stat);
					if (stat)
						Rprintf(" %d", k1);
					else
						Rprintf(" %d-", k1);
				}
				R_FlushConsole();
			}
			if (ko >= 0)
				Rprintf("\n");
			R_FlushConsole();
		}
		Rprintf("\n");
	}
	Rprintf("done\n");
	PutRNGstate();
	CancGlobali();
	CancellaLISTA(globali.simulateprofiles.genenet, false);

	ris = daMATRICE_d(ris1, &nProtected);

	for (i = 0; i < LENGTHm2_d(ext_in1); i++) {
		DeInitCalc(hParsers[i]);
	}
	libera(hParsers);
	DeInitCalc(hParser);

	CancellaLISTA(r1, true);
	CANCELLAm_d(weights1);
	CANCELLAm_d(weights2);
	CANCELLAm_d(ext_in1);
	CancellaLISTA(ext_fun1, true);
	CANCELLAv_d(Xmin1);
	CANCELLAv_d(Xmax1);
	CANCELLAv_d(X01);
	CANCELLAv_d(X02);
	CANCELLAv_d(alpha1);
	CANCELLAv_d(beta1);
	CANCELLAv_d(lambda1);
	CANCELLAv_i(param1);
	CANCELLAv_d(times1);
	CANCELLAv_d(param_orig[0]);
	CANCELLAv_d(param_orig[1]);
	for (j = 0 ; j < MAXPAR; j++) {
		CANCELLAv_d(parm0[j]);
	}
	CANCELLAstr(f_pr_and1);
	CANCELLAstr(act_fun1);
	CANCELLAstr(method1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
