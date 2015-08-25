#include "simulatenet.h"

static const char *nomi_param[MAXPAR] = { "lambda", "alpha", "beta", "Xmin", "Xmax", "X0" };

#define g_M globali.simulatenet.M
#define g_D globali.simulatenet.D
#define g_tmpm1_d globali.simulatenet.tmpm1_d
#define g_Mneg globali.simulatenet.Mneg
#define g_Mdiscr globali.simulatenet.Mdiscr
#define g_reg globali.simulatenet.reg
#define g_aus globali.simulatenet.aus
#define g_genenet globali.simulatenet.genenet
#define g_nulla globali.simulatenet.nulla
#define g_R globali.simulatenet.R
#define g_ris globali.simulatenet.ris

LISTA2 *simulatenet1(LISTA2 *lst, int N, GString *connectivity, int max_reg, double gamma, GString *INdegree, double Cf_cl, VETTOREi *num_subnet, double kappa, GString *f_pr_and, VETTOREd *Xmin, VETTOREd *Xmax, VETTOREd *lambda, VETTOREd *X0, VETTOREd *weight_par, GString *act_fun, VETTOREd *alpha, VETTOREd *beta, VETTOREd *times, GString *method, int save, int ind_itera, int pad_reti, VETTOREd *param_xmin, VETTOREd *param_xzero, muParserHandle_t hParser, double stat_thr, double stat_width, int *stat)
{
	int i, max_lengthR;
	VETTOREd *param[MAXPAR];
	char etich[50], err[256];
	enum TIPO tipi[8];

	_Intestazione("\n*** simulatenet1 ***\n");

	// if (connectivity=="random") g_aus<-connectivityrandom(N=N,max.con=max.g_reg,k=kappa,weight.mean=weight.mean, weight.sd=weight.sd)
	if (!strncmp(connectivity->str, "random", 6))
		g_aus = connectivity_random1(g_aus, N, max_reg, kappa, ACCEDIv_d(weight_par, 1), ACCEDIv_d(weight_par, 2));
	// if (connectivity=="scale free") g_aus<-connectivityscalefree(N=N,gamma=gamma,max.con=max.g_reg,weight.mean=weight.mean, weight.sd=weight.sd,r.tol=0.1,a.tol=1)
	else if (!strncmp(connectivity->str, "scale free", 10))
		g_aus = connectivity_scalefree1(g_aus, N, max_reg, gamma, 0.1, 1.0, ACCEDIv_d(weight_par, 1), ACCEDIv_d(weight_par, 2));
	// if (connectivity=="MTM") g_aus<-connectivitymodular(N=N,gamma=gamma,INdegree=INdegree, Cf.cl=Cf.cl, max.con=max.g_reg,num.subnet=num.subnet, weight.mean=weight.mean, weight.sd=weight.sd,r.tol=0.1,a.tol=1)
	else if (!strncmp(connectivity->str, "MTM", 3))
		g_aus = connectivity_modular1(g_aus, N, max_reg, gamma, INdegree, Cf_cl, num_subnet, 0.1, 1.0, ACCEDIv_d(weight_par, 1), ACCEDIv_d(weight_par, 2));
	// if (connectivity=="geometric") g_aus<-connectivitygeometric(N=N,k=kappa,weight.mean=weight.mean, weight.sd=weight.sd)
	else if (!strncmp(connectivity->str, "geometric", 9))
		g_aus = connectivity_geometric1(g_aus, N, kappa, ACCEDIv_d(weight_par, 1), ACCEDIv_d(weight_par, 2));

	else {
		snprintf(err, 128, "the type of connectivity '%s' is not allowed (the value must be 'random', 'geometric', scale free' or 'MTM')", connectivity->str);
		error(err);
	}
	/*
	g_M<-g_aus[[1]]
	g_Mdiscr<-g_aus[[2]]
	ausR<-createRules(g_M,f.pr.and)
	 	R<-ausR[[1]]
	max.lengthR<-ausR[[2]]
	*/
	CtrlLlst(g_aus, 1);
	g_M = ACCEDIlst(g_aus, 1, md);
	CtrlLlst(g_aus, 2);
	g_Mdiscr = ACCEDIlst(g_aus, 2, mi);
	g_R = createRules1(g_R, g_M, f_pr_and, &max_lengthR, hParser);

	/*
	REG<-matrix(0,ncol=max.lengthR,nrow=N)
	max.L<-0
	for (i in (1:N))
	 {g_aus<-R[[i]]
	  L<-length(g_aus)
	  if (L>0) REG[i,1:L]<-g_aus }
	*/
	CREAm_i(g_reg, N, max_lengthR);
	InitMatr_i(g_reg, 0);
	for (i = 1; i <= N; i++) {
		CtrlLlst(g_R, i);
		if (ACCEDIlst(g_R, i, vi) != NULL) {
			assegna1_mv_riga_i(g_reg, i, ACCEDIlst(g_R, i, vi));
		}
	}

	if (save == 1) {
		snprintf(etich, 50, "Rules%0*d.txt", pad_reti, ind_itera);
		write_m_i(etich, g_reg);
	}

	/*
	g_Mneg<-createNEG(g_Mdiscr)
	 	g_genenet<-list(g_M,R,g_Mneg,lambda,alpha,beta,act.fun,NA,list())
	*/
	g_Mneg = createNEG1(g_Mneg, g_Mdiscr);

	g_tmpm1_d = moltiplica_mm_di(g_tmpm1_d, g_M, g_Mneg);
	if (save == 1) {
		snprintf(etich, 50, "weights%0*d.txt", pad_reti, ind_itera);
		write_m_d(etich, g_tmpm1_d);
	}

	CreaLISTA(g_nulla, tipi, 0);
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
	ASSEGNAlst(g_genenet, 7, vd, NULL);
	CtrlLlst(g_genenet, 8);
	ASSEGNAlst(g_genenet, 8, vd, Xmin);

	if (!strncmp(method->str, "Euler", 5))
		// D<-dinamica(g_genenet,x0,times)*Xmax
		g_D = dinamica1(g_D, g_genenet, g_R, g_nulla, X0, times, NULL, stat_thr, stat_width, 0.0);
	else
		/* out<-lsoda(y=x0, times=times, func=dinamica.lsoda, parms=g_genenet, rtol=1e-3, atol=1e-4)
		out.g_aus<-out[,2:(N+1)]
		D<-t(out.g_aus)*Xmax
		*/
		g_D = lsoda1(g_D, g_genenet, g_R, g_nulla, X0, times, method->str, 1e-4, 1e-3, NULL, stat_thr, stat_width, stat, 0.0);
	moltiplica1_mv_d(g_D, Xmax);

	if (save == 1) {
		snprintf(etich, 50, "SIMdata%0*d.txt", pad_reti, ind_itera);
		write_m_d(etich, g_D);
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
			snprintf(etich, 50, "parameters_%s%0*d.txt", nomi_param[i], pad_reti, ind_itera);
			write_vn_d(etich, param[i], nomi_param[i]);
		}
	}
	// return(list(expr.data=D,weight.matrix=g_M*g_Mneg,Rules=R,parameters=cbind(lambda,alpha,beta)))

	tipi[0] = MATRd;
	tipi[1] = MATRd;
	tipi[2] = VETTd;
	tipi[3] = VETTd;
	tipi[4] = VETTd;
	CreaLISTA(g_ris, tipi, 5);
	CtrlSlst(g_ris, 1);
	ASSEGNAlst(g_ris, 1, md, g_D);
	CtrlSlst(g_ris, 2);
	ASSEGNAlst(g_ris, 2, md, g_tmpm1_d);
	CtrlSlst(g_ris, 3);
	ASSEGNAlst(g_ris, 3, vd, lambda);
	CtrlSlst(g_ris, 4);
	ASSEGNAlst(g_ris, 4, vd, alpha);
	CtrlSlst(g_ris, 5);
	ASSEGNAlst(g_ris, 5, vd, beta);
	//~ CANCELLAm_d(g_M);
	//~ CANCELLAm_i(g_Mneg);
	//~ CANCELLAm_i(g_Mdiscr);
	//~ CANCELLAm_i(g_reg);
	CancellaLISTA(g_aus, false);
	CancellaLISTA(g_genenet, false);
	CancellaLISTA(g_nulla, true);
	if (lst == NULL) {
		lst = mia_alloc(1, LISTA2);
		if (lst == NULL) {
			Rprintf("Not enough memory (simulatenet1 # %d, lst)", __LINE__ - 2);
			error("");
		}
	}
	lst->ris = g_ris;
	lst->r = g_R;

	StrBilanciam();

	return lst;
}

SEXP simulatenet(SEXP N, SEXP connectivity, SEXP max_reg, SEXP gamma, SEXP INdegree, SEXP Cf_cl, SEXP num_subnet, SEXP kappa, SEXP f_pr_and, SEXP act_fun, SEXP alpha, SEXP beta, SEXP lambda, SEXP Xmin, SEXP Xmax, SEXP X0, SEXP weight_par, SEXP param, SEXP times, SEXP stat_thr, SEXP stat_width, SEXP method, SEXP num_reti, SEXP save)
{
	int i, j, n, m, nProtected = 0, save1;
	double min_xzero, max_xzero;
	double stat_thr1, stat_width1;
	int N1, max_reg1, num_reti1, pad_reti;
	double kappa1, Cf_cl1, gamma1, tmp;
	char err[256];
	int stat;
	double p1[MAXPAR], p2[MAXPAR];
	VETTOREd *parm[MAXPAR], *param_orig[2];
	VETTOREi *num_subnet1 = NULL, *param1 = NULL;
	VETTOREd *Xmin1 = NULL, *Xmax1 = NULL, *lambda1 = NULL, *X01 = NULL, *weight_par1 = NULL, *alpha1 = NULL, *beta1 = NULL, *times1 = NULL;
	LISTA2 *lst = NULL;
	muParserHandle_t hParser;
	GString *connectivity1 = NULL, *INdegree1 = NULL, *f_pr_and1 = NULL, *act_fun1 = NULL, *method1 = NULL;
	SEXP R, ris, ret_lista;
	VETTOREd *parm0[MAXPAR] = { NULL, NULL, NULL, NULL, NULL, NULL }; // Xmin e X0 vengono salvati normalizzati

	_InitDbg(false, false, false);

	_Intestazione("\n*** simulatenet ***\n");

	N1 = INTEGER_VALUE(N);
	connectivity1 = inSTRINGA(connectivity, &nProtected, "connectivity");
	max_reg1 = INTEGER_VALUE(max_reg);
	gamma1 = NUMERIC_VALUE(gamma);
	INdegree1 = inSTRINGA(INdegree, &nProtected, "INdegree");
	Cf_cl1 = NUMERIC_VALUE(Cf_cl);
	num_subnet1 = inVETTORE_i(num_subnet, &nProtected);
	kappa1 = NUMERIC_VALUE(kappa);
	f_pr_and1 = inSTRINGA(f_pr_and, &nProtected, "formula");
	Xmin1 = inVETTORE_d(Xmin, &nProtected);
	Xmax1 = inVETTORE_d(Xmax, &nProtected);
	lambda1 = inVETTORE_d(lambda, &nProtected);
	X01 = inVETTORE_d(X0, &nProtected);
	weight_par1 = inVETTORE_d(weight_par, &nProtected);
	if (ACCEDIv_d(weight_par1, 1) < 0 || ACCEDIv_d(weight_par1, 2) < 0)
		error("'weight.par' can not contain values less than zero\n");
	act_fun1 = inSTRINGA(act_fun, &nProtected, "act_fun");
	alpha1 = inVETTORE_d(alpha, &nProtected);
	beta1 = inVETTORE_d(beta, &nProtected);
	times1 = inVETTORE_d(times, &nProtected);
	method1 = inSTRINGA(method, &nProtected, "method");
	param1 = inVETTORE_i(param, &nProtected);
	if (LENGTHv_i(param1) != MAXPAR)
		error("'params' must be a vector of exactly 6 integers\n");
	save1 = INTEGER_VALUE(save);
	num_reti1 = INTEGER_VALUE(num_reti);
	stat_thr1 = NUMERIC_VALUE(stat_thr);
	stat_width1 = NUMERIC_VALUE(stat_width);

	parm[LAMBDA] = lambda1;
	parm0[LAMBDA] = copia_v_d(parm0[LAMBDA], lambda1, 1, LENGTHv_d(lambda1)); // non N1 perche' devo ancora controllare se la lunghezza va bene
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
		if (ACCEDIv_i(param1, j + 1) == ESTERNO) {
			if (LENGTHv_d(parm0[j]) != N1) {
				snprintf(err, 128, "the vector '%s' has %d elements: %d expected", nomi_param[j], LENGTHv_d(parm0[j]), N1);
				error(err);
			}
			if (j == XMAX) {
				for (i = 1; i <= N1; i++) {
					if (ACCEDIv_d(parm0[XMAX], i) <= 0) {
						snprintf(err, 128, "the element %d of vector 'Xmax' (%.16g) is not greater than zero", i, ACCEDIv_d(parm0[XMAX], i));
						warning(err);
					}
					if (ACCEDIv_d(parm0[XMAX], i) <= ACCEDIv_d(parm0[XMIN], i)) {
						snprintf(err, 128, "the element %d of vector 'Xmin' (%.16g) is not lower than the corresponding 'Xmax' (%.16g)", i, ACCEDIv_d(parm[XMIN], i), ACCEDIv_d(parm[XMAX], i));
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

	InitGlobali();
	hParser = InitCalc();

	pad_reti = (int) ceil(log10(num_reti1 + 1));
	Rprintf("Simulatenet: ");
	GetRNGstate();

	// inizializzo i parser
	hParser = InitCalc();

	for (i = 0; i < num_reti1; i++) {
		for (j = 0; j < MAXPAR; j++) {
			switch (ACCEDIv_i(param1, j + 1)) {
				case ESTERNO:
					break;
				case UNIF:
					// nel caso di XZERO ho gia' XMIN e XMAX
					if (j == XZERO) {
						for (n = 1; n <= N1; n++) {
							min_xzero = max_s_d(p1[XZERO], ACCEDIv_d(parm0[XMIN], n));
							max_xzero = min_s_d(p2[XZERO], ACCEDIv_d(parm0[XMAX], n));
							ASSEGNAv_d(parm0[XZERO], n, runif(min_xzero, max_xzero));
						}
					}
					// nel caso di XMAX potrei dover invertire max e min
					else if (j == XMAX) {
						for (n = 1; n <= N1; n++) {
							tmp = runif(p1[XMAX], p2[XMAX]);
							if (tmp >= ACCEDIv_d(parm0[XMIN], n))
								ASSEGNAv_d(parm0[XMAX], n, tmp);
							else {
								ASSEGNAv_d(parm0[XMAX], n, ACCEDIv_d(parm0[XMIN], n));
								ASSEGNAv_d(parm0[XMIN], n, tmp);
							}
						}
					}
					else {
						for (n = 1; n <= N1; n++)
							ASSEGNAv_d((parm0[j]), n, runif(p1[j], p2[j]));
					}
					break;
				case NORM:
					if (j == XZERO)
						error("parameter 'X0' can be only external or uniform, not normal\n");
					else if (j == XMAX) {
						for (n = 1; n <= N1; n++) {
							tmp = rnorm(p1[XMAX], p2[XMAX]);
							if (tmp >= ACCEDIv_d(parm0[XMIN], n))
								ASSEGNAv_d(parm0[XMAX], n, tmp);
							else {
								ASSEGNAv_d(parm0[XMAX], n, ACCEDIv_d(parm0[XMIN], n));
								ASSEGNAv_d(parm0[XMIN], n, tmp);
							}
						}
					}
					else {
						for (n = 1; n <= N1; n++)
							ASSEGNAv_d(parm0[j], n, fabs(rnorm(p1[j], p2[j])));
					}
					break;
				case LOG_NORM:
					if (j == XZERO)
						error("parameter 'X0' can be only external or uniform, not log-normal\n");
					else if (j == XMAX) {
						for (n = 1; n <= N1; n++) {
							tmp = rlnorm(p1[XMAX], p2[XMAX]);
							if (tmp >= ACCEDIv_d(parm0[XMIN], n))
								ASSEGNAv_d(parm0[XMAX], n, tmp);
							else {
								ASSEGNAv_d(parm0[XMAX], n, ACCEDIv_d(parm0[XMIN], n));
								ASSEGNAv_d(parm0[XMIN], n, tmp);
							}
						}
					}
					else {
						for (n = 1; n <= N1; n++)
							ASSEGNAv_d(parm0[j], n, fabs(rlnorm(p1[j], p2[j])));
					}
					break;
				default:
					error("Type of parameter in 'params' not allowed (just 0 = external, 1 = uniform, 2 = normal, 3 = log-normal)\n");
			}
		}
		// normalizzo solo alla fine
		param_orig[0] = copia_v_d(param_orig[0], parm0[XMIN], 1, N1);
		dividi1_vv_d(parm0[XMIN], parm0[XMAX]);

		param_orig[1] = copia_v_d(param_orig[1], parm0[XZERO], 1, N1);
		dividi1_vv_d(parm0[XZERO], parm0[XMAX]);

		for (m = 1; m <= N1; m++) {
			if (ACCEDIv_d(parm0[XZERO], m) > 1 || ACCEDIv_d(parm0[XZERO], m) < ACCEDIv_d(parm0[XMIN], m)) {
				snprintf(err, 256, "the normalized element %d of 'X0' (%.16g) is not between the corresponding normalized 'Xmin' (%.16g) and 1", m, ACCEDIv_d(parm0[XZERO], m), ACCEDIv_d(parm0[XMIN], m));
				warning(err);
			}
		}
		lst = simulatenet1(lst, N1, connectivity1, max_reg1, gamma1, INdegree1, Cf_cl1, num_subnet1, kappa1, f_pr_and1, parm0[XMIN], parm0[XMAX], parm0[LAMBDA], parm0[XZERO], weight_par1, act_fun1, parm0[ALPHA], parm0[BETA], times1, method1, save1, i + 1, pad_reti, param_orig[0], param_orig[1], hParser, stat_thr1, stat_width1, &stat);
		if (stat)
			Rprintf("-");
		else
			Rprintf(".");
		R_FlushConsole();
	}
	Rprintf(" done\n");
	PutRNGstate();

	R = daLISTA(lst->r, &nProtected);
	ris = daLISTA(lst->ris, &nProtected);
	libera(lst);

	// costruisco a mano la lista
	PROTECT(ret_lista = allocVector(VECSXP, 2));
	++nProtected;
	SET_VECTOR_ELT(ret_lista, 0, ris);
	SET_VECTOR_ELT(ret_lista, 1, R);

	CancGlobali();

	CANCELLAv_i(num_subnet1);
	CANCELLAv_d(Xmin1);
	CANCELLAv_d(Xmax1);
	CANCELLAv_d(X01);
	CANCELLAv_d(alpha1);
	CANCELLAv_d(beta1);
	CANCELLAv_d(lambda1);
	CANCELLAv_d(times1);
	// lambda, alpha e beta sono gia` stati cancellati con la lista "ris"
	CANCELLAv_d(parm0[XMIN]);
	CANCELLAv_d(parm0[XMAX]);
	CANCELLAv_d(weight_par1);
	CANCELLAv_i(param1);
	CANCELLAv_d(param_orig[0]);
	CANCELLAv_d(param_orig[1]);
	CANCELLAstr(connectivity1);
	CANCELLAstr(INdegree1);
	CANCELLAstr(f_pr_and1);
	CANCELLAstr(act_fun1);
	CANCELLAstr(method1);

	DeInitCalc(hParser);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ret_lista;
}
