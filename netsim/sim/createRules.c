#include "createRules.h"

#define g_op globali.createRules.op
#define g_tmp globali.createRules.tmp

LISTA *createRules1(LISTA *ris, const MATRICEd *m, const GString *f_pr_and, int *max_lengthL, muParserHandle_t hParser)
{
	int i, ll, lengthL = 0;
	GString *errore;

	_InitDbg(false, false, false);

	_Intestazione("\n***createRules***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tM = ");
	_StampaRawMatr_d(m);
	fprintf(fp_det, "\tf.pr.and = ...\n");
#endif

	ll = LENGTHm1_d(m);
	CreaLISTA(ris, NULL, ll);
	// max.lengthL<-0
	*max_lengthL = 0;
	if (f_pr_and->len > 0) {
#ifdef FDEBUG
		fprintf(fp_fdbg, "Compilo la funzione f_pr_and definita come '%s' ... ", f_pr_and->str);
#endif
		mupSetExpr(hParser, f_pr_and->str);
#ifdef FDEBUG
		fprintf(fp_fdbg, "ok\n");
#endif
	}
	// for(i in 1:n){
	for (i = 1; i <= LENGTHm1_d(m); i++) {
		// g_op <- which(M[i, ] != 0)
		// l_op <- length(g_op)
		g_op = which_m_rowindxne_d(g_op, m, i, 0.0);
		// if (l_op>1)
		if (LENGTHv_i(g_op) > 1) {
			g_tmp = create_logicRule(g_tmp, g_op, (f_pr_and->len > 0), hParser);
			lengthL = LENGTHv_i(g_tmp);
			CtrlLlst(ris, i);
			ACCEDIlst(ris, i, vi) = copia_v_i(ACCEDIlst(ris, i, vi), g_tmp, 1, lengthL);
		}
		else {
			if (LENGTHv_i(g_op) == 1) {
				lengthL = 1;
				CtrlLlst(ris, i);
				ACCEDIlst(ris, i, vi) = copia_v_i(ACCEDIlst(ris, i, vi), g_op, 1, 1);
			}
			else {
				// logic_rule <- NaN
				//ASSEGNAv_i(g_tmp, 1, R_NaN);
				lengthL = 0;
				CtrlSlst(ris, i);
				ASSEGNAlst(ris, i, vi, NULL);
			}
		}
		//~ lengthL = length(logic_rule);
		//~ max.lengthL<-max(max.lengthL,lengthL)
		*max_lengthL = max_s_i(*max_lengthL, lengthL);
	// L[[i]] <- logic_rule
	}
	//~ CANCELLAv_i(g_op);
	//~ CANCELLAv_i(g_tmp);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "createRules output:\n");
	fprintf(fp_det, "\tlist(L, max.lengthL) = ");
	_StampaRawLista(ris);
	fprintf(fp_det, "\tmax.lengthL =  %d\n", *max_lengthL);
#endif

	// return(list(L, max.lengthL))
	return ris;
}

SEXP createRules(SEXP m, SEXP f_pr_and)
{
	int i, nProtected = 0, lengthL;
	VETTOREi *op = NULL, *tmp = NULL, *max_lengthL = NULL;
	MATRICEd *m1;
	muParserHandle_t hParser;
	GString *f_pr_and1, *errore;
	SEXP ret_lista, ret_lista1, mL, lr;

	_InitDbg(false, false, false);

	_Intestazione("\n*** createRules1 ***\n");

	// matrice
	m1 = inMATRICE_d(m, &nProtected);
	f_pr_and1 = inSTRINGA(f_pr_and, &nProtected, "formula");
	CREAv_i(max_lengthL, 1);
	// alloca la lista per il risultato
	PROTECT(ret_lista1 = allocVector(VECSXP,  LENGTHm1_d(m1)));
	++nProtected;

	// max.lengthL<-0
	ASSEGNAv_i(max_lengthL, 1, 0);
	// L <- list()
	PROTECT(ret_lista = allocVector(VECSXP,  2));
	++nProtected;
	if (f_pr_and1->len > 0) {
		muParserHandle_t hParser = InitCalc();
#ifdef FDEBUG
		fprintf(fp_fdbg, "Compilo la funzione f_pr_and definita come '%s' ... ", f_pr_and1->str);
#endif
		mupSetExpr(hParser, f_pr_and1->str);
#ifdef FDEBUG
		fprintf(fp_fdbg, "ok\n");
#endif
	}
	// for(i in 1:n){
	for (i = 1; i <= LENGTHm1_d(m1); i++) {
		// g_op <- which(M[i, ] != 0)
		// l_op <- length(g_op)
		CREAv_i(tmp, 1);
		op = which_m_rowindxne_d(op, m1, i, 0.0);
		// if (l_op>1)
		if (LENGTHv_i(g_op) > 1) {
			tmp = create_logicRule(tmp, op, (f_pr_and1->len > 0), hParser);
			lengthL = LENGTHv_i(tmp);
		}
		else {
			if (LENGTHv_i(op) == 1) {
				tmp = copia_v_i(tmp, op, 1, LENGTHv_i(g_op)); // qui la copio perche' ho comunque due definizioni alternative
				lengthL = 1;
			}
			else {
				// logic_rule <- NaN
				lengthL = 0;
				tmp = NULL;
			}
		}
		//~ lengthL = length(logic_rule);
		//~ max.lengthL<-max(max.lengthL,lengthL)
		ASSEGNAv_i(max_lengthL, 1, max_s_i(ACCEDIv_i(max_lengthL, 1), lengthL));
	// L[[i]] <- logic_rule
		lr = NULL;
		if (lengthL > 0)
			lr = daVETTORE_i(tmp, &nProtected);
		else
			CANCELLAv_i(tmp);
		SET_LENGTH(lr, lengthL);
		SET_VECTOR_ELT(ret_lista1, i - 1, lr);
	}
	// return(list(L, max.lengthL))
	// assegna gli elementi
	if (f_pr_and1->len > 0)
			DeInitCalc(hParser);
	SET_VECTOR_ELT(ret_lista, 0, ret_lista1);
	mL = daVETTORE_i(max_lengthL, &nProtected);
	SET_LENGTH(mL, 1);
	SET_VECTOR_ELT(ret_lista, 1, mL);
	UNPROTECT(nProtected);

	CANCELLAv_i(op);
	CANCELLAm_d(m1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ret_lista;
}
