#include "create_logicRule.h"

#define g_x globali.create_logicRule.x
#define g_s globali.create_logicRule.s
#define g_scalare_d globali.create_logicRule.scalare_d
#define g_nvect globali.create_logicRule.nvect
#define g_e globali.create_logicRule.e
#define g_prob globali.create_logicRule.prob
#define g_tmp_i globali.create_logicRule.tmp_i
#define g_tmp1_i globali.create_logicRule.tmp1_i
#define g_pr_and globali.create_logicRule.pr_and
#define g_pr_or globali.create_logicRule.pr_or
#define g_scalare_i globali.create_logicRule.scalare_i
#define g_blacklist globali.create_logicRule.blacklist
#define g_black_p globali.create_logicRule.black_p
#define g_o globali.create_logicRule.o

/*
VETTOREd *f(int n, VETTOREd *norm_lev)
{
	VETTOREd *ris = NULL;
	VETTOREd *g_scalare_d = NULL;

	if (n == 0) {
		CREAv_d(g_scalare_d, 1);
		ASSEGNAv_d(g_scalare_d, 1, 0.8);
		ris = rep_d(ris, g_scalare_d, LENGTHv_d(norm_lev));
		CANCELLAv_d(g_scalare_d);
	}
	else {
		error("funzione per createRules non definita!\n");
		return NULL;
	}
	return ris;
}*/

// create.logicrule<-function(v, f.pr.and=NULL){
VETTOREi *create_logicRule(VETTOREi *ris, const VETTOREi *v, int f_pr_and, muParserHandle_t hParser)
{
	int livello, n, p, L_p;
	int i, j, i1, l, k, lr, lo;
	double args[1], fVal;
	GString *errore;
#ifdef MDEBUG
	GString *tmp = NULL;
#endif

	_Intestazione("\n***create_logicRule***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tv = ");
	_StampaRawVett_i(v);
	fprintf(fp_det, "\tf.pr.and = ...\n");
#endif

	CREAv_i(g_scalare_i, 1);
	ASSEGNAv_i(g_scalare_i, 1, 0);
	CREAv_d(g_scalare_d, 1);
	// l <- length(v)
	l = LENGTHv_i(v);
	// g_x<-(1:l)/(l)
	g_x = op_ss_seqdiv_d(g_x, 1, l, (double) l);
	// if (is.null(f.pr.and)) pr.and<-rep(0.5,length(g_x))
	if (!f_pr_and) {
		CREAv_d(g_pr_and, LENGTHv_d(g_x));
		InitVett_d(g_pr_and, 0.5);
	}
	// else pr.and<-f.pr.and(g_x)
	else {
		CREAv_d(g_pr_and, LENGTHv_d(g_x));
		mupDefineVar(hParser, "x", &args[0]);
		for (i = 0; i < LENGTHv_d(g_x); i++) {
			args[0] = ACCEDIv_d(g_x, i + 1);
#ifdef FDEBUG
	fprintf(fp_fdbg, "Calcolo la funzione f_pr_and in %.5e: ", args[0]);
#endif
			fVal = mupEval(hParser);
			if (!mupError(hParser))
				g_pr_and->dati[i] = fVal;
#ifdef FDEBUG
	fprintf(fp_fdbg, "%.5e\n", g_pr_and->dati[i]);
#endif
		}
	}
	// pr.or<-1-pr.and
	g_pr_or = complementa_d(g_pr_or, g_pr_and);

	//ris<-c(sample(c(-2,-3),1, g_prob=c(pr.and[1],pr.or[1])),0,0)
	g_prob = vettore2s_d(g_prob, ACCEDIv_d(g_pr_and, 1), ACCEDIv_d(g_pr_or, 1));
	g_tmp_i	= vettore2s_i(g_tmp_i, -2, -3);
	g_s = sample_p(g_s, g_tmp_i, 1, 0, g_prob, "create_logicRule");
	ASSEGNAv_i(g_scalare_i, 1, 0);
	ris = vettore3v_i(ris, g_s, g_scalare_i, g_scalare_i);
	// g_o<-c(2,3)
	g_o = vettore2s_i(g_o, 2, 3);
	// if (ris[1]==(-2)) {
	if (ACCEDIv_i(ris, 1) == -2) {
		// pr.or<-rep(0,length(pr.or))
		CREAv_d(g_pr_or, LENGTHv_d(g_pr_or));
		InitVett_d(g_pr_or, 0.0);
		// pr.and<-rep(1,length(pr.and))
		CREAv_d(g_pr_and, LENGTHv_d(g_pr_and));
		InitVett_d(g_pr_and, 1.0);
	}
	// black.p<-g_blacklist<-vector()
	CREAv_i(g_black_p, 0);
	CREAv_i(g_blacklist, 0);
	// if (l>2){
	if (l > 2) {
		// 	for(i in 1:(l-2)){
		for (i = 1; i <= l - 2; i++) {
			// 		lr<-length(ris)
			// 		lo<-length(g_o)
			// 		g_e<-sample(seq(1,lo,1),1)
			lr = LENGTHv_i(ris);
			lo = LENGTHv_i(g_o);
			g_tmp_i = seq_i(g_tmp_i, 1, lo, 1);
			g_e = sample(g_e, g_tmp_i, 1, 0, "create_logicRule");
			// 		p<-g_o[g_e]
			p = ACCEDIv_i(g_o, ACCEDIv_i(g_e, 1));
			// 		if (2*p>lr){
			if (2 * p > lr) {
				// 			g_nvect<-rep(-1,lr)
				CREAv_i(g_nvect, lr);
				InitVett_i(g_nvect, -1);
				// 			ris<-c(ris,g_nvect)
				ris = accoda1_vv_i(ris, g_nvect);
				// 		}
			}
			// 		g_o[g_e]<-2*p
			ASSEGNAv_i(g_o, ACCEDIv_i(g_e, 1), 2 * p);
			// 		g_o<-c(g_o,(2*p+1))
			g_o = accoda1_vs_i(g_o, 2 * p + 1);
			// 		L.p<-length(black.p)
			L_p = LENGTHv_i(g_black_p);
			// if (L.p>0) { ;
			if (L_p > 0) {
				// for (j in (1:L.p)) {
				for (j = 1; j <= L_p; j++) {
					//g_blacklist <- vector()
					CREAv_i(g_blacklist, 0);
					// g_x<-ceiling(log(max(g_o)/black.p[j],2))
#ifdef MDEBUG
				if (ACCEDIv_i(g_black_p, j) == 0) {
					CREAstr(tmp, "");
					g_string_printf(tmp, "ATTENZIONE (create_logicRule.c, linea 141): divisione per zero!\n");
					warning(tmp->str);
					fprintf(fp_fdbg, tmp->str);
					CANCELLAstr(tmp);
				}
#endif
					n = (int) ceil(log2(max_v_i(g_o) / ACCEDIv_i(g_black_p, j))) ;
					// for (i in (1:g_x))
					for (i1 = 0; i1 < n; i1++) {
						//g_blacklist<-c(g_blacklist,black.p[j]*(2^i)+(0:((2^i)-1)))
						g_tmp_i	= seq_i(g_tmp_i, 0, (2 << i1) - 1, 1);
						somma1_vs_i(g_tmp_i, ACCEDIv_i(g_black_p, j) * (2 << i1));
						g_blacklist = accoda1_vv_i(g_blacklist, g_tmp_i);
					}
				}
			}
			// 		livello<-floor(log(p,2))+1
			livello = (int) floor(log2(p));
			// 		if (p%in%g_blacklist)
			if (esiste_v_i(p, g_blacklist) > 0)
				//  ris[p]<-(-2)
				ASSEGNAv_i(ris, p, -2);
			//     		else {
			else {
				// ris[p]<-sample(c(-2,-3),1,g_prob=c(pr.and[livello],pr.or[livello]))
				g_tmp_i = vettore2s_i(g_tmp_i, -2, -3);
				g_prob = vettore2s_d(g_prob, ACCEDIv_d(g_pr_and, livello), ACCEDIv_d(g_pr_or, livello));
				g_tmp1_i = sample_p(g_tmp1_i, g_tmp_i, 1, 0, g_prob, "create_logicRule");
				ASSEGNAv_i(ris, p, ACCEDIv_i(g_tmp1_i, 1));
			}

			// 		if (ris[p]==(-2))
			if (ACCEDIv_i(ris, p) ==  -2) {
				// black.p<-c(black.p,p)
				g_black_p = accoda1_vs_i(g_black_p, p);
			}

			//     ris[2*p]<-0
			// usare ASSEGNAv_i se non corretto
			ris = assegna_v_i(ris, 2 * p, 0);
			// 		ris[2*p+1]<-0
			// usare ASSEGNAv_i se non corretto
			ris = assegna_v_i(ris, 2 * p + 1, 0);
			// 	}
		}
		// }
	}
	// for(i in 1:length(g_o)){
	for (i = 1; i <= LENGTHv_i(g_o); i++) {
		// 	ris[g_o[i]]<-v[i]
		ASSEGNAv_i(ris, ACCEDIv_i(g_o, i), ACCEDIv_i(v, i));
		// }
	}
	// k<-length(ris)
	k = LENGTHv_i(ris);
	// repeat{
	while (1) {
		// 	if (ris[k]!=-1) break
		if (ACCEDIv_i(ris, k) != -1)
			break;
		// 	ris<-ris[-k]
		elimina1_indx_i(ris, k);
		// 	k<-k-1
		k--;
		// }
	}

	//~ CANCELLAv_d(g_x);
	//~ CANCELLAv_i(g_s);
	//~ CANCELLAv_d(g_scalare_d);
	//~ CANCELLAv_i(g_nvect);
	//~ CANCELLAv_i(g_e);
	//~ CANCELLAv_d(g_prob);
	//~ CANCELLAv_i(g_tmp_i);
	//~ CANCELLAv_i(g_tmp1_i);
	//~ CANCELLAv_d(g_pr_and);
	//~ CANCELLAv_d(g_pr_or);
	//~ CANCELLAv_i(g_scalare_i);
	//~ CANCELLAv_i(g_blacklist);
	//~ CANCELLAv_i(g_black_p);
	//~ CANCELLAv_i(g_o);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "create_logicRule output:\n");
	fprintf(fp_det, "\tr = ");
	_StampaRawVett_i(ris);
#endif

	// return(ris)
	return ris;
	// }
}
