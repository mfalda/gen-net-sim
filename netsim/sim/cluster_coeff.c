#include "cluster_coeff.h"

#define g_tmpm_i globali.cluster_coeff.tmpm_i
#define g_ind globali.cluster_coeff.ind

VETTOREd *cluster_coeff2(VETTOREd *ris, const MATRICEi *W, double *coeff)
{
	int i, n, l, j, ind_j, ind_h, h, ng;
	double Kg;

	_Intestazione("\n***cluster_coeff***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tW = ");
	_StampaRawMatr_i(W);
#endif

	// W<-W+t(W)
	g_tmpm_i = trasponi_i(g_tmpm_i, W);
	somma1_m_i(g_tmpm_i, W);
	// g_ind<-which(W!=0,arr.g_ind=TRUE)
	g_ind = which_m_indxne_i(g_ind, g_tmpm_i, 0);
	// W[g_ind]<-1  # matrice simmetrica con solo 0 e 1
	assegna1_ms_indx_i(g_tmpm_i, g_ind, 1);
	// diag(W)<-0
	assegna1_s_diag_i(g_tmpm_i, 0);
	// N<-dim(W)[1]
	n = LENGTHm1_i(g_tmpm_i);
	// ris<-rep(0,N)
	CREAv_d(ris, n);
	InitVett_d(ris, 0.0);
	// for (i in (1:N)) {
	for (i = 1; i <= n; i++) {
		// Kg <- sum(W[i,])
		Kg = somma_riga_i(g_tmpm_i, i);
		// neighbours<-which(W[i,]!=0)
		g_ind = which_m_rowindxne_i(g_ind, g_tmpm_i, i, 0);
		// ng<-0
		ng = 0;
		// L<-length(neighbours)
		l = LENGTHv_i(g_ind);
		// if (L>1) {
		if (l > 1) {
			// for (j in (1:L)) {
			for (j = 1; j <= l; j++) {
				// g_ind.j<-neighbours[j]
				ind_j = ACCEDIv_i(g_ind, j);
				// for (h in (j:L))
				for (h = j; h <= l; h++)	{
					// g_ind.h<-neighbours[h]
					ind_h =ACCEDIv_i(g_ind, h);
					// if (W[g_ind.j,g_ind.h]==1)
					if (ACCEDIm_i(g_tmpm_i, ind_j, ind_h) == 1)
						// ng<-ng+1
						ng++;
				}
			}
			// Cg[i]<-2*ng/(Kg*(Kg-1))
			ASSEGNAv_d(ris, i, (double) 2.0 * ng / (Kg * (Kg - 1)));
		}
	}

	//~ CANCELLAm_i(g_tmpm_i);
	//~ CANCELLAv_i(g_ind);

	// coeff<-mean(Cg,na.rm=TRUE)
	*coeff = media_v_d(ris);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "cluster_coeff output:\n");
	fprintf(fp_det, "\tlist(coeff,Cg) =  %.16g\n", *coeff);
	_StampaRawVett_d(ris);
#endif

	// return(list(coeff,Cg))
	return ris;
}

LISTA *cluster_coeff1(LISTA *ris, const MATRICEi *W)
{
	MATRICEi *tmpm_i = NULL;
	VETTOREi *ind = NULL, *neighbours = NULL;
	VETTOREd *Cg = NULL, *scalare_d = NULL;
	int i, n, l, j, ind_j, ind_h, h, ng, Kg;
	double coeff;
	enum TIPO tipi[2];

	_Intestazione("\n*** cluster_coeff1 ***\n");

	CREAv_d(scalare_d, 1);
	// W<-W+t(W)
	tmpm_i = trasponi_i(tmpm_i, W);
	somma1_m_i(tmpm_i, W);
	// ind<-which(W!=0,arr.ind=TRUE)
	ind = which_m_indxne_i(ind, tmpm_i, 0);
	// W[ind]<-1  # matrice simmetrica con solo 0 e 1
	assegna1_ms_indx_i(tmpm_i, ind, 1);
	// diag(W)<-0
	assegna1_s_diag_i(tmpm_i, 0);
	// N<-dim(W)[1]
	n = LENGTHm1_i(tmpm_i);
	// Cg<-rep(0,N)
	CREAv_d(Cg, n);
	InitVett_d(Cg, 0.0);
	// for (i in (1:N)) {
	for (i = 1; i <= n; i++) {
		// Kg <- sum(W[i,])
		Kg = somma_riga_i(tmpm_i, i);
		// neighbours<-which(W[i,]!=0)
		neighbours = which_m_rowindxne_i(neighbours, tmpm_i, i, 0);
		// ng<-0
		ng = 0;
		// L<-length(neighbours)
		l = LENGTHv_i(neighbours);
		// if (L>1) {
		if (l > 1) {
			// for (j in (1:L)) {
			for (j = 1; j <= l; j++) {
				// ind.j<-neighbours[j]
				ind_j = ACCEDIv_i(neighbours, j);
				// for (h in (j:L))
				for (h = j; h <= l; h++)	{
					// ind.h<-neighbours[h]
					ind_h =ACCEDIv_i(neighbours, h);
					// if (W[ind.j,ind.h]==1)
					if (ACCEDIm_i(tmpm_i, ind_j, ind_h) == 1)
						// ng<-ng+1
						ng++;
				}
			}
			// Cg[i]<-2*ng/(Kg*(Kg-1))
			ASSEGNAv_d(Cg, i, (double) 2.0 * ng / (Kg * (Kg - 1)));
		}
	}

	CANCELLAm_i(tmpm_i);
	CANCELLAv_i(ind);
	if (n > 0)
		CANCELLAv_i(neighbours);

	// coeff<-mean(Cg,na.rm=TRUE)
	coeff = media_v_d(Cg);
	// return(list(coeff,Cg))

	tipi[0] = REALE;
	tipi[1] = VETTd;
	CreaLISTA(ris, tipi, 2);
	ASSEGNAv_d(scalare_d, 1, coeff);
	CtrlSlst(ris, 1);
	ASSEGNAlst(ris, 1, vd, scalare_d);
	CtrlSlst(ris, 2);
	ASSEGNAlst(ris, 2, vd, Cg);

	StrBilanciam();

	return ris;
}

SEXP cluster_coeff(SEXP W)
{
	int nProtected = 0;
	MATRICEi *w1;
	LISTA *l = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** cluster_coeff ***\n");

	w1 = inMATRICE_i(W, &nProtected);

	l = cluster_coeff1(l, w1);

	ris = daLISTA(l, &nProtected);

	CANCELLAm_i(w1);

	UNPROTECT(nProtected);

	controllaCanc_i();
	controllaCanc_d();

	return ris;
}
