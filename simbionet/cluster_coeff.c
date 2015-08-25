#include "cluster_coeff.h"

#define g_tmpm_i globali.cluster_coeff.tmpm_i
#define g_ind globali.cluster_coeff.ind
#define g_neighbours globali.cluster_coeff.neighbours

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
	// ind<-which(W!=0,arr.ind=TRUE)
	g_ind = which_m_indxne_i(g_ind, g_tmpm_i, 0);
	// W[ind]<-1  # matrice simmetrica con solo 0 e 1
	assegna1_ms_indx_i(g_tmpm_i, g_ind, 1);
	// diag(W)<-0
	assegna1_s_diag_i(g_tmpm_i, 0);
	// N<-dim(W)[1]
	n = LENGTHm1_i(g_tmpm_i);
	// Cg<-rep(0,N)
	CREAv_d(ris, n);
	InitVett_d(ris, 0.0);
	// for (i in (1:N)) {
	for (i = 1; i <= n; i++) {
		// Kg <- sum(W[i,])
		Kg = somma_riga_i(g_tmpm_i, i);
		// neighbours<-which(W[i,]!=0)
		g_neighbours = which_m_rowindxne_i(g_neighbours, g_tmpm_i, i, 0);
		// ng<-0
		ng = 0;
		// L<-length(neighbours)
		l = LENGTHv_i(g_neighbours);
		// if (L>1) {
		if (l > 1) {
			// for (j in (1:L)) {
			for (j = 1; j <= l; j++) {
				// ind.j<-neighbours[j]
				ind_j = ACCEDIv_i(g_neighbours, j);
				// for (h in (j:L))
				for (h = j; h <= l; h++)	{
					// ind.h<-neighbours[h]
					ind_h =ACCEDIv_i(g_neighbours, h);
					// if (W[ind.j,ind.h]==1)
					if (ACCEDIm_i(g_tmpm_i, ind_j, ind_h) == 1)
						// ng<-ng+1
						ng++;
				}
			}
			// Cg[i]<-2*ng/(Kg*(Kg-1))
			ASSEGNAv_d(ris, i, (double) 2.0 * ng / (Kg * (Kg - 1)));
		}
	}

	// coeff<-mean(Cg,na.rm=TRUE)
	*coeff = media_v_d(ris);

	StrBilanciam();

	_Intestazione("\n*** Esco da cluster_coeff2 ***\n");
#ifdef DET
	fprintf(fp_det, "cluster_coeff output:\n");
	fprintf(fp_det, "\tlist(coeff,Cg) =  %.7g\n", *coeff);
	_StampaRawVett_d(ris);
#endif

	// return(list(coeff,Cg))
	return ris;
}
