#include "connectivity_geometric.h"

#define g_M globali.connectivity_geometric.M
#define g_Mdiscr globali.connectivity_geometric.Mdiscr
#define g_x globali.connectivity_geometric.x
#define g_y globali.connectivity_geometric.y
#define g_xy globali.connectivity_geometric.xy
#define g_d globali.connectivity_geometric.d
#define g_s globali.connectivity_geometric.s
#define g_regulatedind globali.connectivity_geometric.regulatedind
#define g_indL globali.connectivity_geometric.indL
#define g_Sr globali.connectivity_geometric.Sr
#define g_aus globali.connectivity_geometric.aus
#define g_ind globali.connectivity_geometric.ind
#define g_ind1 globali.connectivity_geometric.ind1
#define g_ind0 globali.connectivity_geometric.ind0
#define g_tmp1_i globali.connectivity_geometric.tmp1_i
#define g_tmp2_i globali.connectivity_geometric.tmp2_i
#define g_tmp1_d globali.connectivity_geometric.tmp1_d

// connectivitygeometric<-function(N=50,k=3, weight.mean=1, weight.sd=0.1)
LISTA *connectivity_geometric1(LISTA *ris, int N, double k, double weight_mean, double weight_sd)
{
	int i, L, L1, ri, n;
	double r;
	enum TIPO tipi[2];

	_Intestazione("\n***connectivity_geometric***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tN =  %d\n", N);
	fprintf(fp_det, "\tk =  %.16g\n", k);
	fprintf(fp_det, "\tweight_mean =  %.16g\n", weight_mean);
	fprintf(fp_det, "\tweight_sd =  %.16g\n", weight_sd);
#endif

	GetRNGstate();
	//  Mdiscr<-matrix(0,ncol=N,nrow=N)
	CREAm_i(g_Mdiscr, N, N);
	InitMatr_i(g_Mdiscr, 0);
	//  r<-sqrt(k/(N*pi))*1.1*sqrt(2)
	r = sqrt((double) k / (N * M_PI)) * 1.1 * sqrt(2.0);
	//  g_x<-runif_s(N)
	CREAv_d(g_x, N);
	CREAv_d(g_y, N);
	CREAv_d(g_xy, 1);
	g_x = runif_s(g_x, N, 0.0, 1.0, "cg");
	//  g_y<-runif_s(N)
	g_y = runif_s(g_y, N, 0.0, 1.0, "cg");
	//  for (i in (1:N))
	for (i = 1; i <= N; i++) {
		//   {g_d<-sqrt((g_x[i]-g_x)^2+(g_y[i]-g_y)^2)
		g_d = distanza_2dvs_d(g_d, g_x, g_y, ACCEDIv_d(g_x, i), ACCEDIv_d(g_y, i));
		//    g_ind<-setdiff(which((g_d<r)&(g_d>0)),seq(0,(i-1),1))
		g_tmp1_i = which_v_andglt_d(g_tmp1_i, g_d, r, 0.0);
		g_tmp2_i = seq_i(g_tmp2_i, 0, (i - 1), 1);
		g_ind = setdiff_i(g_ind, g_tmp1_i, g_tmp2_i);
		//    g_s<-sampleB(c(0,1),length(g_ind),replace=TRUE)
		g_tmp1_i = vettore2s_i(g_tmp1_i, 0, 1);
		g_s = sampleB(g_s, g_tmp1_i, LENGTHv_i(g_ind), 1);
		//    g_ind1<-which(g_s==1)
		g_ind1 = which_v_indxeq_i(g_ind1, g_s, 1);
		//    if (length(g_ind1)>0) Mdiscr[i,g_ind[g_ind1]]<-1
		if (LENGTHv_i(g_ind1) > 0) {
			g_tmp1_i = assegna_v_indx_i(g_tmp1_i, g_ind, g_ind1);
			assegna1_ms_rigaindx_i(g_Mdiscr, i, g_tmp1_i, 1);
		}
		//    g_ind0<-which(g_s==0)
		g_ind0 = which_v_indxeq_i(g_ind0, g_s, 0);
		//    if (length(g_ind0)>0)  Mdiscr[g_ind[g_ind0],i]<-1
		if (LENGTHv_i(g_ind0) > 0) {
			g_tmp1_i = assegna_v_indx_i(g_tmp1_i, g_ind, g_ind0);
			assegna1_ms_indxcol_i(g_Mdiscr, g_tmp1_i, i, 1);
		}
	}
	// ###check that every gene has at least 1 regulator
	//##check that every gene has at least 1 regulator;
	// g_Sr<-apply(Mdiscr,1,sum)
	g_Sr = somma_righe_i(g_Sr, g_Mdiscr);
	// g_indL<-which(g_Sr==0)
	g_indL = which_v_indxeq_i(g_indL, g_Sr, 0);
	// L<-length(g_indL)
	L = LENGTHv_i(g_indL);
	L1 = L;
	n = 0;
	// while (L>0)
	while (L > 0) {
		//  {ri<-g_indL[1]
		ri = ACCEDIv_i(g_indL, 1);
		//   g_regulatedind<-which(Mdiscr[,ri]==1)
		g_regulatedind = which_m_colindxeq_i(g_regulatedind, g_Mdiscr, ri, 1);
		//   if (length(g_regulatedind)==0)
		if (LENGTHv_i(g_regulatedind) == 0) {
			//     {#cat("\n k is low with respect to N. It is suggested to increase k")
		 //cat("\n k is low with respect to N_ It is suggested to increase k");
			//      g_x[ri]<-runif_s(1)
			g_xy = runif_s(g_xy, 1, 0.0, 1.0, "cg");
			ASSEGNAv_d(g_x, ri, ACCEDIv_d(g_xy, 1));
			//      g_y[ri]<-runif_s(1)
			g_xy = runif_s(g_xy, 1, 0.0, 1.0, "cg");
			ASSEGNAv_d(g_y, ri, ACCEDIv_d(g_xy, 1));
			//      g_d<-sqrt((g_x[ri]-g_x)^2+(g_y[ri]-g_y)^2)
			g_d = distanza_2dvs_d(g_d, g_x, g_y, ACCEDIv_d(g_x, ri),ACCEDIv_d(g_y, ri));
			//      g_ind<-which((g_d<r)&(g_d>0))
			g_ind = which_v_andglt_d(g_ind, g_d, r, 0.0);
			//      g_s<-sampleB(c(0,1),length(g_ind),replace=TRUE)
			g_tmp1_i = vettore2s_i(g_tmp1_i, 0, 1);
			g_s = sampleB(g_s, g_tmp1_i, LENGTHv_i(g_ind), 1);
			//      g_ind1<-which(g_s==1)
			g_ind1 = which_v_indxeq_i(g_ind1, g_s, 1);
			//      if (length(g_ind1)>0) Mdiscr[ri,g_ind[g_ind1]]<-1
			if (LENGTHv_i(g_ind1) > 0) {
				g_tmp1_i = assegna_v_indx_i(g_tmp1_i, g_ind, g_ind1);
				assegna1_ms_rigaindx_i(g_Mdiscr, ri, g_tmp1_i, 1);
			}
			//      g_ind0<-which(g_s==0)
			g_ind0 = which_v_indxeq_i(g_ind0, g_s, 0);
			//      if (length(g_ind0)>0)  Mdiscr[g_ind[g_ind0],ri]<-1
			if (LENGTHv_i(g_ind0) > 0) {
				g_tmp1_i = assegna_v_indx_i(g_tmp1_i, g_ind, g_ind0);
				assegna1_ms_indxcol_i(g_Mdiscr, g_tmp1_i, ri, 1);
			}
		}
		//   else {g_s<-sampleB(g_regulatedind,1);  Mdiscr[ri,g_s]<-1}
		else {
			g_s = sampleB(g_s, g_regulatedind, 1, 0);
			ASSEGNAm_i(g_Mdiscr, ri, ACCEDIv_i(g_s, 1), 1);
		}
		//   g_Sr<-apply(Mdiscr,1,sum)
		g_Sr = somma_righe_i(g_Sr, g_Mdiscr);
		//   g_indL<-which(g_Sr==0)
		g_indL = which_v_indxeq_i(g_indL, g_Sr, 0);
		//   L<-length(g_indL)
		L = LENGTHv_i(g_indL);
		if (n > 2 * N && L == L1) {
			warning("the algorithm does not converge: exiting!\n");
			L = 0;
		}
		L1 = L;
		n++;
		//  }
	}
	// #######assign weights to the connectivity matrix##########
	//######assign weights to the connectivity matrix##########;
	// g_ind<-which(Mdiscr==1,arr.g_ind=TRUE)
	g_ind = which_m_indxeq_i(g_ind, g_Mdiscr, 1);
	// L<-dim(g_ind)[1]
	L = LENGTHv_i(g_ind);
	// g_aus<-abs(rnorm_s(L,weight.mean, weight.sd))
	CREAv_d(g_tmp1_d, L);
	g_tmp1_d = rnorm_s(g_tmp1_d, L, weight_mean, weight_sd, "cg");
	g_aus = abs_v_d(g_aus, g_tmp1_d);
	// g_M<-matrix(0,ncol=N,nrow=N)
	CREAm_d(g_M, N, N);
	InitMatr_d(g_M, 0.0);
	// g_M[g_ind]<-aus
	assegna1_mv_indx_d(g_M, g_ind, g_aus);
	PutRNGstate();
	// return(list(g_M,Mdiscr))
	tipi[0] = MATRd;
	tipi[1] = MATRi;
	CreaLISTA(ris, tipi, 2);
	ris->dati[0].md = g_M;
	ris->dati[1].mi = g_Mdiscr;
	//~ CANCELLAv_d(g_x);
	//~ CANCELLAv_d(g_y);
	//~ CANCELLAv_d(g_d);
	//~ CANCELLAv_i(g_s);
	//~ CANCELLAv_i(g_regulatedind);
	//~ CANCELLAv_i(g_indL);
	//~ CANCELLAv_i(g_Sr);
	//~ CANCELLAv_d(g_aus);
	//~ CANCELLAv_i(g_ind);
	//~ CANCELLAv_i(g_ind1);
	//~ CANCELLAv_i(g_ind0);
	//~ CANCELLAv_i(g_tmp1_i);
	//~ CANCELLAv_i(g_tmp2_i);
	//~ CANCELLAv_d(g_tmp1_d);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "cg output:\n");
	fprintf(fp_det, "\tlist(M,Mdiscr) = ");
	_StampaRawMatr_d(g_M);
	_StampaRawMatr_i(g_Mdiscr);
#endif

	return ris;
}

SEXP connectivity_geometric(SEXP N, SEXP k, SEXP weight_mean, SEXP weight_sd)
{
	int nProtected = 0;
	int N1;
	double k1, weight_mean1, weight_sd1;
	LISTA *l = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** connectivity_geometric ***\n");

	N1 = INTEGER_VALUE(N);
	k1 = NUMERIC_VALUE(k);
	weight_mean1 = NUMERIC_VALUE(weight_mean);
	weight_sd1 = NUMERIC_VALUE(weight_sd);

	l = connectivity_geometric1(l, N1, k1 ,weight_mean1, weight_sd1);
	ris = daLISTA(l, &nProtected);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
