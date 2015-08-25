#include "score_sf.h"

#define g_tmpm_d globali.score_sf.tmpm_d
#define g_scalare_d globali.score_sf.scalare_d
#define g_tmp_i1 globali.score_sf.tmp_i1
#define g_tmp_d1 globali.score_sf.tmp_d1
#define g_T1 globali.score_sf.T1
#define g_T2 globali.score_sf.T2
#define g_old1 globali.score_sf.old1
#define g_old2 globali.score_sf.old2
#define g_new1 globali.score_sf.new1
#define g_new2 globali.score_sf.new2
#define g_toll1 globali.score_sf.toll1
#define g_toll2 globali.score_sf.toll2
#define g_a globali.score_sf.a
#define g_b globali.score_sf.b
#define g_ind globali.score_sf.ind
#define g_m globali.score_sf.m
#define g_S1 globali.score_sf.S1
#define g_S2 globali.score_sf.S2
#define g_ind1 globali.score_sf.ind1
#define g_ind2 globali.score_sf.ind2
#define g_indinf globali.score_sf.indinf
#define g_tmp_i2 globali.score_sf.tmp_i2
#define g_indbad globali.score_sf.indbad
#define g_ind0 globali.score_sf.ind0
#define g_ind3 globali.score_sf.ind3

// Score.sf<-function(S,ST,Freq,n,toll)
VETTOREd *score_sf1(VETTOREd *ris, const VETTOREi *S, const VETTOREd *ST, const VETTOREd *Freq, int n, const VETTOREd *toll)
{
	double somma;

	_Intestazione("\n***score_sf***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tS = ");
	_StampaRawVett_i(S);
	fprintf(fp_det, "\tST = ");
	_StampaRawVett_d(ST);
	fprintf(fp_det, "\tFreq = ");
	_StampaRawVett_d(Freq);
	fprintf(fp_det, "\tn =  %d\n", n);
	fprintf(fp_det, "\ttoll = ");
	_StampaRawVett_d(toll);
#endif

	CREAv_d(g_scalare_d, 1);
	//  S.new<-S+n
	g_ind = somma_vs_i(g_ind, S, n);
	//  g_T1<-ST[S+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, S, 1);
	g_T1 = assegna_v_indxNA_d(g_T1, ST, g_tmp_i1);
	//  g_T2<-ST[S.new+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, g_ind, 1);
	g_T2 = assegna_v_indxNA_d(g_T2, ST, g_tmp_i1);
	//  OLD1<-Freq[S+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, S, 1);
	g_old1 = assegna_v_indxNA_d(g_old1, Freq, g_tmp_i1);
	//  OLD2<-Freq[S.new+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, g_ind, 1);
	g_old2 = assegna_v_indxNA_d(g_old2, Freq, g_tmp_i1);
	//  NEW1<-OLD1-1
	g_new1 = somma_vs_d(g_new1, g_old1, -1.0);
	//  NEW2<-OLD2+1
	g_new2 = somma_vs_d(g_new2, g_old2, 1.0);
	//  g_toll1<-toll[S+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, S, 1);
	g_toll1 = assegna_v_indx_d(g_toll1, toll, g_tmp_i1);
	//  g_toll2<-toll[S.new+1]
	g_tmp_i2 = somma_vs_i(g_tmp_i2, g_ind, 1);
	g_toll2 = assegna_v_indx_d(g_toll2, toll, g_tmp_i2);
	//  g_a<-abs(OLD1-g_T1)
	g_a = diff_vv_d(g_a, g_old1, g_T1);
	abs1_v_d(g_a);
	//  g_b<-abs(NEW1-g_T1)
	g_b = diff_vv_d(g_b, g_new1, g_T1);
	abs1_v_d(g_b);
	//  g_ind<-which(is.na(g_T1))
	g_ind = which_v_indxNA_d(g_ind, g_T1, false);
	//  if (length(g_ind)<length(g_a))
	if (LENGTHv_i(g_ind) < LENGTHv_d(g_a)) {
		//   {g_a[g_ind]<-0
		assegna1_vs_indx_d(g_a, g_ind, 0.0);
		//    g_b[g_ind]<-0
		assegna1_vs_indx_d(g_b, g_ind, 0.0);
		//    g_m<-apply(cbind(g_a,g_b),1,max)
		g_tmpm_d = cbind2v_d(g_tmpm_d, g_a, g_b);
		g_m = max_righe_d(g_m, g_tmpm_d);
		//    g_S1<-(sign(g_a-g_b))*g_m/T1
		g_S1 = f_aux_d(g_S1, g_a, g_b, g_m, g_T1);
		//    g_indbad<-which((g_b-g_toll1)>0)
		g_tmp_d1 = diff_vv_d(g_tmp_d1, g_b, g_toll1);
		g_indbad = which_v_indxgt_d(g_indbad, g_tmp_d1, 0.0);
		//    if (length(g_indbad)>0)
		if (LENGTHv_i(g_indbad) > 0) {
			//     {g_ind1<-which((g_a-g_toll1)<=0)
			g_tmp_d1 = diff_vv_d(g_tmp_d1, g_a, g_toll1);
			g_ind1 = which_v_indxlt_d(g_ind1, g_tmp_d1, 0.0);
			//      g_ind2<-which( ((sign(OLD1-g_T1)-sign(NEW1-g_T1))!=0)&(sign(NEW1-g_T1)!=0)&(sign(OLD1-g_T1)!=0) )
			g_ind2 = f_aux2_d(g_ind2, g_old1, g_new1, g_T1);
			//      g_ind3<-which( ((g_a-g_toll1)>0)&(g_a<g_b) )
			g_ind3 = f_aux3_d(g_ind3, g_a, g_toll1, g_b);
			//      g_ind0<-union(g_ind3,union(g_ind1,g_ind2))
			g_ind1 = unione1_i(g_ind1, g_ind2);
			g_ind0 = unione_i(g_ind0, g_ind3, g_ind1);
			//      g_indinf<-union(intersect(g_ind0,g_indbad),which(NEW1<0))
			g_tmp_i1 = which_v_indxlt_d(g_tmp_i1, g_new1, 0.0);
			g_tmp_i2 = interseca_i(g_tmp_i2, g_ind0, g_indbad);
			g_indinf = unione_i(g_indinf, g_tmp_i2, g_tmp_i1);			//      g_S1[g_indinf]<-(-Inf)
			assegna1_vs_indx_d(g_S1, g_indinf, R_NegInf);
			//     }
		}
		//    g_S1[g_ind]<-0
		assegna1_vs_indx_d(g_S1, g_ind, 0.0);
		//   }
	}
	//  else g_S1<-rep(0,length(g_a))
	else {
		CREAv_d(g_S1, LENGTHv_d(g_a));
		InitVett_d(g_S1, 0.0);
	}
	//  g_S1[which(g_S1==Inf)]<-sum(ST,na.rm=TRUE)
	somma = somma_v_d(ST, true);
	assegna1_v_indxeq_d(g_S1, g_S1, R_PosInf, somma);
	//
	//  g_a<-abs(OLD2-g_T2)
	g_a = diff_vv_d(g_a, g_old2, g_T2);
	abs1_v_d(g_a);
	//  g_b<-abs(NEW2-g_T2)
	g_b = diff_vv_d(g_b, g_new2, g_T2);
	abs1_v_d(g_b);
	//  g_m<-apply(cbind(g_a,g_b),1,max)
	g_tmpm_d = cbind2v_d(g_tmpm_d, g_a, g_b);
	g_m = max_righe_d(g_m, g_tmpm_d);
	//  g_ind<-which(is.na(g_T2))
	g_ind = which_v_indxNA_d(g_ind, g_T2, false);
	//  if (length(g_ind)<length(g_a))
	if (LENGTHv_i(g_ind) < LENGTHv_d(g_a)) {
		//   {g_a[g_ind]<-0
		assegna1_vs_indx_d(g_a, g_ind, 0.0);
		//    g_b[g_ind]<-0
		assegna1_vs_indx_d(g_b, g_ind, 0.0);
		//    g_m<-apply(cbind(g_a,g_b),1,max)
		g_tmpm_d = cbind2v_d(g_tmpm_d, g_a, g_b);
		g_m = max_righe_d(g_m, g_tmpm_d);
		//    g_S2<-(sign(g_a-g_b))*g_m/T2
		g_S2 = f_aux_d(g_S2, g_a, g_b, g_m , g_T2);
		//    g_indbad<-which((g_b-g_toll2)>0)
		g_tmp_d1 = diff_vv_d(g_tmp_d1, g_b, g_toll2);
		g_indbad = which_v_indxgt_d(g_indbad, g_tmp_d1, 0.0);
		//    if (length(g_indbad)>0)
		if (LENGTHv_i(g_indbad) > 0) {
			//     {g_ind1<-which((g_a-g_toll2)<=0)
			g_tmp_d1 = diff_vv_d(g_tmp_d1, g_a, g_toll2);
			g_ind1 = which_v_indxle_d(g_ind1, g_tmp_d1, 0.0);
			//      g_ind2<-which( ((sign(OLD2-g_T2)-sign(NEW2-g_T2))!=0)&(sign(NEW2-g_T2)!=0)&(sign(OLD2-g_T2)!=0) )
			g_ind2 = f_aux2_d(g_ind2, g_old2, g_new2, g_T2);
			//      g_ind3<-which( ((g_a-g_toll2)>0)&(g_a<g_b) )
			g_ind3 = f_aux3_d(g_ind3, g_a, g_toll2, g_b);
			//      g_ind0<-union(g_ind3,union(g_ind1,g_ind2))
			g_ind1 = unione1_i(g_ind1, g_ind2);
			g_ind0 = unione_i(g_ind0, g_ind3, g_ind1);
			//      g_indinf<-doubleersect(g_ind0,g_indbad)
			g_indinf = interseca_i(g_indinf, g_ind0, g_indbad);
			//      g_S2[g_indinf]<-(-Inf)
			assegna1_vs_indx_d(g_S2, g_indinf, R_NegInf);
			//     }
		}
		//    g_S2[g_ind]<-0
		assegna1_vs_indx_d(g_S2, g_ind, 0.0);
		//   }
	}
	//  else g_S2<-rep(0,length(g_a))
	else {
		CREAv_d(g_S2, LENGTHv_d(g_a));
		InitVett_d(g_S2, 0.0);
	}
	//  g_ind<-which(g_T2==0)
	g_ind = which_v_indxeq_d(g_ind, g_T2, 0.0) ;
	//  g_S2[g_ind]<-(-Inf)
	assegna1_vs_indx_d(g_S2, g_ind, R_NegInf);
	//  Sc<-(g_S1+g_S2)/2
	ris = somma_vv_d(ris, g_S1, g_S2);
	dividi1_vs_d(ris, 2.0);

	//~ CANCELLAm_d(g_tmpm_d);
	//~ CANCELLAv_d(g_scalare_d);
	//~ CANCELLAv_i(g_tmp_i1);
	//~ CANCELLAv_d(g_tmp_d1);
	//~ CANCELLAv_d(g_T1);
	//~ CANCELLAv_d(g_T2);
	//~ CANCELLAv_d(g_old1);
	//~ CANCELLAv_d(g_old2);
	//~ CANCELLAv_d(g_new1);
	//~ CANCELLAv_d(g_new2);
	//~ CANCELLAv_d(g_toll1);
	//~ CANCELLAv_d(g_toll2);
	//~ CANCELLAv_d(g_a);
	//~ CANCELLAv_d(g_b);
	//~ CANCELLAv_i(g_ind);
	//~ CANCELLAv_d(g_m);
	//~ CANCELLAv_d(g_S1);
	//~ CANCELLAv_d(g_S2);
	//~ CANCELLAv_i(g_ind1);
	//~ CANCELLAv_i(g_ind2);
	//~ CANCELLAv_i(g_indinf);
	//~ CANCELLAv_i(g_tmp_i2);
	//~ CANCELLAv_i(g_indbad);
	//~ CANCELLAv_i(g_ind0);
	//~ CANCELLAv_i(g_ind3);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "score_sf output:\n\tSc = ");
	_StampaRawVett_d(ris);
#endif

	//  return(Sc)
	return ris;
	// }
}

SEXP score_sf(SEXP S, SEXP ST, SEXP Freq, SEXP n, SEXP toll)
{
	int n1, nProtected = 0;
	VETTOREi *S1 = NULL;
	VETTOREd *ST1 = NULL, *Freq1 = NULL, *toll1 = NULL, *ris1 = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** score_sf ***\n");

	S1 = inVETTORE_i(S, &nProtected);
	ST1 = inVETTORE_d(ST, &nProtected);
	Freq1 = inVETTORE_d(Freq, &nProtected);
	toll1 = inVETTORE_d(toll, &nProtected);
	//  int
	n1 = INTEGER_VALUE(n);

	ris1 = score_sf1(ris1, S1, ST1, Freq1, n1, toll1);

	ris = daVETTORE_d(ris1, &nProtected);

	CANCELLAv_i(S1);
	CANCELLAv_d(ST1);
	CANCELLAv_d(Freq1);
	CANCELLAv_d(toll1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
