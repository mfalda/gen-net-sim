#include "scoremodular.h"


#define g_tmpm_d globali.scoremodular.tmpm_d
#define g_scalare_d globali.scoremodular.scalare_d
#define g_S_new globali.scoremodular.S_new
#define g_tmp_i1 globali.scoremodular.tmp_i1
#define g_tmp1_d globali.scoremodular.tmp1_d
#define g_T1 globali.scoremodular.T1
#define g_T2 globali.scoremodular.T2
#define g_old1 globali.scoremodular.old1
#define g_old2 globali.scoremodular.old2
#define g_new1 globali.scoremodular.new1
#define g_new2 globali.scoremodular.new2
#define g_toll1 globali.scoremodular.toll1
#define g_toll2 globali.scoremodular.toll2
#define g_a globali.scoremodular.a
#define g_b globali.scoremodular.b
#define g_ind globali.scoremodular.ind
#define g_m globali.scoremodular.m
#define g_S1 globali.scoremodular.S1
#define g_S2 globali.scoremodular.S2
#define g_ind1 globali.scoremodular.ind1
#define g_ind2 globali.scoremodular.ind2
#define g_indinf globali.scoremodular.indinf

// Score<-function(S,ST,Freq,n,toll)
VETTOREd *score1(VETTOREd *Sc, const VETTOREi *S, const VETTOREd *ST, const VETTOREd *Freq, int n, const VETTOREd *toll)
{
	_Intestazione("\n***score***\n");
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
	g_S_new = somma_vs_i(g_S_new, S, n);
	//  g_T1<-ST[S+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, S, 1);
	// usare assegna_vindx_d se non corretto
	g_T1 = assegna_v_indxNA_d(g_T1, ST, g_tmp_i1);
	//  g_T2<-ST[S.new+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, g_S_new, 1);
	g_T2 = assegna_v_indxNA_d(g_T2, ST, g_tmp_i1);
	//  OLD1<-Freq[S+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, S, 1);
	g_old1 = assegna_v_indxNA_d(g_old1, Freq, g_tmp_i1);
	//  OLD2<-Freq[S.new+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, g_S_new, 1);
	g_old2 = assegna_v_indxNA_d(g_old2, Freq, g_tmp_i1);
	//  NEW1<-OLD1-1
	g_new1 = somma_vs_d(g_new1, g_old1, -1.0);
	//  NEW2<-OLD2+1
	g_new2 = somma_vs_d(g_new2, g_old2, 1.0);
	//  g_toll1<-toll[S+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, S, 1);
	// usare assegna_vindx_d se non corretto
	g_toll1 = assegna_v_indxNA_d(g_toll1, toll, g_tmp_i1);
	//  g_toll2<-toll[S.new+1]
	g_tmp_i1 = somma_vs_i(g_tmp_i1, g_S_new, 1);
	// usare assegna_vindx_d se non corretto
	g_toll2 = assegna_v_indxNA_d(g_toll2, toll, g_tmp_i1);
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
		/* ESEMPIO:
			g_m = cbind(seq(1,5),seq(1,3))
			warning: le righe del ris. non sono multiple
			g_m = [ [1, 1], [2,2], [3, 3], [4, 1], [5, 2] ]
			poi calcola il massimo di ciascuna riga
		*/
		g_tmpm_d = cbind2v_d(g_tmpm_d, g_a, g_b);
		g_m = max_righe_d(g_m, g_tmpm_d);
		//    g_S1<-(sign(g_a-g_b))*g_m/T1
		g_S1 = f_aux_d(g_S1, g_a, g_b, g_m, g_T1);
		//    g_ind1<-which((g_a-g_toll1)<0)
		g_tmp1_d = diff_vv_d(g_tmp1_d, g_a, g_toll1);
		g_ind1 = which_v_indxlt_d(g_ind1, g_tmp1_d, 0.0);
		//    g_ind2<-which((g_b-g_toll1)>0)
		g_tmp1_d = diff_vv_d(g_tmp1_d, g_b, g_toll1);
		g_ind2 = which_v_indxgt_d(g_ind2, g_tmp1_d, 0.0);
		//    g_indinf<-intersect(g_ind1,g_ind2)
		g_indinf = interseca_i(g_indinf, g_ind1, g_ind2);
		//    g_S1[g_indinf]<-(-Inf)
		assegna1_vs_indx_d(g_S1, g_indinf, R_NegInf);
		//    g_S1[g_ind]<-0
		assegna1_vs_indx_d(g_S1, g_ind, 0.0);
		//	g_ind<-which(is.na(g_S1))                          #added 19 09 2008   g_a questo punto g_S1==NA se g_T1==0
		g_ind = which_v_indxNA_d(g_ind, g_S1, false);
		//	g_S1[g_ind]<-Inf*sign(g_a[g_ind]-g_b[g_ind])                         #added 19 09 2008
		g_tmp1_d = rep_s_d(g_tmp1_d, R_PosInf, LENGTHv_d(g_a));
		f_aux10_d(g_S1, g_ind, g_a, g_b, g_tmp1_d);
		//	g_ind<-which(g_S1==Inf)                            #added 19 09 2008
		g_ind = which_v_indxeq_d(g_ind, g_S1, R_PosInf);
		//	g_S1[g_ind]<-(sign(g_a[g_ind]-g_b[g_ind]))*g_m[g_ind]/1 // ERRORE?
		f_aux10_d(g_S1, g_ind, g_a, g_b, g_m);
		//   }
	}
	//  else g_S1<-rep(0,length(g_a))
	else {
		CREAv_d(g_S1, LENGTHv_d(g_a));
		InitVett_d(g_S1, 0.0);
	}
	//  g_a<-abs(OLD2-g_T2)
	g_a = diff_vv_d(g_a, g_old2, g_T2);
	abs1_v_d(g_a);
	//  g_b<-abs(NEW2-g_T2)
	g_b = diff_vv_d(g_b, g_new2, g_T2);
	abs1_v_d(g_b);
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
		g_S2 = f_aux_d(g_S2, g_a, g_b, g_m ,g_T2);
		//    g_ind1<-which((g_a-g_toll2)<0)
		g_tmp1_d = diff_vv_d(g_tmp1_d, g_a, g_toll2);
		g_ind1 = which_v_indxlt_d(g_ind1, g_tmp1_d, 0.0);
		//    g_ind2<-which((g_b-g_toll2)>0)
		g_tmp1_d = diff_vv_d(g_tmp1_d, g_b, g_toll2);
		g_ind2 = which_v_indxgt_d(g_ind2, g_tmp1_d, 0.0);
		//    g_indinf<-intersect(g_ind1,g_ind2)
		g_indinf = interseca_i(g_indinf, g_ind1, g_ind2);
		//    g_S2[g_indinf]<-(-Inf)
		assegna1_vs_indx_d(g_S2, g_indinf, R_NegInf);
		//    g_S2[g_ind]<-0
		assegna1_vs_indx_d(g_S2, g_ind, 0.0);
		//	g_ind<-which(is.na(g_S2))                          #added 19 09 2008  g_a questo punto  g_S2==NA se g_T2==0
		g_ind = which_v_indxNA_d(g_ind, g_S2, false);
		//	g_S2[g_ind]<-Inf*sign(g_a[g_ind]-g_b[g_ind])                         #added 19 09 2008
		g_tmp1_d = rep_s_d(g_tmp1_d, R_PosInf, LENGTHv_d(g_a));
		f_aux10_d(g_S2, g_ind, g_a, g_b, g_tmp1_d);
		//	g_ind<-which(g_S2==Inf)                            #added 19 09 2008
		g_ind = which_v_indxeq_d(g_ind, g_S2, R_PosInf);
		//	g_S2[g_ind]<-(sign(g_a[g_ind]-g_b[g_ind]))*g_m[g_ind]/1        #added 19 09 2008
		g_tmp1_d = rep_s_d(g_tmp1_d, R_PosInf, LENGTHv_d(g_a));
		f_aux10_d(g_S2, g_ind, g_a, g_b, g_m);
		//   }
	}
	//  else g_S2<-rep(0,length(g_a))
	else {
		CREAv_d(g_S2, LENGTHv_d(g_a));
		InitVett_d(g_S2, 0.0);
	}
	//  Sc<-(g_S1+g_S2)/2
	Sc = somma_vv_d(Sc, g_S1, g_S2);
	dividi1_vs_d(Sc, 2.0);
	//  return(Sc)

	StrBilanciam();

	_Intestazione("\n*** Esco da scoremodular ***\n\n");

#ifdef DET
	fprintf(fp_det, "scoremodular output:\n\tSc = ");
	_StampaRawVett_d(Sc);
#endif

	return Sc;
	// }
}
