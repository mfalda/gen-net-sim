#include "connetti_scalefree.h"

#define g_ind_out globali.connetti_scalefree.ind_out
#define g_co_out globali.connetti_scalefree.co_out
#define g_Sc2_out globali.connetti_scalefree.Sc2_out
#define g_Sc2_in globali.connetti_scalefree.Sc2_in
#define g_toll globali.connetti_scalefree.toll
#define g_Sc1_out globali.connetti_scalefree.Sc1_out
#define g_non_connessi globali.connetti_scalefree.non_connessi
#define g_p_out globali.connetti_scalefree.p_out
#define g_Sc1_in globali.connetti_scalefree.Sc1_in
#define g_co_in globali.connetti_scalefree.co_in
#define g_Freq_in globali.connetti_scalefree.Freq_in
#define g_Freq_in1 globali.connetti_scalefree.Freq_in1
#define g_Sc_out globali.connetti_scalefree.Sc_out
#define g_ind_in globali.connetti_scalefree.ind_in
#define g_Freq_out globali.connetti_scalefree.Freq_out
#define g_Freq_out1 globali.connetti_scalefree.Freq_out1
#define g_Sin globali.connetti_scalefree.Sin
#define g_Sout globali.connetti_scalefree.Sout
#define g_tmp1_i globali.connetti_scalefree.tmp1_i
#define g_tmp2_i globali.connetti_scalefree.tmp2_i
#define g_tmp1_d globali.connetti_scalefree.tmp1_d
#define g_connessi globali.connetti_scalefree.connessi
#define g_nonco_out globali.connetti_scalefree.nonco_out
#define g_p globali.connetti_scalefree.p
#define g_p_in globali.connetti_scalefree.p_in
#define g_Sc_in globali.connetti_scalefree.Sc_in
#define g_nonco_in globali.connetti_scalefree.nonco_in

// connetti.scalefree<-function(Mdiscr,STout,STin,dist,g_toll,max.con,und=FALSE)
MATRICEi *connetti_scalefree(MATRICEi *ris, const VETTOREd *STout, const VETTOREd *STin, const VETTOREd *dist, const VETTOREd *toll1, int max_con, bool und/*=FALSE*/)
{
	int choice, m, N;
	double p1, p2;

	_Intestazione("\n*** connetti_scalefree ***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tMdiscr = ");
	_StampaRawMatr_i(ris);
	fprintf(fp_det, "\tSTout = ");
	_StampaRawVett_d(STout);
	fprintf(fp_det, "\tSTin = ");
	_StampaRawVett_d(STin);
	fprintf(fp_det, "\tdist = ");
	_StampaRawVett_d(dist);
	fprintf(fp_det, "\ttoll1 = ");
	_StampaRawVett_d(toll1);
	fprintf(fp_det, "\tmax_con = %d\n", max_con);
	if (und)
		fprintf(fp_det, "\tund = TRUE\n");
	else
		fprintf(fp_det, "\tund = FALSE\n");
#endif

	// g_toll verra` modificato, quindi lo copio, mentre Mdiscr, anche se modificata, serve per il risultato
	g_toll = copia_v_d(g_toll, toll1, 1, LENGTHv_d(toll1));
	//	non.g_connessi<-which(dist==Inf)
	g_non_connessi = which_v_indxeq_d(g_non_connessi, dist, R_PosInf);
//   g_connessi<-which(dist!=Inf)
	g_connessi = which_v_indxne_d(g_connessi, dist, R_PosInf);
//   N<-length(dist)
	N = LENGTHv_d(dist);
//   g_Sout<-apply(Mdiscr,2,sum)
	g_Sout = somma_colonne_i(g_Sout, ris);
//   m<-max(g_Sout)+1
	m = max_v_i(g_Sout) + 1;
//   Freq.out<-hist(g_Sout,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
	g_tmp1_i = seq_i(g_tmp1_i, 0, m, 1);
	g_Freq_out1 = hist1(g_Freq_out1, g_Sout, g_tmp1_i, 0, 1, 0);
	g_Freq_out = promuovi_i(g_Freq_out, g_Freq_out1);
//   Sc.out<-Score(S=g_Sout,ST=STout,Freq=Freq.out,n=1,g_toll)
	g_Sc_out = score1(g_Sc_out, g_Sout, STout, g_Freq_out, 1, g_toll);
//   if (und==TRUE)
	if (und) {
		//   Sc.out[which(g_Sout==max.con)]<-(-Inf)
		g_tmp1_i = which_v_indxeq_i(g_tmp1_i, g_Sout, max_con);
		assegna1_v_indx_d(g_Sc_out, g_tmp1_i, R_NegInf);
	}
//   g_Sin<-apply(Mdiscr,1,sum)
	g_Sin = somma_righe_i(g_Sin, ris);
//   m<-max(g_Sin)+1
	m = max_v_i(g_Sin) + 1;
//   Freq.in<-hist(g_Sin,breaks=seq(0,m,1),right=FALSE,plot=FALSE)$counts
	g_tmp1_i = seq_i(g_tmp1_i, 0, m, 1);
	g_Freq_in1 = hist1(g_Freq_in1, g_Sin, g_tmp1_i, 0, 1, 0);
	g_Freq_in = promuovi_i(g_Freq_in, g_Freq_in1);
//   Sc.in<-Score(S=g_Sin,ST=STin,Freq=Freq.in,n=1,g_toll)
	g_Sc_in = score1(g_Sc_in, g_Sin, STin, g_Freq_in, 1, g_toll);
//   Sc.in[which(g_Sin==max.con)]<-(-Inf)
	g_tmp1_i = which_v_indxeq_i(g_tmp1_i, g_Sin, max_con);
	assegna1_v_indx_d(g_Sc_in, g_tmp1_i, R_NegInf);
//   co.in<-intersect(g_connessi,which(Sc.in!=(-Inf)))
	g_tmp2_i = which_v_indxne_d(g_tmp2_i, g_Sc_in, R_NegInf);
	g_co_in = interseca_i(g_co_in, g_connessi, g_tmp2_i);
//   nonco.in<-intersect(non.g_connessi,which(Sc.in!=(-Inf)))
	g_nonco_in = interseca_i(g_nonco_in, g_non_connessi, g_tmp2_i);
//   co.out<-intersect(g_connessi,which(Sc.out!=(-Inf)))
	g_tmp1_i = which_v_indxne_d(g_tmp1_i, g_Sc_out, R_NegInf);
	g_co_out = interseca_i(g_co_out, g_connessi, g_tmp1_i);
//   nonco.out<-intersect(non.g_connessi,which(Sc.out!=(-Inf)))
	g_nonco_out = interseca_i(g_nonco_out, g_non_connessi, g_tmp1_i);
//   while ( ((length(nonco.in)==0)|(length(co.out)==0))&((length(co.in)==0)|(length(nonco.out)==0)) )
	while (((LENGTHv_i(g_nonco_in) == 0) || (LENGTHv_i(g_co_out) == 0)) && ((LENGTHv_i(g_co_in) == 0) || (LENGTHv_i(g_nonco_out) == 0))) {
//     {cat("\n WARNING: tolerance parameters were relaxed to allow connectivity of the graph\n")
		warning("tolerance parameters were relaxed to allow connectivity of the graph\n");
//      g_toll<-g_toll+1
		incr1_v_d(g_toll, 1.0);
		g_tmp1_d = rep_s_d(g_tmp1_d, R_PosInf, N);
//      Sc.out<-Score(S=g_Sout,ST=STout,Freq=Freq.out,n=1,g_toll=rep(Inf,N))
		g_Sc_out = score1(g_Sc_out, g_Sout, STout, g_Freq_out, 1, g_tmp1_d);
		g_tmp1_i = which_v_indxne_d(g_tmp1_i, g_Sc_out, R_NegInf);
//      Sc.in<-Score(S=g_Sin,ST=STin,Freq=Freq.in,n=1,g_toll=rep(Inf,N))
		g_Sc_in = score1(g_Sc_in, g_Sin, STin, g_Freq_in, 1, g_tmp1_d);
		g_tmp2_i = which_v_indxne_d(g_tmp2_i, g_Sc_in, R_NegInf);
//      co.in<-intersect(g_connessi,which(Sc.in!=(-Inf)))
		g_co_in = interseca_i(g_co_in, g_connessi, g_tmp2_i);
//      nonco.in<-intersect(non.g_connessi,which(Sc.in!=(-Inf)))
		g_nonco_in = interseca_i(g_nonco_in, g_non_connessi, g_tmp2_i);
//      co.out<-intersect(g_connessi,which(Sc.out!=(-Inf)))
		g_co_out = interseca_i(g_co_out, g_connessi, g_tmp1_i);
//      nonco.out<-intersect(non.g_connessi,which(Sc.out!=(-Inf)))
		g_nonco_out = interseca_i(g_nonco_out, g_non_connessi, g_tmp1_i);
//     }
	}
//   choice<-0
	choice = 0;
//   Sc1.in<-Sc.in[nonco.in]
	g_Sc1_in = copia_v_indx_d(g_Sc1_in, g_Sc_in, g_nonco_in);
//   Sc1.out<-Sc.out[co.out]
	g_Sc1_out = copia_v_indx_d(g_Sc1_out, g_Sc_out, g_co_out);
//   Sc2.in<-Sc.in[co.in]
	g_Sc2_in = copia_v_indx_d(g_Sc2_in, g_Sc_in, g_co_in);
//   Sc2.out<-Sc.out[nonco.out]
	g_Sc2_out = copia_v_indx_d(g_Sc2_out, g_Sc_out, g_nonco_out);
//   if ( (length(nonco.in)!=0)&(length(co.out)!=0)&(length(co.in)!=0)&(length(nonco.out)!=0) )
	if ((LENGTHv_i(g_nonco_in) != 0) && (LENGTHv_i(g_co_out) != 0) && (LENGTHv_i(g_co_in) != 0) && (LENGTHv_i(g_nonco_out) != 0)) {
// 	{p1<-mean(Sc1.in)/2+mean(Sc1.out)/2
		p1 = media_v_d(g_Sc1_in) / 2.0 + media_v_d(g_Sc1_out) / 2.0;
//          p2<-mean(Sc2.in)/2+mean(Sc2.out)/2
		p2 = media_v_d(g_Sc2_in) / 2.0 + media_v_d(g_Sc2_out) / 2.0;
//          g_p<-c(p1,p2)
		g_p = vettore2s_d(g_p, p1, p2);
//          if (length(which(g_p<=0))>0)
		g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_p, 0.0);
		if (LENGTHv_i(g_tmp1_i) > 0)
			// g_p<-g_p-min(g_p)+1/(N^2)
			somma1_vs_d(g_p, 1 / (N * N) - min_v_d(g_p));
// 	 g_p<-g_p/sum(g_p)
		dividi1_vs_d(g_p, somma_v_d(g_p, false));
//          choice<-sampleB(c(1,2),1,prob=g_p)
		g_tmp1_i = vettore2s_i(g_tmp1_i, 1, 2);
		g_tmp2_i = sampleB_p(g_tmp2_i, g_tmp1_i, 1, false, g_p);
		choice = ACCEDIv_i(g_tmp2_i, 1);
// 	}
	}
//   if ( ((length(nonco.in)==0)|(length(co.out)==0)) | (choice==2) )
	if (((LENGTHv_i(g_nonco_in) == 0) || (LENGTHv_i(g_co_out) == 0)) || (choice == 2)) {
// 	{if (length(which(Sc2.in<=0))>0) Sc2.in<-Sc2.in-min(Sc2.in)+1/(N^2)
		g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc2_in, 0.0);
		if (LENGTHv_i(g_tmp1_i) > 0)
			dividi1_vs_d(g_Sc2_in, 1 / (N * N) - min_v_d(g_Sc2_in));
// 	 g_p.in<-Sc2.in/sum(Sc2.in)
		g_p_in = dividi_vs_d(g_p_in, g_Sc2_in, somma_v_d(g_Sc2_in, false));
// 	 ind.in<-sampleB(co.in,1,prob=g_p.in)
		g_ind_in = sampleB_p(g_ind_in, g_co_in, 1, false, g_p_in);
//          if (length(which(Sc2.out<=0))>0) Sc2.out<-Sc2.out-min(Sc2.out)+1/(N^2)
		g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc2_out, 0.0);
		if (LENGTHv_i(g_tmp1_i) > 0)
			somma1_vs_d(g_Sc2_out, 1 / (N * N) - min_v_d(g_Sc2_out));
// 	 g_p.out<-Sc2.out/sum(Sc2.out)
		g_p_out = dividi_vs_d(g_p_out, g_Sc2_out, somma_v_d(g_Sc2_out, false));
// 	 ind.out<-sampleB(nonco.out,1,prob=g_p.out)
		g_ind_out = sampleB_p(g_ind_out, g_nonco_out, 1, false, g_p_out);
// 	}
	}
//   if ( ((length(co.in)==0)|(length(nonco.out)==0)) | (choice==1) )
	if ( ((LENGTHv_i(g_co_in) == 0) || (LENGTHv_i(g_nonco_out) == 0)) || (choice == 1) )
// 	{if (length(which(Sc1.in<=0))>0)
	{
		// Sc1.in<-Sc1.in-min(Sc1.in)+1/(N^2)
		g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc1_in, 0.0);
		if (LENGTHv_i(g_tmp1_i) > 0)
			somma1_vs_d(g_Sc1_in, 1 / (N * N) - min_v_d(g_Sc1_in));
// 	 g_p.in<-Sc1.in/sum(Sc1.in)
		g_p_in = dividi_vs_d(g_p_in, g_Sc1_in, somma_v_d(g_Sc1_in, false));
// 	 ind.in<-sampleB(nonco.in,1,prob=g_p.in)
		g_ind_in = sampleB_p(g_ind_in, g_nonco_in, 1, false, g_p_in);
//          if (length(which(Sc1.out<=0))>0) Sc1.out<-Sc1.out-min(Sc1.out)+1/(N^2)
		g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc1_out, 0.0);
		if (LENGTHv_i(g_tmp1_i) > 0)
			somma1_vs_d(g_Sc1_out, 1 / (N * N) - min_v_d(g_Sc1_out));
// 	 g_p.out<-Sc1.out/sum(Sc1.out)
		g_p_out = dividi_vs_d(g_p_out, g_Sc1_out, somma_v_d(g_Sc1_out, false));
// 	 ind.out<-sampleB(co.out,1,prob=g_p.out)
		g_ind_out = sampleB_p(g_ind_out, g_co_out, 1, false, g_p_out);
// 	}
	}
//   Mdiscr[ind.in,ind.out]<-1
	assegna1_m_vv_i(ris, g_ind_in, g_ind_out, 1);
//   if (und==TRUE)  Mdiscr[ind.out,ind.in]<-1
	if (und)
		assegna1_m_vv_i(ris, g_ind_out, g_ind_in, 1);

	StrBilanciam();

	_Intestazione("\n*** Esco da connetti_scalefree ***\n\n");

#ifdef DET
	fprintf(fp_det, "connetti_scalefree output:\n\tris = ");
	_StampaRawMatr_i(ris);
#endif

//  return(Mdiscr)
	return ris;
// }
}
