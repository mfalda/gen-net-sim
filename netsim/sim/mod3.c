#include "mod3.h"

#define g_M_out globali.mod3.M_out
#define g_indok globali.mod3.indok
#define g_tmp3_i globali.mod3.tmp3_i
#define g_tmp1_d globali.mod3.tmp1_d
#define g_tmpSTin globali.mod3.tmpSTin
#define g_tmpSTout globali.mod3.tmpSTout
#define g_Freq_in globali.mod3.Freq_in
#define g_Freq_out globali.mod3.Freq_out
#define g_tmp3_d globali.mod3.tmp3_d
#define g_p globali.mod3.p
#define g_Sc globali.mod3.Sc
#define g_tmp2_i globali.mod3.tmp2_i
#define g_ind_M globali.mod3.ind_M
#define g_ind1 globali.mod3.ind1
#define g_indS globali.mod3.indS
#define g_indBS globali.mod3.indBS
#define g_Sin globali.mod3.Sin
#define g_ind globali.mod3.ind
#define g_indInf globali.mod3.indInf
#define g_tmp1_i globali.mod3.tmp1_i
#define g_scalare_i globali.mod3.scalare_i
#define g_scalare_d globali.mod3.scalare_d
#define g_p_sc globali.mod3.p_sc

// MOD3<-function(Ng.out,Ng.in,STin,STout,max.con,Cg=0)
MATRICEi *mod3(MATRICEi *ris, int Ng_out, int Ng_in, const VETTOREd *STin, const VETTOREd *STout, int max_con, double Cg)
{
	int N, j, n_reg, L, Ls, ns;
	double S, tmpd;
	int aus;
#ifdef MDEBUG
	GString *tmp = NULL;
#endif

	_Intestazione("\n***mod3***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tNg_out =  %d\n", Ng_out);
	fprintf(fp_det, "\tNg_in =  %d\n", Ng_in);
	fprintf(fp_det, "\tSTin = ");
	_StampaRawVett_d(STin);
	fprintf(fp_det, "\tSTout = ");
	_StampaRawVett_d(STout);
	fprintf(fp_det, "\tmax_con =  %d\n", max_con);
	fprintf(fp_det, "\tCg =  %.16g\n", Cg);
#endif

/*
#ifdef FDEBUG
	fprintf(fp_fdbg, "Ng_out = %d\n", Ng_out);
	fprintf(fp_fdbg, "Ng_in = %d\n", Ng_in);
	_StampaVett_d(STin);
	_StampaVett_d(STout);
	fprintf(fp_fdbg, "max_con = %d\nCg = %.5e\n", max_con, Cg);
	fprintf(fp_fdbg, "\n");
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif*/

	CREAv_d(g_scalare_d, 1);
	CREAv_i(g_scalare_i, 1);
	// ris<-matrix(0,Ng.out+Ng.in,Ng.out+Ng.in)
	CREAm_i(ris, Ng_out + Ng_in, Ng_out + Ng_in);
	InitMatr_i(ris, 0);
	// N<-Ng.in+Ng.out
	N = Ng_in + Ng_out;
	// S<-sum(STout[2:(Ng.in+1)],na.rm=TRUE)
	g_tmp1_d = segmento_v_d(g_tmp1_d, STout, 2, Ng_in + 1);
	S = somma_v_d(g_tmp1_d, true);
	// if (S>0) g_p<-STout[2:(Ng.in+1)]/S
	if (S > 0.0) {
		g_p = segmento_v_d(g_p, STout, 2, Ng_in + 1);
		dividi1_vs_d(g_p, S);
	}
	// else g_p<-rep(1/Ng.in,Ng.in)
	else {
#ifdef MDEBUG
		if (Ng_in == 0) {
			CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (mod3.c, linea 53): divisione per zero!\n");
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			CANCELLAstr(tmp);
		}
#endif
		CREAv_d(g_p, Ng_in);
		InitVett_d(g_p, (double) 1.0 / Ng_in);
	}
	// M.out<-sampleB(seq(1,Ng.in,1),Ng.out,prob=g_p,replace=TRUE)
	g_tmp1_i = seq_i(g_tmp1_i, 1, Ng_in, 1);
	g_M_out = sampleB_p(g_M_out, g_tmp1_i, Ng_out, 1, g_p);
	g_tmpSTin = copia_v_d(g_tmpSTin, STin, 1, LENGTHv_d(STin));
	g_tmpSTout = copia_v_d(g_tmpSTout, STout, 1, LENGTHv_d(STout));
	// while (sum(M.out)<(N-1)) # rm.na == FALSE
	while (somma_v_i(g_M_out, false) < N - 1) {
		//   {if (S>0)
		if (S > 0.0) {
			// 	{STout<-(STout[1:(Ng.in+1)]/S)*Ng.in
			segmento1_v_d(g_tmpSTout, 1, Ng_in + 1);
#ifdef MDEBUG
		if (Ng_in == 0) {
			CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (mod3.c, linea 74): divisione per zero!\n");
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			CANCELLAstr(tmp);
		}
#endif
			dividi1_vs_d(g_tmpSTout, (double) S / Ng_in);
	//          aus<-max(M.out)+1
			aus = max_v_i(g_M_out) + 1;
			// Freq.out<-hist(M.out,breaks=seq(0,Ng.in+1,1),right=FALSE,plot=FALSE)$counts/N
			g_tmp2_i = seq_i(g_tmp2_i, 0, Ng_in + 1, 1);
			g_tmp1_i = hist1(g_tmp1_i, g_M_out, g_tmp2_i, 0, 1, 0);
			g_Freq_out = dividi_vs_i(g_Freq_out, g_tmp1_i, (double) N);
			//          g_Sc<-Score(S=M.out,ST=STout,Freq=Freq.out,n=1,toll=rep(Inf,(Ng.out+1)))
			CREAv_d(g_tmp3_d, Ng_out + 1);
			InitVett_d(g_tmp3_d, R_PosInf);
			g_Sc = score1(g_Sc, g_M_out, g_tmpSTout, g_Freq_out, 1, g_tmp3_d);
			//          g_indok<-which(g_Sc!=-Inf)
			g_indok = which_v_indxne_d(g_indok, g_Sc, R_NegInf);
			// 	 g_indInf <- setdiff(seq(1, Ng.out, 1), g_indok)
			g_tmp3_i = seq_i(g_tmp3_i, 1, Ng_out, 1);
			g_indInf = setdiff_i(g_indInf, g_tmp3_i, g_indok);
			// 	 g_Sc[g_indInf]<-min(c(0,g_Sc[g_indok]))-1
			g_tmp1_d = assegna_v_indx_d(g_tmp1_d, g_Sc, g_indok);
			tmpd = min_v_d(g_tmp1_d);
			if (tmpd < 0.0)
				tmpd = 0.0;
			assegna1_v_indx_d(g_Sc, g_indInf, tmpd - 1);
			//          g_p<-g_Sc/sum(g_Sc)
			g_p = dividi_vs_d(g_p, g_Sc, somma_v_d(g_Sc, false));
			// 	}
		}
		//    else g_p<-rep(1/Ng.in,Ng.in)
		else {
#ifdef MDEBUG
		if (Ng_in == 0) {
			CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (mod3.c, linea 111): divisione per zero!\n");
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			CANCELLAstr(tmp);
		}
#endif
			CREAv_d(g_p, Ng_in);
			InitVett_d(g_p, (double) 1.0 / Ng_in);
		}
		//    g_ind.M<-sampleB(seq(1,Ng.out,1),1,prob=g_p)
		g_tmp2_i = seq_i(g_tmp2_i, 1, Ng_out, 1);
		g_ind_M = sampleB_p(g_ind_M, g_tmp2_i, 1, 0, g_p);
		//    M.out[g_ind.M]<-M.out[g_ind.M]+1
		ASSEGNAv_i(g_M_out, ACCEDIv_i(g_ind_M, 1), ACCEDIv_i(g_M_out, ACCEDIv_i(g_ind_M, 1)) + 1);
	}
	// for (j in (1:Ng.out))
	for (j = 1; j <= Ng_out; j++) {
		//   {n.reg<-M.out[j]
		n_reg = ACCEDIv_i(g_M_out, j);
		//    if (j==1)
		if (j == 1) {
			//      {ris[1:n.reg,(j+Ng.in)]<-1
			g_tmp1_i = seq_i(g_tmp1_i, 1, n_reg, 1);
			ASSEGNAv_i(g_scalare_i, 1, j + Ng_in);
			assegna1_m_vv_i(ris, g_tmp1_i, g_scalare_i, 1);
			//       g_indS<-1:n.reg
			g_indS = seq_i(g_indS, 1, n_reg, 1);
			//       g_indBS<-setdiff(1:Ng.in,g_indS)
			g_tmp2_i = seq_i(g_tmp2_i, 1, Ng_in, 1);
			g_indBS = setdiff_i(g_indBS, g_tmp2_i, g_indS);
			//      }
		}
		//     else
		else {
			//      {L<-length(g_indBS)
			L = LENGTHv_i(g_indBS);
			//       Ls<-min(L,n.reg-1)
			Ls = min_s_i(L, n_reg - 1);
			//       ns<-n.reg-Ls
			ns = n_reg - Ls;
			//       S<-sum(STin[1:(Ng.out+1)],na.rm=TRUE)
			g_tmp1_d = segmento_v_d(g_tmp1_d, g_tmpSTin, 1, Ng_out + 1);
			S = somma_v_d(g_tmp1_d, true);
			//       g_Sin<-apply(ris,1,sum)
			g_Sin = somma_righe_i(g_Sin, ris);
			//       if (S>0)
			if (S > 0.0) {
				// 	{STin<-(STin[1:(Ng.out+1)]/S)*Ng.out
				segmento1_v_d(g_tmpSTin, 1, Ng_out + 1);
#ifdef MDEBUG
		if (Ng_out == 0) {
			CREAstr(tmp, "");
			g_string_printf(tmp, "ATTENZIONE (mod3.c, linea 163): divisione per zero!\n");
			warning(tmp->str);
			fprintf(fp_fdbg, tmp->str);
			CANCELLAstr(tmp);
		}
#endif
				dividi1_vs_d(g_tmpSTin, (double) S / Ng_out);
				//          aus<-max(g_Sin)+1
				aus = max_v_i(g_Sin) + 1;
				//          Freq.in<-hist(g_Sin,breaks=seq(0,aus,1),right=FALSE,plot=FALSE)$counts/N
				g_tmp2_i = seq_i(g_tmp2_i, 0, aus, 1);
				g_tmp1_i = hist1(g_tmp1_i, g_Sin, g_tmp2_i, 0, 1, 0);
				g_Freq_in = dividi_vs_i(g_Freq_in, g_tmp1_i, (double) N);
				//          g_Sc<-Score(S=g_Sin[g_indS],ST=STin,Freq=Freq.in,n=1,toll=rep(Inf,(Ng.out+1)))
				g_tmp2_i = assegna_v_indx_i(g_tmp2_i, g_Sin, g_indS);
				CREAv_d(g_tmp3_d, Ng_out + 1);
				InitVett_d(g_tmp3_d, R_PosInf);
				g_Sc = score1(g_Sc, g_tmp2_i, g_tmpSTin, g_Freq_in, 1, g_tmp3_d);
				//          g_indok<-which(g_Sc!=-Inf)
				g_indok = which_v_indxne_d(g_indok, g_Sc, R_NegInf) ;
				// 	 g_indInf <- setdiff(seq(1, length(g_indS), 1), g_indok)
				g_tmp1_i = seq_i(g_tmp1_i, 1, LENGTHv_i(g_indS), 1);
				g_indInf = setdiff_i(g_indInf, g_tmp1_i, g_indok);
				// 	 g_Sc[g_indInf]<-min(c(0,g_Sc[g_indok]))-1
				g_tmp1_d = assegna_v_indx_d(g_tmp1_d, g_Sc, g_indok);
				tmpd = min_v_d(g_tmp1_d);
				if (tmpd < 0.0)
					tmpd = 0.0;
				assegna1_v_indx_d(g_Sc, g_indInf, tmpd - 1);
				//          g_p.sc<-g_Sc/sum(g_Sc)
				g_p_sc = dividi_vs_d(g_p_sc, g_Sc, somma_v_d(g_Sc, false));
				// 	}
			}
			//        else g_p.sc<-rep(1/length(g_indS),length(g_indS))
			else {
#ifdef MDEBUG
				if (LENGTHv_i(g_indS) == 0) {
					CREAstr(tmp, "");
					g_string_printf(tmp, "ATTENZIONE (mod3.c, linea 201): divisione per zero!\n");
					warning(tmp->str);
					fprintf(fp_fdbg, tmp->str);
					CANCELLAstr(tmp);
				}
#endif
				CREAv_d(g_p_sc, LENGTHv_i(g_indS));
				InitVett_d(g_p_sc, (double) 1.0 / LENGTHv_i(g_indS));
			}
			//
			//        g_ind<-which(g_Sin[g_indS]==max.con)

			g_ind = which_v_indxeq_i(g_ind, g_tmp1_i, max_con);
			//        g_p.sc[g_ind]<-0
			assegna1_vs_indx_d(g_p_sc, g_ind, 0.0);
			//        g_p.sc<-g_p.sc/sum(g_p.sc)
			g_p = dividi_vs_d(g_p, g_p_sc, somma_v_d(g_p_sc, false));
			//       g_ind1<-sampleB(g_indS,ns,prob=g_p.sc)
			g_ind1 = sampleB_p(g_ind1, g_indS, ns, 0, g_p_sc);
			//       if (Ls>0)
			if (Ls > 0) {
				// 	{g_ind<-g_indBS[1:Ls]
				g_ind = segmento_v_i(g_ind, g_indBS, 1, Ls);
				//          g_indS<-c(g_indS,g_ind)
				g_indS = accoda1_vv_i(g_indS, g_ind);
				//          g_indBS<-setdiff(1:Ng.in,g_indS)
				g_tmp1_i = seq_i(g_tmp1_i, 1, Ng_in, 1);
				g_indBS = setdiff_i(g_indBS, g_tmp1_i, g_indS);
				//          g_ind1<-c(g_ind1,g_ind)
				g_ind1 = accoda1_vv_i(g_ind1, g_ind);
				// 	}
			}
			//       m[g_ind1,(j+Ng.in)]<-1
			ASSEGNAv_i(g_scalare_i, 1, j + Ng_in);
			assegna1_m_vv_i(ris, g_ind1, g_scalare_i, 1);
			//      }
		}
		//   }
	}
	// if ( (Cg>0)&((Ng.in+Ng.out)>2) ) ris<-triangola(M=ris,Cg=Cg,max.con)
	if (Cg > 0.0 && Ng_in + Ng_out > 2)
		ris = triangola(ris, Cg, max_con);
	// return(ris)
	//~ CANCELLAv_i(g_M_out);
	//~ CANCELLAv_i(g_indok);
	//~ CANCELLAv_i(g_tmp3_i);
	//~ CANCELLAv_d(g_tmp1_d);
	//~ CANCELLAv_d(g_Freq_out);
	//~ CANCELLAv_d(g_tmp3_d);
	//~ CANCELLAv_d(g_p);
	//~ CANCELLAv_d(g_Sc);
	//~ CANCELLAv_i(g_tmp2_i);
	//~ CANCELLAv_i(g_ind_M);
	//~ CANCELLAv_i(g_ind1);
	//~ CANCELLAv_i(g_indS);
	//~ CANCELLAv_i(g_indBS);
	//~ CANCELLAv_i(g_Sin);
	//~ CANCELLAv_i(g_ind);
	//~ CANCELLAv_i(g_indInf);
	//~ CANCELLAv_i(g_tmp1_i);
	//~ CANCELLAv_i(g_scalare_i);
	//~ CANCELLAv_d(g_scalare_d);
	//~ CANCELLAv_d(g_p_sc);
#ifdef FDEBUG
	fprintf(fp_fdbg, "\n\n");
	fprintf(fp_fdbg, "*****************************************\n");
	fprintf(fp_fdbg, "*****************************************\n\n");
#endif

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "mod3 output:\n");
	fprintf(fp_det, "\tm = ");
	_StampaRawMatr_i(ris);
#endif

	return ris;
	// }
}
