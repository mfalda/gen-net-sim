#include "probmod_und.h"

#define g_Sc globali.probmod_und.Sc
#define g_scalare_d globali.probmod_und.scalare_d
#define g_checkIN globali.probmod_und.checkIN
#define g_checkOUT globali.probmod_und.checkOUT
#define g_memory globali.probmod_und.memory
#define g_M_in globali.probmod_und.M_in
#define g_M_out globali.probmod_und.M_out
#define g_tmp_i1 globali.probmod_und.tmp_i1
#define g_tmp_i2 globali.probmod_und.tmp_i2
#define g_indok globali.probmod_und.indok
#define g_indInf globali.probmod_und.indInf
#define g_I globali.probmod_und.I
#define g_ord_ind globali.probmod_und.ord_ind
#define g_aus globali.probmod_und.aus
#define g_rs globali.probmod_und.rs
#define g_ind1 globali.probmod_und.ind1
#define g_scalare_i globali.probmod_und.scalare_i
#define g_I_add globali.probmod_und.I_add
#define g_S_out globali.probmod_und.S_out
#define g_tmp_d1 globali.probmod_und.tmp_d1

MATRICEd *probmod2_und(MATRICEd *ris, MATRICEi *M, const VETTOREi *h, const VETTOREi *Sout, const VETTOREd *STout, const VETTOREd *Freq_out, const VETTOREd *toll, struct RisProbM *ris1)
{
	int nh, ng, i, j , L, ind;
	int n_out;
	char lab[3];
	double mn, Pm, somma;

	_Intestazione("\n***probmod_und***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tM = ");
	_StampaRawMatr_i(M);
	fprintf(fp_det, "\th = ");
	_StampaRawVett_i(h);
	fprintf(fp_det, "\tSout = ");
	_StampaRawVett_i(Sout);
	fprintf(fp_det, "\tSTout = ");
	_StampaRawVett_d(STout);
	fprintf(fp_det, "\tFreq_out = ");
	_StampaRawVett_d(Freq_out);
	fprintf(fp_det, "\ttoll = ");
	_StampaRawVett_d(toll);
#endif

	CREAv_d(g_scalare_d, 1);
	CREAv_i(g_scalare_i, 1);
//   nh<-length(h)
	nh = LENGTHv_i(h);
//   ng<-dim(M)[1]
	ng = LENGTHm1_i(M);
//   g_Sc<-g_checkOUT<-matrix(0,ng,nh)
	CREAm_i(g_checkOUT, ng, nh);
	InitMatr_i(g_checkOUT, 0);
	CREAm_d(ris, ng, nh);
	InitMatr_d(ris, 0.0);
//   M.out<-apply(M,2,sum)
	g_M_out = somma_colonne_i(g_M_out, M);
//   g_memory<-matrix(NA,ng,3)
	CREAm_i(g_memory, ng, 3);
	InitMatr_i(g_memory, NA_INTEGER);
//   Pm<-1
	Pm = 1.0;
//   lab<-"NO"
	strncpy(lab, "NO", 3);
//   i<-0
	i = 0;
//   while (i<ng)
	while (i < ng) {
//    {i<-i+1
		i++;
//     n.out<-M.out[i]
		n_out = ACCEDIv_i(g_M_out, i);
//	g_ind1<-which(g_memory[,2]==n.out)
		g_ind1 = which_m_colindxeq_i(g_ind1, g_memory, 2, n_out);
//     if (length(g_ind1)>0) {
		if (LENGTHv_i(g_ind1) > 0) {
			//ind<-g_memory[g_ind1,3]
			ind = copia_m_colindx_i(g_memory, g_ind1, 3);
			// g_Sc[i,]<-g_Sc[ind]
			ris = aggiungi_ms_riga_d(ris, i, ACCEDImv_d(ris, ind));
			// g_checkOUT[i,]<-g_checkOUT[ind,]
			g_checkOUT = aggiungi_riga_i(g_checkOUT, i, g_checkOUT, ind);
		}
	//      else
		else {
//			if (n.out!=0) {
			if (n_out != 0) {
		//				S.out<-Score(S=Sout[h],ST=STout,Freq=Freq.out,n=n.out,toll=toll)
				g_tmp_i1 = assegna_v_indxNA_i(g_tmp_i1, Sout, h);
				g_S_out = score1(g_S_out, g_tmp_i1, STout, Freq_out, n_out, toll);
//				indok<-which(S.out!=-Inf)
				g_indok = which_v_indxne_d(g_indok, g_S_out, R_NegInf);
//				checkOUT[i,indok]<-1
				g_checkOUT = aggiungi_ms_rigaindx_i(g_checkOUT, i, g_indok, 1);
	//           if (length(g_indok)<1) {Pm<-NA; i<-ng}
				if (LENGTHv_i(g_indok) < 1) {
					Pm = NA_REAL;
					i = ng;
				}
//				indInf <- setdiff(seq(1, nh, 1), indok)
				g_tmp_i1 = seq_i(g_tmp_i1, 1, nh, 1);
				g_indInf = setdiff_i(g_indInf, g_tmp_i1, g_indok);
//				S.out[indInf]<-min(c(0,S.out[indok]))-1
				g_tmp_d1 = copia_v_indx_d(g_tmp_d1, g_S_out, g_indok);
				if (LENGTHv_d(g_tmp_d1) == 0)
					mn =-1.0;
				else
					mn = min_s_d(0.0, min_v_d(g_tmp_d1)) - 1;
				assegna1_vs_indx_d(g_S_out, g_indInf, mn);
	//  	}
			}
	//        else {
			else {
				// S.out<-rep(0,nh)
				CREAv_d(g_S_out, nh);
				InitVett_d(g_S_out, 0.0);
				// g_checkOUT[i,]<-rep(1,nh)
				CREAv_i(g_tmp_i1, nh);
				InitVett_i(g_tmp_i1, 1);
				g_checkOUT = aggiungi_mv_riga_i(g_checkOUT, i, g_tmp_i1);
			}
//			Sc[i,]<-S.out
			ris = aggiungi_mv_riga_d(ris, i, g_S_out);
	//       }
		}
	//     }#end while (i<ng)
	}//end while (i<ng);
	//   if (!is.na(Pm))
	if (!ISNA(Pm)) {
	//    {Pm<-sum(g_Sc)/ng
		somma = somma_m_d(ris);
		Pm = somma / ng;
	//     g_checkOUT[which(g_checkOUT<0)]<-0
		assegna1_m_indxlt_i(g_checkOUT, g_checkOUT, 0, 0);
	//     g_checkOUT[which(g_checkOUT>0)]<-1
		assegna1_m_indxgt_i(g_checkOUT, g_checkOUT, 0, 1);
	//     g_aus<-apply(g_checkOUT,2,sum)
		g_aus = somma_colonne_i(g_aus, g_checkOUT);
	//     L<-length(which(g_aus!=0))
		g_tmp_i1 = which_v_indxne_i(g_tmp_i1, g_aus, 0);
		L = LENGTHv_i(g_tmp_i1);
	//     if (L<ng) {Pm<-NA; i<-ng}
		if (L < ng) {
			Pm = NA_REAL;
			i = ng;
		}
	//      else
		else {
	//       {g_rs<-apply(g_checkOUT,1,sum)
			g_rs = somma_righe_i(g_rs, g_checkOUT);
	//        ord.ind<-order(g_rs)
			g_ord_ind = ordine_i(g_ord_ind, g_rs, false);
	//        g_I<-vector()
			CREAv_i(g_I, 0);
	//        j<-0
			j = 0;
	//        while (j<ng)
			while (j < ng) {
	//    	{j<-j+1
				j++;
// 	 ind<-ord.ind[j]
				ind = ACCEDIv_i(g_ord_ind, j);
// 	 g_I.add<-which(g_checkOUT[ind,]!=0)
				g_I_add = which_m_rowindxne_i(g_I_add, g_checkOUT, ind, 0);
//          g_I<-union(g_I,g_I.add)
				g_I = unione1_i(g_I, g_I_add);
//          if (length(g_I)<j) {Pm<-NA; j<-ng}
				if (LENGTHv_i(g_I) < j) {
					Pm = NA_REAL;
					j = ng;
				}
	// 	}
			}
	//       }
		}
	//    }# end if (!is.na(Pm))
	}// end if (!is_na(Pm))
	// return(list(Pm,M,g_Sc,lab))

	ris1->score = Pm;

	strncpy(ris1->label, lab, 3);
	// ora cancello cio' che non mi serve piu`

	StrBilanciam();

	_Intestazione("\n*** Esco da probmod_und ***\n\n");

#ifdef DET
	fprintf(fp_det, "probmod output:\n");
	fprintf(fp_det, "\tlist(Pm,M,g_Sc,lab) =  %.7g\n", ris1->score);
	_StampaRawMatr_i(M);
	_StampaRawMatr_d(ris);
	fprintf(fp_det, "%s\n", ris1->label);
#endif

	return ris;
	//  }
}
