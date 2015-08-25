#include "probmod.h"

#define g_Sc globali.probmod.Sc
#define g_checkIN globali.probmod.checkIN
#define g_checkOUT globali.probmod.checkOUT
#define g_memory globali.probmod.memory
#define g_M_in globali.probmod.M_in
#define g_M_out globali.probmod.M_out
#define g_tmp_i1 globali.probmod.tmp_i1
#define g_tmp_i2 globali.probmod.tmp_i2
#define g_tmp_i3 globali.probmod.tmp_i3
#define g_indInf globali.probmod.indInf
#define g_I globali.probmod.I
#define g_ord_ind globali.probmod.ord_ind
#define g_rs globali.probmod.rs
#define g_ind1 globali.probmod.ind1
#define g_scalare_i globali.probmod.scalare_i
#define g_I_add globali.probmod.I_add
#define g_S_in globali.probmod.S_in
#define g_S_out globali.probmod.S_out
#define g_tmp_d1 globali.probmod.tmp_d1

MATRICEd *probmod2(MATRICEd *ris, const MATRICEi *M, const VETTOREi *h, const VETTOREi *Sin, const VETTOREi *Sout, const VETTOREd *STin, const VETTOREd *STout, const VETTOREd *Freq_in, const VETTOREd *Freq_out, const VETTOREd *toll, struct RisProbM *ris1)
{
	int nh, ng, i, j, L, ind2, ind;
	int n_out, n_in;
	char lab[3];
	double mn, Pm, somma;

	_Intestazione("\n***probmod***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tM = ");
	_StampaRawMatr_i(M);
	fprintf(fp_det, "\th = ");
	_StampaRawVett_i(h);
	fprintf(fp_det, "\tSin = ");
	_StampaRawVett_i(Sin);
	fprintf(fp_det, "\tSout = ");
	_StampaRawVett_i(Sout);
	fprintf(fp_det, "\tSTin = ");
	_StampaRawVett_d(STin);
	fprintf(fp_det, "\tSTout = ");
	_StampaRawVett_d(STout);
	fprintf(fp_det, "\tFreq_in = ");
	_StampaRawVett_d(Freq_in);
	fprintf(fp_det, "\tFreq_out = ");
	_StampaRawVett_d(Freq_out);
	fprintf(fp_det, "\ttoll = ");
	_StampaRawVett_d(toll);
#endif

	CREAv_i(g_scalare_i, 1);
//   nh<-length(h)
	nh = LENGTHv_i(h);
//   ng<-dim(M)[1]
	ng = LENGTHm1_i(M);
//   g_Sc<-g_checkIN<-g_checkOUT<-matrix(0,ng,nh)
	CREAm_i(g_checkOUT, ng, nh);
	InitMatr_i(g_checkOUT, 0);
	CREAm_i(g_checkIN, ng, nh);
	InitMatr_i(g_checkIN, 0);
	CREAm_d(ris, ng, nh);
	InitMatr_d(ris, 0.0);
//   M.in<-apply(M,1,sum)
	g_M_in = somma_righe_i(g_M_in, M);
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
//     n.in<-M.in[i]
		n_in = ACCEDIv_i(g_M_in, i);
//     g_ind1<-intersect(which(g_memory[,1]==n.in),which(g_memory[,2]==n.out))
		g_tmp_i1 = which_m_colindxeq_i(g_tmp_i1, g_memory, 1, n_in);
		g_tmp_i2 = which_m_colindxeq_i(g_tmp_i2, g_memory, 2, n_out);
		g_ind1 = interseca_i(g_ind1, g_tmp_i1, g_tmp_i2);
//     if (length(g_ind1)>0) {
		if (LENGTHv_i(g_ind1) > 0) {
			//ind<-g_memory[g_ind1,3]
			ind = copia_m_colindx_i(g_memory, g_ind1, 3);
			// g_Sc[i,]<-g_Sc[ind]
			ris = aggiungi_ms_riga_d(ris, i, ACCEDImv_d(ris, ind));
			// g_checkIN[i,]<-g_checkIN[ind,]
			g_checkIN = aggiungi_riga_i(g_checkIN, i, g_checkIN, ind);
			// g_checkOUT[i,]<-g_checkOUT[ind,]
			g_checkOUT = aggiungi_riga_i(g_checkOUT, i, g_checkOUT, ind);
		}
	//      else
		else {
	//       {if (n.out!=0)
			if (n_out != 0) {
	// 	{S.out<-Score(S=Sout[h],ST=STout,Freq=Freq.out,n=n.out,toll=toll)
				g_tmp_i1 = assegna_v_indxNA_i(g_tmp_i1, Sout, h);
				g_S_out = score1(g_S_out, g_tmp_i1, STout, Freq_out, n_out, toll);
	// 	 indok<-which(S.out!=-Inf)
				g_tmp_i3 = which_v_indxne_d(g_tmp_i3, g_S_out, R_NegInf);
	//           g_checkOUT[i,indok]<-1
				g_checkOUT = aggiungi_ms_rigaindx_i(g_checkOUT, i, g_tmp_i3, 1);
	//           if (length(indok)<1) {Pm<-NA; i<-ng}
				if (LENGTHv_i(g_tmp_i3) < 1) {
					Pm = NA_REAL;
					i = ng;
				}
	//           g_indInf <- setdiff(seq(1, nh, 1), indok)
				g_tmp_i1 = seq_i(g_tmp_i1, 1, nh, 1);
				g_indInf = setdiff_i(g_indInf, g_tmp_i1, g_tmp_i3);
	// 	  S.out[g_indInf]<-min(c(0,S.out[indok]))-1
				g_tmp_d1 = copia_v_indx_d(g_tmp_d1, g_S_out, g_tmp_i3);
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
	//        if (n.in!=0)
			if (n_in != 0) {
	// 	 {S.in<-Score(S=Sin[h],ST=STin,Freq=Freq.in,n=n.in,toll=toll)
				g_tmp_i1 = assegna_v_indxNA_i(g_tmp_i1, Sin, h);
				g_S_in = score1(g_S_in, g_tmp_i1, STin, Freq_in, n_in, toll);
	// 	  indok<-which(S.in!=-Inf)
				g_tmp_i3 = which_v_indxne_d(g_tmp_i3, g_S_in, R_NegInf);
	//           g_checkIN[i,indok]<-1
				g_checkIN = aggiungi_ms_rigaindx_i(g_checkIN, i, g_tmp_i3, 1);
	//           if (length(indok)<1) {Pm<-NA; i<-ng; lab<-"in"}
				if (LENGTHv_i(g_tmp_i3) < 1) {
					Pm = NA_REAL;
					i = ng;
					strncpy(lab, "in", 3);
				}
	//          g_indInf <- setdiff(seq(1, nh, 1), indok)
				g_tmp_i1 = seq_i(g_tmp_i1, 1, nh, 1);
				g_indInf = setdiff_i(g_indInf, g_tmp_i1, g_tmp_i3);
	// 	  S.in[g_indInf]<-min(c(0,S.in[indok]))-1
				g_tmp_d1 = copia_v_indx_d(g_tmp_d1, g_S_in, g_tmp_i3);
				if (LENGTHv_d(g_tmp_d1) == 0)
					mn =-1.0;
				else
					mn = min_s_d(0.0, min_v_d(g_tmp_d1)) - 1;
				assegna1_vs_indx_d(g_S_in, g_indInf, mn);
	//       	 }
			}
	//        else {S.in<-rep(0,nh); g_checkIN[i,]<-1}
			else {
				CREAv_d(g_S_in, nh);
				InitVett_d(g_S_in, 0.0);
				g_checkIN = aggiungi_ms_riga_i(g_checkIN, i, 1);
			}
	//        g_Sc[i,]<-S.out+S.in
			g_tmp_d1 = somma_vv_d(g_tmp_d1, g_S_out, g_S_in);
			ris = aggiungi_mv_riga_d(ris, i, g_tmp_d1);
	//       }
		}
	//     }#end while (i<ng)
	}//end while (i<ng);
	//   if (!is.na(Pm))
	if (!ISNA(Pm)) {
	//    {Pm<-sum(g_Sc)/ng
		somma = somma_m_d(ris);
		Pm = somma / ng;
	//     #check in
	 //check in;
	//     aus<-apply(g_checkIN,2,sum)
		g_tmp_i3 = somma_colonne_i(g_tmp_i3, g_checkIN);
	//     L<-length(which(aus!=0))
		g_tmp_i1 = which_v_indxne_i(g_tmp_i1, g_tmp_i3, 0);
		L = LENGTHv_i(g_tmp_i1);
	//     if (L<ng) {Pm<-NA; i<-ng; lab<-"in"}
		if (L < ng) {
			Pm = NA_REAL;
			i = ng;
			strncpy(lab,"in", 3);
		}
	//      else
		else {
	//       {g_rs<-apply(g_checkIN,1,sum)
			g_rs = somma_righe_i(g_rs, g_checkIN);
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
				ind2 = ACCEDIv_i(g_ord_ind, j);
	// 	 g_I.add<-which(g_checkIN[ind,]!=0)

				g_I_add = which_m_rowindxne_i(g_I_add, g_checkIN, ind2, 0);
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
	//     #check out
		 //check out;
	//     g_checkOUT[which(g_checkOUT<0)]<-0
		assegna1_m_indxlt_i(g_checkOUT, g_checkOUT, 0, 0);
	//     g_checkOUT[which(g_checkOUT>0)]<-1
		assegna1_m_indxgt_i(g_checkOUT, g_checkOUT, 0, 1);
	//     aus<-apply(g_checkOUT,2,sum)
		g_tmp_i3 = somma_colonne_i(g_tmp_i3, g_checkOUT);
	//     L<-length(which(aus!=0))
		g_tmp_i1 = which_v_indxne_i(g_tmp_i1, g_tmp_i3, 0);
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
	// ho assegnato prima il viceversa, quindi non dovrebbe servire
	//~ g_score_matr = g_Sc;
	strncpy(ris1->label, lab, 3);
	// ora cancello cio' che non mi serve piu`
	//~ CANCELLAv_d(g_S_out);
	//~ CANCELLAv_d(g_S_in);
	//~ //CANCELLAm_d(g_Sc);
	//~ CANCELLAm_i(g_checkIN);
	//~ CANCELLAm_i(g_checkOUT);
	//~ CANCELLAm_i(g_memory);
	//~ CANCELLAv_i(g_M_in);
	//~ CANCELLAv_i(g_M_out);
	//~ CANCELLAv_i(g_tmp_i1);
	//~ CANCELLAv_i(g_tmp_i2);
	//~ CANCELLAv_i(g_tmp_i3);
	//~ CANCELLAv_i(g_indInf);
	//~ CANCELLAv_i(g_I);
	//~ CANCELLAv_i(g_ord_ind);
	//~ CANCELLAv_i(g_rs);
	//~ CANCELLAv_i(g_ind1);
	//~ CANCELLAv_i(g_scalare_i);
	//~ CANCELLAv_i(g_I_add);
	//~ CANCELLAv_d(g_S_in);
	//~ CANCELLAv_d(g_S_out);
	//~ CANCELLAv_d(g_tmp_d1);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "probmod output:\n");
	fprintf(fp_det, "\tlist(Pm,M,Sc,lab) =  %.16g\n", ris1->score);
	_StampaRawMatr_i(M);
	_StampaRawMatr_d(ris);
	fprintf(fp_det, "%s\n", ris1->label);
#endif

	return ris;
	//  }
}

// probmod<-function(M,h,Sin,Sout,STin,STout,Freq.in,Freq.out,toll)
LISTA *probmod1(LISTA *ris, MATRICEi *M, VETTOREi *h, VETTOREi *Sin, VETTOREi *Sout, VETTOREd *STin, VETTOREd *STout, VETTOREd *Freq_in, VETTOREd *Freq_out, VETTOREd *toll)
{
	int nh, ng, i, j, L, ind2, ind;
	int n_out, n_in;
	GString *lab = NULL;
	double mn, Pm, somma;
	MATRICEi *checkOUT = NULL, *checkIN = NULL, *memory = NULL;
	VETTOREi *M_in = NULL, *M_out = NULL;
	VETTOREi *tmp_i1 = NULL, *tmp_i2 = NULL, *indok = NULL, *indInf = NULL,
		*I = NULL, *ord_ind = NULL, *aus = NULL, *rs = NULL, *ind1 = NULL, *scalare_i = NULL, *I_add = NULL;
	VETTOREd *S_in = NULL, *S_out = NULL, *tmp_d1 = NULL;
	MATRICEd *Sc = NULL;
	enum TIPO tipi[5];

	_Intestazione("\n*** probmod1 ***\n");

	CREAv_i(scalare_i, 1);
//   nh<-length(h)
	nh = LENGTHv_i(h);
//   ng<-dim(M)[1]
	ng = LENGTHm1_i(M);
//   Sc<-checkIN<-checkOUT<-matrix(0,ng,nh)
	CREAm_i(checkOUT, ng, nh);
	InitMatr_i(checkOUT, 0);
	CREAm_i(checkIN, ng, nh);
	InitMatr_i(checkIN, 0);
	// Sc potrebbe gia` essere stato allocato!
	if (ris != NULL) {
		CtrlLlst(ris, 3);
		Sc = ACCEDIlst(ris, 3, md);
	}
	CREAm_d(Sc, ng, nh);
	InitMatr_d(Sc, 0.0);
//   M.in<-apply(M,1,sum)
	M_in = somma_righe_i(M_in, M);
//   M.out<-apply(M,2,sum)
	M_out = somma_colonne_i(M_out, M);
//   memory<-matrix(NA,ng,3)
	CREAm_i(memory, ng, 3);
	InitMatr_i(memory, NA_INTEGER);
//   Pm<-1
	Pm = 1.0;
//   lab<-"NO"
	CREAstr(lab, "NO");
//   i<-0
	i = 0;
//   while (i<ng)
	while (i < ng) {
//    {i<-i+1
		i++;
//     n.out<-M.out[i]
		n_out = ACCEDIv_i(M_out, i);
//     n.in<-M.in[i]
		n_in = ACCEDIv_i(M_in, i);
//     ind1<-intersect(which(memory[,1]==n.in),which(memory[,2]==n.out))
		tmp_i1 = which_m_colindxeq_i(tmp_i1, memory, 1, n_in);
		tmp_i2 = which_m_colindxeq_i(tmp_i2, memory, 2, n_out);
		ind1 = interseca_i(ind1, tmp_i1, tmp_i2);
//     if (length(ind1)>0) {
		if (LENGTHv_i(ind1) > 0) {
			//ind<-memory[ind1,3]
			ind = copia_m_colindx_i(memory, ind1, 3);
			// Sc[i,]<-Sc[ind]
			Sc = aggiungi_ms_riga_d(Sc, i, ACCEDImv_d(Sc, ind));
			// checkIN[i,]<-checkIN[ind,]
			checkIN = aggiungi_riga_i(checkIN, i, checkIN, ind);
			// checkOUT[i,]<-checkOUT[ind,]
			checkOUT = aggiungi_riga_i(checkOUT, i, checkOUT, ind);
		}
	//      else
		else {
	//       {if (n.out!=0)
			if (n_out != 0) {
	// 	{S.out<-Score(S=Sout[h],ST=STout,Freq=Freq.out,n=n.out,toll=toll)
				tmp_i1 = copia_v_indx_i(tmp_i1, Sout, h);
				S_out = score1(S_out, tmp_i1, STout, Freq_out, n_out, toll);
	// 	 indok<-which(S.out!=-Inf)
				indok = which_v_indxne_d(indok, S_out, R_NegInf);
	//           checkOUT[i,indok]<-1
				checkOUT = aggiungi_ms_rigaindx_i(checkOUT, i, indok, 1);
	//           if (length(indok)<1) {Pm<-NA; i<-ng}
				if (LENGTHv_i(indok) < 1) {
					Pm = NA_REAL;
					i = ng;
				}
	//           indInf <- setdiff(seq(1, nh, 1), indok)
				tmp_i1 = seq_i(tmp_i1, 1, nh, 1);
				indInf = setdiff_i(indInf, tmp_i1, indok);
	// 	  S.out[indInf]<-min(c(0,S.out[indok]))-1
				tmp_d1 = copia_v_indx_d(tmp_d1, S_out, indok);
				if (LENGTHv_d(tmp_d1) == 0)
					mn =-1.0;
				else
					mn = min_s_d(0.0, min_v_d(tmp_d1)) - 1;
				assegna1_vs_indx_d(S_out, indInf, mn);
	//  	}
			}
	//        else {
			else {
				// S.out<-rep(0,nh)
				CREAv_d(S_out, nh);
				InitVett_d(S_out, 0.0);
				// checkOUT[i,]<-rep(1,nh)
				CREAv_i(tmp_i1, nh);
				InitVett_i(tmp_i1, 1);
				checkOUT = aggiungi_mv_riga_i(checkOUT, i, tmp_i1);
			}
	//        if (n.in!=0)
			if (n_in != 0) {
	// 	 {S.in<-Score(S=Sin[h],ST=STin,Freq=Freq.in,n=n.in,toll=toll)
				tmp_i1 = assegna_v_indxNA_i(tmp_i1, Sin, h);
				S_in = score1(S_in, tmp_i1, STin, Freq_in, n_in, toll);
	// 	  indok<-which(S.in!=-Inf)
				indok = which_v_indxne_d(indok, S_in, R_NegInf);
	//           checkIN[i,indok]<-1
				checkIN = aggiungi_ms_rigaindx_i(checkIN, i, indok, 1);
	//           if (length(indok)<1) {Pm<-NA; i<-ng; lab<-"in"}
				if (LENGTHv_i(indok) < 1) {
					Pm = NA_REAL;
					i = ng;
					lab = g_string_assign(lab, "in");
				}
	//          indInf <- setdiff(seq(1, nh, 1), indok)
				tmp_i1 = seq_i(tmp_i1, 1, nh, 1);
				indInf = setdiff_i(indInf, tmp_i1, indok);
	// 	  S.in[indInf]<-min(c(0,S.in[indok]))-1
				tmp_d1 = copia_v_indx_d(tmp_d1, S_in, indok);
				if (LENGTHv_d(tmp_d1) == 0)
					mn =-1.0;
				else
					mn = min_s_d(0.0, min_v_d(tmp_d1)) - 1;
				assegna1_vs_indx_d(S_in, indInf, mn);
	//       	 }
			}
	//        else {S.in<-rep(0,nh); checkIN[i,]<-1}
			else {
				CREAv_d(S_in, nh);
				InitVett_d(S_in, 0.0);
				checkIN = aggiungi_ms_riga_i(checkIN, i, 1);
			}
	//        Sc[i,]<-S.out+S.in
			tmp_d1 = somma_vv_d(tmp_d1, S_out, S_in);
			Sc = aggiungi_mv_riga_d(Sc, i, tmp_d1);
	//       }
		}
	//     }#end while (i<ng)
	}//end while (i<ng);
	//   if (!is.na(Pm))
	if (Pm != NA_REAL) {
	//    {Pm<-sum(Sc)/ng
		somma = somma_m_d(Sc);
		Pm = somma / ng;
	//     #check in
	 //check in;
	//     aus<-apply(checkIN,2,sum)
		aus = somma_colonne_i(aus, checkIN);
	//     L<-length(which(aus!=0))
		tmp_i1 = which_v_indxne_i(tmp_i1, aus, 0);
		L = LENGTHv_i(tmp_i1);
	//     if (L<ng) {Pm<-NA; i<-ng; lab<-"in"}
		if (L < ng) {
			Pm = NA_REAL;
			i = ng;
			lab = g_string_assign(lab, "in");
		}
	//      else
		else {
	//       {rs<-apply(checkIN,1,sum)
			rs = somma_righe_i(rs, checkIN);
	//        ord.ind<-order(rs)
			ord_ind = ordine_i(ord_ind, rs, false);
	//        I<-vector()
			CREAv_i(I, 0);
	//        j<-0
			j = 0;
	//        while (j<ng)
			while (j < ng) {
	//    	{j<-j+1
				j++;
	// 	 ind<-ord.ind[j]
				ind2 = ACCEDIv_i(ord_ind, j);
	// 	 I.add<-which(checkIN[ind,]!=0)

				I_add = which_m_rowindxne_i(I_add, checkIN, ind2, 0);
	//          I<-union(I,I.add)
				I = unione1_i(I, I_add);
	//          if (length(I)<j) {Pm<-NA; j<-ng}
				if (LENGTHv_i(I) < j) {
					Pm = NA_REAL;
					j = ng;
				}
	// 	}
			}
	//       }
		}
	//     #check out
		 //check out;
	//     checkOUT[which(checkOUT<0)]<-0
		assegna1_m_indxlt_i(g_checkOUT, g_checkOUT, 0, 0);
	//     checkOUT[which(checkOUT>0)]<-1
		assegna1_m_indxgt_i(g_checkOUT, g_checkOUT, 0, 1);
	//     aus<-apply(checkOUT,2,sum)
		aus = somma_colonne_i(aus, checkOUT);
	//     L<-length(which(aus!=0))
		tmp_i1 = which_v_indxne_i(tmp_i1, aus, 0);
		L = LENGTHv_i(tmp_i1);
	//     if (L<ng) {Pm<-NA; i<-ng}
		if (L < ng) {
			Pm = NA_REAL;
			i = ng;
		}
	//      else
		else {
	//       {rs<-apply(checkOUT,1,sum)
			rs = somma_righe_i(rs, checkOUT);
	//        ord.ind<-order(rs)
			ord_ind = ordine_i(ord_ind, rs, false);
	//        I<-vector()
			CREAv_i(I, 0);
	//        j<-0
			j = 0;
	//        while (j<ng)
			while (j < ng) {
	//    	{j<-j+1
				j++;
// 	 ind<-ord.ind[j]
				ind = ACCEDIv_i(ord_ind, j);
// 	 I.add<-which(checkOUT[ind,]!=0)
				I_add = which_m_rowindxne_i(I_add, checkOUT, ind, 0);
//          I<-union(I,I.add)
				I = unione1_i(I, I_add);
//          if (length(I)<j) {Pm<-NA; j<-ng}
				if (LENGTHv_i(I) < j) {
					Pm = NA_REAL;
					j = ng;
				}
	// 	}
			}
	//       }
		}
	//    }# end if (!is.na(Pm))
	}// end if (!is_na(Pm))
	// return(list(Pm,M,Sc,lab))

	tipi[0] = VETTd;
	tipi[1] = MATRi;
	tipi[2] = MATRd;
	tipi[3] = STRINGA;
	// creo già qui il 5° elemento
	tipi[4] = VETTi;
	// in realta` CreaLISTA potrebbe anche fare a meno di questo controllo: e` solo per inizializzare alcuni suoi elementi
	if (ris == NULL) {
		CreaLISTA(ris, tipi, 5);
		CtrlLlst(ris, 1);
		CREAv_d(ACCEDIlst(ris, 1, vd), 1);
		CtrlLlst(ris, 3);
		CREAstr(ACCEDIlst(ris, 3, str), "");
		// non occorre creare il 5' elemento, perche´ verra` semplicemente assegnato (quindi usato come un puntatore)
	}
	CtrlLlst(ris, 1);
	ASSEGNAv_d(ACCEDIlst(ris, 1, vd), 1, Pm);
	CtrlSlst(ris, 2);
	ASSEGNAlst(ris, 2, mi, M);
	CtrlSlst(ris, 3);
	ASSEGNAlst(ris, 3, md, Sc);
	CtrlLlst(ris, 4);
	g_string_assign(ACCEDIlst(ris, 4, str), lab->str);
	// ora cancello cio' che non mi serve piu'
	CANCELLAv_d(S_out);
	CANCELLAv_d(S_in);
	// Sc viene passata in uscita, NON devo cancellarla!
	//CANCELLAm_d(Sc);
	CANCELLAm_i(checkIN);
	CANCELLAm_i(checkOUT);
	CANCELLAm_i(memory);
	CANCELLAv_i(M_in);
	CANCELLAv_i(M_out);
	CANCELLAv_i(tmp_i1);
	CANCELLAv_i(tmp_i2);
	CANCELLAv_i(indok);
	CANCELLAv_i(indInf);
	CANCELLAv_i(I);
	CANCELLAv_i(ord_ind);
	CANCELLAv_i(aus);
	CANCELLAv_i(rs);
	CANCELLAv_i(ind1);
	CANCELLAv_i(scalare_i);
	CANCELLAv_i(I_add);
	CANCELLAv_d(S_in);
	CANCELLAv_d(S_out);
	CANCELLAv_d(tmp_d1);

	StrBilanciam();

	return ris;
	//  }
}

SEXP probmod(SEXP m, SEXP h, SEXP Sin, SEXP Sout, SEXP STin, SEXP STout, SEXP Freq_in, SEXP Freq_out, SEXP toll)
{
	int nProtected = 0;
	MATRICEi *m1;
	VETTOREi *h1, *Sin1, *Sout1;
	VETTOREd *STin1, *STout1, *Freq_in1, *Freq_out1, *toll1;
	LISTA *l = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** probmod ***\n");

	m1 = inMATRICE_i(m, &nProtected);
	h1 = inVETTORE_i(h, &nProtected);
	Sin1 = inVETTORE_i(Sin, &nProtected);
	Sout1 = inVETTORE_i(Sout, &nProtected);
	Freq_in1 = inVETTORE_d(Freq_in, &nProtected);
	Freq_out1 = inVETTORE_d(Freq_out, &nProtected);
	toll1 = inVETTORE_d(toll, &nProtected);
	STin1 = inVETTORE_d(STin, &nProtected);
	STout1 = inVETTORE_d(STout, &nProtected);

	l = probmod1(l, m1, h1 ,Sin1 , Sout1, STin1, STout1, Freq_in1, Freq_out1, toll1);
	ris = daLISTA(l, &nProtected);

	CANCELLAv_i(h1);
	CANCELLAv_i(Sin1);
	CANCELLAm_i(m1);
	CANCELLAv_i(Sout1);
	CANCELLAv_d(STin1);
	CANCELLAv_d(Freq_in1);
	CANCELLAv_d(Freq_out1);
	CANCELLAv_d(STout1);
	CANCELLAv_d(toll1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
