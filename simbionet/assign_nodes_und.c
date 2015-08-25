#include "assign_nodes_und.h"

#define g_Sout_h globali.assign_nodes_und.Sout_h
#define g_tmp_i globali.assign_nodes_und.tmp_i
#define g_tmp_d globali.assign_nodes_und.tmp_d
#define g_or_h globali.assign_nodes_und.or_h
#define g_aus_h globali.assign_nodes_und.aus_h
#define g_M_out globali.assign_nodes_und.M_out
#define g_Ord globali.assign_nodes_und.Ord
#define g_p globali.assign_nodes_und.p
#define g_ind_h globali.assign_nodes_und.ind_h
#define g_ri globali.assign_nodes_und.ri
#define g_co globali.assign_nodes_und.co
#define g_ind globali.assign_nodes_und.ind
#define g_index globali.assign_nodes_und.index

// NB: e` GIUSTO che h e Mdiscr vengano modificati!
VETTOREi *assign_nodes2_und(VETTOREi *ris, const MATRICEi *M, MATRICEi *Mdiscr, VETTOREi *h, const VETTOREi *hubs, const MATRICEd *Sc, const VETTOREi *Sout, int max_con)
{
	int i, j, n, ng, nh;
	double mn, somma;

	_Intestazione("\n***assign_nodes_und***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tM = ");
	_StampaRawMatr_i(M);
	fprintf(fp_det, "\tMdiscr = ");
	_StampaRawMatr_i(Mdiscr);
	fprintf(fp_det, "\th = ");
	_StampaRawVett_i(h);
	fprintf(fp_det, "\thubs = ");
	_StampaRawVett_i(hubs);
	fprintf(fp_det, "\tSc = ");
	_StampaRawMatr_d(Sc);
	fprintf(fp_det, "\tmax_con =  %d\n", max_con);
#endif

	 //~ nh<-length(h)
	nh = LENGTHv_i(h);
	 //~ ng<-dim(M)[1]
	ng = LENGTHm1_i(M);
	 //~ or.h<-h
	g_or_h = copia_v_i(g_or_h, h, 1, nh); // copia, come nell'instruzione di R
	 //~ aus.h<-seq(1,nh,1)
	g_aus_h = seq_i(g_aus_h, 1, nh, 1);
	 //~ new.hubs<-vector() // in realta` ha sempre almeno un elemento
	CREAv_i(ris, 1);
	ris->dim = 0;
	//  g_index<-rep(0,ng)
	CREAv_i(g_index, ng);
	InitVett_i(g_index, 0);
	// M.out<-apply(M,2,sum)
	g_M_out = somma_righe_i(g_M_out, M);
	// ORD<-order(M.in,decreasing=TRUE)
	g_Ord = ordine_i(g_Ord, g_M_out, true);
	//  for (j in (1:ng))
	for (j = 1; j <= ng; j++) {
		i = ACCEDIv_i(g_Ord, j);
		// g_p<-Sc[i,]
		g_p = riga_d(g_p, Sc, i);
		// if (length(which(g_p<=0))>0)
		g_tmp_i = which_v_indxle_d(g_tmp_i, g_p, 0.0);
		if (LENGTHv_i(g_tmp_i) > 0) {
			// g_p<-g_p-min(g_p)+1/(nh^2)
			mn = min_v_d(g_p);
			somma1_vs_d(g_p, -mn + (double) 1.0 / (nh * nh));
		}
		// Sout.h<-Sout[h]
		g_Sout_h = assegna_v_indxNA_i(g_Sout_h, Sout, h);
		// n<-M.out[i]
		if (i > 0) {
			n = ACCEDIv_i(g_M_out, i);
			// g_ind<-which(Sout.h>(max.con-n))
			g_ind = which_v_indxgt_i(g_ind, g_Sout_h, max_con - n);
			if (LENGTHv_i(g_ind) > 0)
				// g_p[g_ind]<-0
				assegna1_vs_indx_d(g_p, g_ind, 0.0);
		}
		// g_p<-g_p->dati[aus.h]/sum(g_p->dati[aus.h])
		g_tmp_d = assegna_v_indxNA_d(g_tmp_d, g_p, g_aus_h);
		somma = somma_v_d(g_tmp_d, false);
		g_p = dividi_vs_d(g_p, g_tmp_d, somma);
		// g_ind.h<-sampleB(h,size=1,prob=g_p)
		g_ind_h = sampleB_p(g_ind_h, h, 1, 0, g_p);
		// g_index[i]<-g_ind.h
		if (i > 0)
			ASSEGNAv_i(g_index, i, ACCEDIv_i(g_ind_h, 1));
		// h <- setdiff(h, g_ind.h)
	   setdiff1_i(h, g_ind_h);
		// aus.h <-setdiff(aus.h, which(or.h == g_ind.h))
		g_tmp_i = which_indx_vv_eq_i(g_tmp_i, g_or_h, g_ind_h);
		setdiff1_i(g_aus_h, g_tmp_i);
	}
	CREAv_i(g_ri, 1);
	//  for (j in (1:ng))
	for (j = 1; j <= ng; j++) {
		// g_ri<-g_index[j]
		ASSEGNAv_i(g_ri, 1, ACCEDIv_i(g_index, j));
		// g_ind <- which(M[j, ] == 1)
		g_ind = which_m_rowindxeq_i(g_ind, M, j, 1);
		// g_co<-g_index[g_ind]
		g_co = copia_v_indx_i(g_co, g_index, g_ind);
		// Mdiscr[g_ri,g_co]<-1
		assegna1_m_vv_i(Mdiscr, g_ri, g_co, 1);
		// Mdiscr[g_co,g_ri]<-1
		assegna1_m_vv_i(Mdiscr, g_co, g_ri, 1);
		// if (j % in % hubs)
		if (esiste_v_i(j, hubs) > 0)
			// new.hubs < -c(new.hubs, g_ri)
			ris = accoda1_vv_i(ris, g_ri);
	}

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "assign_nodes_und output:\n");
	fprintf(fp_det, "\tlist(Mdiscr,h,new_hubs) = ");
	_StampaRawMatr_i(Mdiscr);
	_StampaRawVett_i(h);
	_StampaRawVett_i(ris);
#endif

	// return(list(Mdiscr,h,new.hubs));
	return ris;
}
