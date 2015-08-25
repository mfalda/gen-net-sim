#include "assign_nodes.h"

#define g_Sin_h globali.assign_nodes.Sin_h
#define g_tmp_i globali.assign_nodes.tmp_i
#define g_tmp_d globali.assign_nodes.tmp_d
#define g_or_h globali.assign_nodes.or_h
#define g_aus_h globali.assign_nodes.aus_h
#define g_M_in globali.assign_nodes.M_in
#define g_Ord globali.assign_nodes.Ord
#define g_p globali.assign_nodes.p
#define g_ind_h globali.assign_nodes.ind_h
#define g_ri globali.assign_nodes.ri
#define g_co globali.assign_nodes.co
#define g_ind globali.assign_nodes.ind
#define g_index globali.assign_nodes.index

// Mdiscr e h verranno modificate, per cui non le restituisco!
VETTOREi *assign_nodes2(VETTOREi *ris, const MATRICEi *M, MATRICEi *Mdiscr, VETTOREi *h, const VETTOREi *hubs, const MATRICEd *Sc, const VETTOREi *Sin, int max_con)
{
	int i, j, n, ng, nh;
	double mn, somma;

	_Intestazione("\n***assign_nodes***\n");
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
	fprintf(fp_det, "\tSin = ");
	_StampaRawVett_i(Sin);
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
	// M.in<-apply(M,1,sum)
	g_M_in = somma_righe_i(g_M_in, M);
	// ORD<-order(M.in,decreasing=TRUE)
	g_Ord = ordine_i(g_Ord, g_M_in, true);
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
		// Sin.h<-Sin[h]
		g_Sin_h = assegna_v_indxNA_i(g_Sin_h, Sin, h);
		// n<-M.in[i]
		if (i > 0) {
			n = ACCEDIv_i(g_M_in, i);
			// g_ind<-which(Sin.h>(max.con-n))
			g_ind = which_v_indxgt_i(g_ind, g_Sin_h, max_con - n);
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
		// if (j % in % hubs)
		if (esiste_v_i(j, hubs) > 0)
			// new.hubs < -c(new.hubs, g_ri)
			ris = accoda1_vv_i(ris, g_ri);
	}
	//~ CANCELLAv_i(g_Sin_h);
	//~ CANCELLAv_i(g_tmp_i);
	//~ CANCELLAv_d(g_tmp_d);
	//~ CANCELLAv_i(g_or_h);
	//~ CANCELLAv_i(g_aus_h);
	//~ CANCELLAv_i(g_M_in);
	//~ CANCELLAv_i(g_Ord);
	//~ CANCELLAv_d(g_p);
	//~ CANCELLAv_i(g_ind_h);
	//~ CANCELLAv_i(g_ri);
	//~ CANCELLAv_i(g_co);
	//~ CANCELLAv_i(g_ind);
   //~ CANCELLAv_i(g_index);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "assign_nodes output:\n");
	fprintf(fp_det, "\tlist(Mdiscr,h,new_hubs) = ");
	_StampaRawMatr_i(Mdiscr);
	_StampaRawVett_i(h);
	_StampaRawVett_i(ris);
#endif
	// return(list(Mdiscr,h,new.hubs));
	return ris;
}

LISTA *assign_nodes1(LISTA *ris, MATRICEi *M, MATRICEi *Mdiscr, VETTOREi *h, VETTOREi *hubs, MATRICEd *Sc, VETTOREi *Sin, int max_con)
{
	int i, j, n, ng, nh;
	double mn, somma;
	VETTOREi *or_h = NULL, *new_hubs= NULL, *Sin_h= NULL, *index= NULL, *Ord= NULL, *ind= NULL, *tmp_i= NULL, *tmp1_i= NULL, *aus_h= NULL, *ind_h= NULL, *M_in= NULL, *co= NULL, *ri = NULL;
	VETTOREd *p = NULL, *tmp_d = NULL;
	enum TIPO tipi[3];

	_Intestazione("\n*** assign_nodes1 ***\n");

	 //~ nh<-length(h)
	nh = LENGTHv_i(h);
	 //~ ng<-dim(M)[1]
	ng = LENGTHm1_i(M);
	 //~ or.h<-h
	or_h = copia_v_i(or_h, h, 1, nh); // copia, come nell'instruzione di R
	 //~ aus.h<-seq(1,nh,1)
	aus_h = seq_i(aus_h, 1, nh, 1);
	 //~ new.hubs<-vector()
	CREAv_i(new_hubs, 0);
	//  index<-rep(0,ng)
	CREAv_i(index, ng);
	InitVett_i(index, 0);
	// M.in<-apply(M,1,sum)
	M_in = somma_righe_i(M_in, M);
	// ORD<-order(M.in,decreasing=TRUE)
	Ord = ordine_i(Ord, M_in, true);
	//  for (j in (1:ng))
	for (j = 1; j <= ng; j++) {
		i = ACCEDIv_i(Ord, j);
		// p<-Sc[i,]
		p = riga_d(p, Sc, i);
		// if (length(which(p<=0))>0)
		tmp_i = which_v_indxle_d(tmp_i, p, 0.0);
		if (LENGTHv_i(tmp_i) > 0) {
			// p<-p-min(p)+1/(nh^2)
			mn = min_v_d(p);
			somma1_vs_d(p, -mn + nh * nh);
		}
		// Sin.h<-Sin[h]
		Sin_h = assegna_v_indxNA_i(Sin_h, Sin, h);
		// n<-M.in[i]
		if (i > 0) {
			n = ACCEDIv_i(M_in, i);
			// ind<-which(Sin.h>(max.con-n))
			ind = which_v_indxgt_i(ind, Sin_h, max_con - n);
			if (LENGTHv_i(ind) > 0)
				// p[ind]<-0
				assegna1_vs_indx_d(p, ind, 0.0);
		}
		// p<-p->dati[aus.h]/sum(p->dati[aus.h])
		tmp_d = assegna_v_indxNA_d(tmp_d, p, aus_h);
		somma = somma_v_d(tmp_d, false);
		p = dividi_vs_d(p, tmp_d, somma);
		// ind.h<-sampleB(h,size=1,prob=p)
		ind_h = sampleB_p(ind_h, h, 1, 0, p);
		// index[i]<-ind.h
		if (i > 0)
			ASSEGNAv_i(index, i, ACCEDIv_i(ind_h, 1));
		// h <- setdiff(h, ind.h)
	   setdiff1_i(h, ind_h);
		// aus.h <-setdiff(aus.h, which(or.h == ind.h))
		tmp1_i = which_indx_vv_eq_i(tmp_i, or_h, ind_h);
		setdiff1_i(aus_h, tmp1_i);
	}
	CANCELLAv_i(Sin_h);
	CANCELLAv_i(tmp_i);
	CANCELLAv_d(tmp_d);
	CANCELLAv_i(or_h);
	CANCELLAv_i(aus_h);
	CANCELLAv_i(M_in);
	CANCELLAv_i(Ord);
	if (ng < 1)
		CANCELLAv_i(tmp1_i);
	CANCELLAv_d(p);
	CANCELLAv_i(ind_h);

	CREAv_i(ri, 1);
	//  for (j in (1:ng))
	for (j = 1; j <= ng; j++) {
		// ri<-index[j]
		ASSEGNAv_i(ri, 1, ACCEDIv_i(index, j));
		// ind <- which(M[j, ] == 1)
		ind = which_m_rowindxeq_i(ind, M, j, 1);
		// co<-index[ind]
		co = copia_v_indx_i(co, index, ind);
		// Mdiscr[ri,co]<-1
		assegna1_m_vv_i(Mdiscr, ri, co, 1);
		// if (j % in % hubs)
		if (esiste_v_i(j, hubs) > 0)
			// new.hubs < -c(new.hubs, ri)
			new_hubs = accoda1_vv_i(new_hubs, ri);
	}
	CANCELLAv_i(ri);
	CANCELLAv_i(co);
	CANCELLAv_i(ind);
	// return(list(Mdiscr,h,new.hubs));
	tipi[0] = MATRi;
	tipi[1] = VETTi;
	tipi[2] = VETTi;
	CreaLISTA(ris, tipi, 3);
	CtrlSlst(ris, 1);
	ASSEGNAlst(ris, 1, mi, Mdiscr);
	CtrlSlst(ris, 2);
	ASSEGNAlst(ris, 2, vi, h);
	CtrlSlst(ris, 3);
	ASSEGNAlst(ris, 3, vi, new_hubs);

	StrBilanciam();

	return ris;
}

SEXP assign_nodes(SEXP M, SEXP Mdiscr, SEXP h, SEXP hubs, SEXP Sc, SEXP Sin, SEXP max_con)
{
	int nProtected = 0;
	MATRICEi *M1, *Mdiscr1;
	MATRICEd *Sc1;
	VETTOREi *h1, *hubs1, *Sin1;
	int max_con1;
	LISTA *l = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** assign_nodes ***\n");

	M1 = inMATRICE_i(M, &nProtected);
	Mdiscr1 = inMATRICE_i(Mdiscr, &nProtected);
	h1 = inVETTORE_i(h, &nProtected);
	hubs1 = inVETTORE_i(hubs, &nProtected);
	Sc1 = inMATRICE_d(Sc, &nProtected);
	Sin1 = inVETTORE_i(Sin, &nProtected);
	max_con1 = INTEGER_VALUE(max_con);

	l = assign_nodes1(l, M1, Mdiscr1, h1, hubs1, Sc1, Sin1, max_con1);
	ris = daLISTA(l, &nProtected);

	CANCELLAm_i(M1);
	CANCELLAv_i(hubs1);
	CANCELLAm_d(Sc1);
	CANCELLAv_i(Sin1);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
