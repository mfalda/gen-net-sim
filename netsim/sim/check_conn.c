#include "check_conn.h"

#define g_tmp_i globali.check_conn.tmp_i
#define g_colore globali.check_conn.colore
#define g_grigi globali.check_conn.grigi
#define g_ind globali.check_conn.ind
#define g_adj globali.check_conn.adj
#define g_tmp1_d globali.check_conn.tmp1_d
#define g_scalare_d globali.check_conn.scalare_d

VETTOREd *check_conn1(VETTOREd *ris, MATRICEi *Mdiscr)
{
	int N, vert, u, i, v, L, lg;

	_Intestazione("\n*** check_conn1 ***\n");

	CREAv_d(g_scalare_d, 1);
	CREAv_i(g_grigi, LENGTHm1_i(Mdiscr));
	N = LENGTHm1_i(Mdiscr);
	g_Mt = trasponi_i(g_Mt, Mdiscr);
	g_Maus = somma_mm_i(g_Maus, Mdiscr, Mt);
	// ind<-which(Maus!=0,arr.ind=TRUE)
	g_ind = which_m_indxne_i(g_ind, g_Maus, 0);

	// Maus[ind]<-1
	assegna1_ms_indx_i(g_Maus, g_ind, 1);
	vert = 1;
	// colore <- c(1, rep(0, N - 1))
	CREAv_i(g_tmp_i, N - 1);
	InitVett_i(g_tmp_i, 0);
	ASSEGNAv_i(g_scalare_i, 1, 1);
	g_colore = vettore2v_i(g_colore, scalare_i, tmp_i);
	// dist <- c(0, rep(Inf, N - 1))
	CREAv_d(g_tmp1_d, N - 1);
	InitVett_d(g_tmp1_d, R_PosInf);
	ASSEGNAv_d(g_scalare_d, 1, 0.0);
	ris = vettore2v_d(ris, g_scalare_d, g_tmp1_d);
	ASSEGNAv_i(g_grigi, 1, 1);
	while (LENGTHv_i(g_grigi) != 0) {
		u = ACCEDIv_i(g_grigi, 1);
		//~ adj<-which(Maus[u,]!=0)
		g_adj = which_m_rowindxne_i(g_adj, g_Maus, u, 0);
		L = LENGTHv_i(g_adj);
		if (L > 0) {
			for (i = 0; i < L; i++) {
				v = ACCEDIv_i(g_adj, i);
				if (ACCEDIv_i(g_colore, v) == 0)	{
					ASSEGNAv_i(g_colore, v, 1);
					ASSEGNAv_d(ris, v, ACCEDIv_d(ris, u) + 1);
					// grigi<-c(grigi,v)
					g_grigi = accoda1_vs_i(g_grigi, v);
				}
			}
		}
		lg = LENGTHv_i(g_grigi);
		if (lg >= 2)
			// grigi<-grigi[2:lg]
			segmento1_v_i(g_grigi, 2, lg);
		else
			InitVett_i(g_grigi, 0);
	}
	//~ CANCELLAv_i(tmp_i);
	//~ CANCELLAv_i(colore);
	//~ CANCELLAv_i(grigi);
	//~ CANCELLAv_i(ind);
	//~ CANCELLAv_i(adj);
	//~ CANCELLAv_d(tmp1_d);
	//~ CANCELLAv_d(scalare_d);

	StrBilanciam();

	return ris;
}

SEXP check_conn(SEXP Mdiscr)
{
	int nProtected = 0;
	MATRICEi *Mdiscr1;
	VETTOREd *v1;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** check_conn ***\n");

	Mdiscr1 = inMATRICE_i(Mdiscr, &nProtected);

	v1 = check_conn1(Mdiscr1);

	ris = daVETTORE_d(v1, &nProtected);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
