#include "check_conn.h"

#define g_Mt globali.check_conn.Mt
#define g_Maus globali.check_conn.Maus
#define g_tmp_i globali.check_conn.tmp_i
#define g_colore globali.check_conn.colore
#define g_grigi globali.check_conn.grigi
#define g_ind globali.check_conn.ind
#define g_adj globali.check_conn.adj
#define g_tmp_i globali.check_conn.tmp_i
#define g_tmp1_d globali.check_conn.tmp1_d
#define g_scalare_i globali.check_conn.scalare_i
#define g_scalare_d globali.check_conn.scalare_d

VETTOREd *check_conn1(VETTOREd *ris, const MATRICEi *Mdiscr)
{
	int N, vert, u, i, v, L, lg;

	_Intestazione("\n***check_conn***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tMdiscr = ");
	_StampaRawMatr_i(Mdiscr);
#endif

	CREAv_i(g_scalare_i, 1);
	CREAv_d(g_scalare_d, 1);
	// g_grigi<-1
	CREAv_i(g_grigi, 1);
	ASSEGNAv_i(g_grigi, 1, 1);
	N = LENGTHm1_i(Mdiscr);
	g_Mt = trasponi_i(g_Mt, Mdiscr);
	g_Maus = somma_mm_i(g_Maus, Mdiscr, g_Mt);
	// g_ind<-which(g_Maus!=0,arr.g_ind=TRUE)
	g_ind = which_m_indxne_i(g_ind, g_Maus, 0);

	// g_Maus[g_ind]<-1
	assegna1_ms_indx_i(g_Maus, g_ind, 1);
	vert = 1;
	// g_colore <- c(1, rep(0, N - 1))
	CREAv_i(g_tmp_i, N - 1);
	InitVett_i(g_tmp_i, 0);
	ASSEGNAv_i(g_scalare_i, 1, 1);
	g_colore = vettore2v_i(g_colore, g_scalare_i, g_tmp_i);
	// dist <- c(0, rep(Inf, N - 1))
	CREAv_d(g_tmp1_d, N - 1);
	InitVett_d(g_tmp1_d, R_PosInf);
	ASSEGNAv_d(g_scalare_d, 1, 0.0);
	ris = vettore2v_d(ris, g_scalare_d, g_tmp1_d);
	ASSEGNAv_i(g_grigi, 1, 1);
	while (LENGTHv_i(g_grigi) != 0) {
		u = ACCEDIv_i(g_grigi, 1);
		//~ g_adj<-which(g_Maus[u,]!=0)
		g_adj = which_m_rowindxne_i(g_adj, g_Maus, u, 0);
		L = LENGTHv_i(g_adj);
		if (L > 0) {
			for (i = 1; i <= L; i++) {
				v = ACCEDIv_i(g_adj, i);
				if (ACCEDIv_i(g_colore, v) == 0)	{
					ASSEGNAv_i(g_colore, v, 1);
					ASSEGNAv_d(ris, v, ACCEDIv_d(ris, u) + 1);
					// g_grigi<-c(g_grigi,v)
					g_grigi = accoda1_vs_i(g_grigi, v);
				}
			}
		}
		lg = LENGTHv_i(g_grigi);
		if (lg >= 2)
			// g_grigi<-g_grigi[2:lg]
			segmento1_v_i(g_grigi, 2, lg);
		else
			break;
	}

	StrBilanciam();

	_Intestazione("\n*** Esco da check_conn ***\n\n");

#ifdef DET
	fprintf(fp_det, "check_conn output:\n");
	fprintf(fp_det, "\tdist = ");
	_StampaRawVett_d(ris);
#endif

	return ris;
}

SEXP check_conn(SEXP Mdiscr)
{
	int nProtected = 0;
	MATRICEi *Mdiscr1;
	VETTOREd *v1 = NULL;
	SEXP ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** check_conn ***\n");

	Mdiscr1 = inMATRICE_i(Mdiscr, &nProtected);

	v1 = check_conn1(v1, Mdiscr1);

	ris = daVETTORE_d(v1, &nProtected);

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	return ris;
}
