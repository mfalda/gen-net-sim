#include "hubs.h"

#define g_S globali.hubs.S
#define g_P globali.hubs.P
#define g_ind globali.hubs.ind
#define g_ind2 globali.hubs.ind2
#define g_tmp1_i globali.hubs.tmp1_i
#define g_tmp1_d globali.hubs.tmp1_d
#define g_Hio globali.hubs.Hio
#define g_H globali.hubs.H
#define g_pl globali.hubs.pl

// hubs<-function(mod)
LISTA *hubs(LISTA *ris, const MATRICEi *mod, bool *fb)
{
//hubs defined based on out-degree and ability to reach other nodes
	bool feedback;
	enum TIPO tipi[2];

	_Intestazione("\n***hubs***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tmod = ");
	_StampaRawMatr_i(mod);
#endif

//  g_S<-apply(mod,2,sum)
	g_S = somma_colonne_i(g_S, mod);
//  g_ind<-which(g_S==max(g_S))
	g_ind = which_v_indxeq_i(g_ind, g_S, max_v_i(g_S));
//  g_pl<-pathlength(mod)
	g_pl = pathlength(g_pl, mod, false);
//  if ( min(diag(g_pl))!=Inf )  feedback<-TRUE
	g_tmp1_d = diag_d(g_tmp1_d, g_pl);
	if (min_v_d(g_tmp1_d) != R_PosInf)
		feedback = true;
//  else feedback<-FALSE
	else
		feedback = false;
//  if (length(g_ind)==1) g_H<-ind
	if (LENGTHv_i(g_ind) == 1)
		g_H = copia_v_i(g_H, g_ind, 1, LENGTHv_i(g_ind));
//  else
	else {
//    {diag(g_pl)<-0
		assegna1_s_diag_d(g_pl, 0.0);
//     g_P<-apply(g_pl,2,sum)[g_ind]
		g_tmp1_d = somma_colonne_d(g_tmp1_d, g_pl);
		g_P = copia_v_indx_d(g_P, g_tmp1_d, g_ind);
//     if (min(g_P)!=Inf) {g_ind2<-which(g_P==min(g_P)); g_H<-g_ind[g_ind2]}
		if (min_v_d(g_P) != R_PosInf) {
			g_ind2 = which_v_indxeq_d(g_ind2, g_P, min_v_d(g_P));
			g_H = copia_v_indx_i(g_H, g_ind, g_ind2);
		}
//     else g_H<- ind
		else
			g_H = copia_v_i(g_H, g_ind, 1, LENGTHv_i(g_ind));
//    }
	}
//  #hubs defined based on global-degree and ability to reach other nodes
//hubs defined based on global-degree and ability to reach other nodes;
//  g_S<-g_S+apply(mod,1,sum)
	g_tmp1_i = somma_righe_i(g_tmp1_i, mod);
	somma1_vv_i(g_S, g_tmp1_i);
//  g_ind<-which(g_S==max(g_S))
	g_ind = which_v_indxeq_i(g_ind, g_S, max_v_i(g_S));
//  g_Hio<-ind
	g_Hio = copia_v_i(g_Hio, g_ind, 1, LENGTHv_i(g_ind));

//  return(list(g_H,feedback,g_Hio))
	tipi[0] = VETTi;
	tipi[1] = VETTi;
	CreaLISTA(ris, tipi, 2);
	ris->dati[0].vi = g_H;
	ris->dati[1].vi = g_Hio;
	*fb = feedback;

	StrBilanciam();

	_Intestazione("\n*** Esco da hubs ***\n\n");

#ifdef DET
	fprintf(fp_det, "hubs output:\n");
	fprintf(fp_det, "\tlist(g_H,feedback,g_Hio) = ");
	_StampaRawVett_i(g_H);
	if (*fb)
		fprintf(fp_det, "feedback =  TRUE\n");
	else
		fprintf(fp_det, "feedback =  FALSE\n");
	_StampaRawVett_i(g_Hio);
#endif

	return ris;
// }
}
