#include "pathlength.h"

#define g_dist_tot globali.pathlength.dist_tot
#define g_new_neighbours globali.pathlength.new_neighbours
#define g_distanza globali.pathlength.distanza
#define g_visitati globali.pathlength.visitati
#define g_neighbours globali.pathlength.neighbours
#define g_tmp1_i globali.pathlength.tmp1_i
#define g_tmpm_i globali.pathlength.tmpm_i

// pathlength<-function(W,und=FALSE)
MATRICEd *pathlength(MATRICEd *ris, const MATRICEi *W, bool und/* = FALSE*/)
{
	int N, j, i, ind, L;
	double conta;

	_Intestazione("\n*** pathlength ***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tW = ");
	_StampaRawMatr_i(W);
	if (und)
		fprintf(fp_det, "\tund = TRUE\n");
	else
		fprintf(fp_det, "\tund = FALSE\n");
#endif

// N<-dim(W)[1]
	N = LENGTHm1_i(W);
//  if (und==TRUE)
	if (und) {
		// W<-W+t(W)
		g_tmpm_i = trasponi_i(g_tmpm_i, W);
		somma1_m_i(g_tmpm_i, W);
	}
	else {
		g_tmpm_i = copia_m_i(g_tmpm_i, W);
	}
//  W[which(W!=0,arr.ind=TRUE)]<-1
	g_tmp1_i = which_m_indxne_i(g_tmp1_i, g_tmpm_i, 0);
	assegna1_ms_indx_i(g_tmpm_i, g_tmp1_i, 1);
//  M<-matrix(NA,N,N)
	CREAm_d(ris, N, N);
	InitMatr_d(ris, NA_REAL);
//  g_dist.tot<-vector()
	CREAv_i(g_dist_tot, 0);
//  for (i in (1:N))
	for (i = 1; i <= N; i++) {
//   g_distanza<-rep(NA,N)
		g_distanza = rep_s_d(g_distanza, NA_REAL, N);
//    #distanza[i]<-0
		//g_distanza[i]=0;
//    g_visitati<-vector()  #visitati<-i
		CREAv_i(g_visitati, 0);
//    conta<-1
		conta = 1.0;
//    g_neighbours<-which(W[,i]!=0)#setdiff(which(W[,i]!=0),i)
		g_neighbours = which_m_colindxne_i(g_neighbours, g_tmpm_i, i, 0);
//setdiff(which(W[,i]!=0),i);
//    L<-length(g_neighbours)
		L = LENGTHv_i(g_neighbours);
//    while (L>0)
		while (L > 0) {
//      g_distanza[g_neighbours]<-conta
			assegna1_vs_indx_d(g_distanza, g_neighbours, conta);
//       #distanza[i]<-0
			//g_distanza[i]=0;
//       g_visitati<-union(g_visitati,g_neighbours)
			g_visitati = unione_i(g_visitati, g_visitati, g_neighbours);
//       conta<-conta+1
			conta++;
//       new.g_neighbours<-vector()
			CREAv_i(g_new_neighbours, 0);
//       for (j in (1:L))
			for (j = 1; j <= L; j++) {
				//          {ind<-g_neighbours[j]
				ind = ACCEDIv_i(g_neighbours, j);
// 	        new.g_neighbours<-c(new.g_neighbours,which(W[,ind]!=0))
				g_tmp1_i = which_m_colindxne_i(g_tmp1_i, W, ind, 0);
				g_new_neighbours = unione1_i(g_new_neighbours, g_tmp1_i);
//          }
			}
//       new.g_neighbours<-union(new.g_neighbours,new.g_neighbours)
			g_new_neighbours = unione1_i(g_new_neighbours, g_new_neighbours);
//       new.g_neighbours<-setdiff(new.g_neighbours,g_visitati)
			setdiff1_i(g_new_neighbours, g_visitati);
//       g_neighbours<-new.g_neighbours
			g_neighbours = copia_v_i(g_neighbours, g_new_neighbours, 1, LENGTHv_i(g_new_neighbours));
//       L<-length(g_neighbours)
			L = LENGTHv_i(g_neighbours);
//     }
		}
//   M[,i]<-g_distanza
		assegna1_mv_colonna_d(ris, i, g_distanza);
//  }
	}
// if (und==TRUE) diag(M)<-0
	if (und)
		assegna1_s_diag_i(g_tmpm_i, 0);
// M[which(is.na(M),arr.ind=TRUE)]<-Inf
	assegna1_m_indxNA_d(ris, ris, R_PosInf, false);

	StrBilanciam();

	_Intestazione("\n*** Esco da pathlength ***\n");

#ifdef DET
	fprintf(fp_det, "pathlength output:\n\tM = ");
	_StampaRawMatr_d(ris);
#endif

// return(M)
	return ris;
// }
}
