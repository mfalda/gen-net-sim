#include "triangola.h"

#define g_tmp1_i globali.triangola.tmp1_i
#define g_tmp2_i globali.triangola.tmp2_i
#define g_tmp3_d globali.triangola.tmp3_d
#define g_M_in globali.triangola.M_in
#define g_sk globali.triangola.sk
#define g_coord_ind globali.triangola.coord_ind
#define g_ind_aus globali.triangola.ind_aus
#define g_ind globali.triangola.ind
#define g_coord globali.triangola.coord
#define g_tmp_coord globali.triangola.tmp_coord
#define g_tmp_coord1 globali.triangola.tmp_coord1
#define g_Dmem globali.triangola.Dmem

// triangola<-function(M,Cg,max.con)
MATRICEi *triangola(MATRICEi *ris, double Cg, int max_con)
// {
{
	int i, S_i, L, Ng, conta, k, aus1, aus2, ng_s, Lr, j, h;
	double ng;

	_Intestazione("\n***triangola***\n");
#ifdef DET
	fprintf(fp_det, "input:\n");
	fprintf(fp_det, "\tM = ");
	_StampaRawMatr_i(ris);
	fprintf(fp_det, "\tCg =  %.16g\n", Cg);
	fprintf(fp_det, "\tmax_con =  %d\n", max_con);
#endif

	CREAv_d(g_tmp3_d, 1);
	//  Ng<-dim(M)[1]
	Ng = LENGTHm1_i(ris);
	//  g_Dmem<-diag(M)
	g_Dmem = diag_i(g_Dmem, ris);
	//  diag(M)<-0
	assegna1_s_diag_i(ris, 0);
	//  for (i in (1:Ng))
	for (i = 1; i <= Ng; i++) {
		//     {g_ind<-union(which(M[i,]!=0),which(M[,i]!=0))  #neighbours
		g_tmp1_i= which_m_rowindxne_i(g_tmp1_i, ris, i, 0);
		g_tmp2_i = which_m_colindxne_i(g_tmp2_i, ris, i, 0);
		g_ind = unione_i(g_ind, g_tmp1_i, g_tmp2_i);  //neighbours;
		//      S.i<-length(g_ind)
		S_i = LENGTHv_i(g_ind);
		//      if (S.i>=2)
		if (S_i >= 2) {
			// 	{ng<-Cg*S.i*(S.i-1)/2  #average number of links among neighbours
			ng = (double) Cg * S_i * (S_i - 1) / 2.0;  //average number of links among neighbours;
			// 	 ng.s<-round(rnorm_s(1,ng,1)) # cioè troncamento
			g_tmp3_d = rnorm_s(g_tmp3_d, 1, ng, 1, "triangola");
			ng_s = rround(ACCEDIv_d(g_tmp3_d, 1), 0);
			// 	 if (ng.s>0)
			if (ng_s > 0) {
				// 	   {L<-length(g_ind)
				L = LENGTHv_i(g_ind);
				// 	    conta<-0
				conta = 0;
				// 	    g_coord<-matrix(0,1,2) // in realta` ha sempre almeno due righe
				CREAm_i(g_coord, 2, 2);
				g_coord->nr = 1;
				InitMatr_i(g_coord, 0);
				//             k<-0
				k = 0;
				// 	    for (j in (1:(L-1)))
				for (j = 1; j <= L - 1; j++) {
					// 	     {for (h in ((j+1):L))
					for (h = j + 1; h <= L; h++) {
						// 		{aus1<-M[g_ind[j],g_ind[h]]
						aus1 = ACCEDIm_i(ris, ACCEDIv_i(g_ind, j), ACCEDIv_i(g_ind, h));
						//                  aus2<-M[g_ind[h],g_ind[j]]
						aus2 = ACCEDIm_i(ris, ACCEDIv_i(g_ind, h), ACCEDIv_i(g_ind, j));
						// 		 if ( (aus1==1)|(aus2==1) ) conta<-conta+1
						if ( (aus1 == 1) || (aus2 == 1) )
							conta++;
						//                   else g_coord<-rbind(g_coord,c(g_ind[j],g_ind[h]),c(g_ind[h],g_ind[j]))
						else {
							// g_coord = rbind(g_coord, vettore2(g_ind[j], g_ind[h]), vettore2(g_ind[h], g_ind[j]));
							g_tmp1_i = vettore2s_i(g_tmp1_i, ACCEDIv_i(g_ind, j), ACCEDIv_i(g_ind, h));
							g_coord = aggiungi_mv_riga_i(g_coord, LENGTHm1_i(g_coord) + 1, g_tmp1_i);
							g_tmp1_i = vettore2s_i(g_tmp1_i, ACCEDIv_i(g_ind, h), ACCEDIv_i(g_ind, j));
							// a questo punto ha già ingrandito la matrice
							g_coord = aggiungi_mv_riga_i(g_coord, LENGTHm1_i(g_coord) + 1, g_tmp1_i);
						}
						// 		}
					}
					// 	     }
				}
				// 	    ng.s<-ng.s-conta
				ng_s -= conta;
				//             Lr<-ng.s
				Lr = ng_s;
				//             while (Lr>0)
				while (Lr > 0) {
					//              {M.in<-apply(M,1,sum)
					g_M_in = somma_righe_i(g_M_in, ris);
					// 	      g_ind<-which(M.in<max.con)
					g_ind = which_v_indxlt_i(g_ind, g_M_in, max_con);
					//               if (length(g_ind)>0)
					if (LENGTHv_i(g_ind) > 0) {
						// 		{g_ind.aus<-which(g_coord[,1]%in%g_ind)
						g_ind_aus = which_m_colindxin_i(g_ind_aus, g_coord, 1, g_ind);
						// 		 k<-length(g_ind.aus)
						k = LENGTHv_i(g_ind_aus);
						// 		 if (k>0)
						if (k > 0) {
							// 		   {g_coord<-matrix(g_coord[g_ind.aus,],ncol=2)
							g_tmp_coord1 = righe_i(g_tmp_coord1, g_coord, g_ind_aus);
							g_coord = copia_m_i(g_coord, g_tmp_coord1);
							// 		    g_sk<-sampleB(seq(1,k,1),1)
							g_tmp1_i = seq_i(g_tmp1_i, 1, k, 1);
							g_sk = sampleB(g_sk, g_tmp1_i, 1, 0);
							//                     g_coord.g_ind<-g_coord[g_sk,]
							g_coord_ind = riga_i(g_coord_ind, g_coord, ACCEDIv_i(g_sk, 1));
							// 		    M[g_coord.g_ind[1],g_coord.g_ind[2]]<-1
							ASSEGNAm_i(ris, ACCEDIv_i(g_coord_ind, 1), ACCEDIv_i(g_coord_ind, 2), 1);
							//                     g_ind<-which((g_coord[, 1] != g_coord.g_ind[1])&(g_coord[, 1] != g_coord.g_ind[2]))
							g_ind = which_m_colneand2_i(g_ind, g_coord, 1, ACCEDIv_i(g_coord_ind, 1), ACCEDIv_i(g_coord_ind, 2));
							// 		    if (length(g_ind)>0) {g_coord<-matrix(g_coord[g_ind,],ncol=2); Lr<-Lr-1}
							if (LENGTHv_i(g_ind) > 0) {
								g_tmp_coord1 = righe_i(g_tmp_coord1, g_coord, g_ind);
								g_coord = copia_m_i(g_coord, g_tmp_coord1);
								Lr--;
							}
							// 		     else Lr<-0
							else
								Lr = 0;
							// 		   }
						}
						//                   else Lr<-0
						else
							Lr = 0;
						// 		}
					}
					// 	       else Lr<-0
					else
						Lr = 0;
					//              }
				}
				// 	   }
			}
			// 	}
		}
		//     }
	}
	//~ CANCELLAv_i(g_tmp1_i);
	//~ CANCELLAv_i(g_tmp2_i);
	//~ CANCELLAv_i(g_M_in);
	//~ CANCELLAv_i(g_sk);
	//~ CANCELLAv_i(g_coord_ind);
	//~ CANCELLAv_i(g_ind_aus);
	//~ CANCELLAv_i(g_ind);
	//~ CANCELLAm_i(g_coord);
	//~ CANCELLAm_i(g_tmp_coord);
	//~ CANCELLAm_i(g_tmp_coord1);

	// diag(M)<-Dmem
	assegna1_v_diag_i(ris, g_Dmem);
	//~ CANCELLAv_i(g_Dmem);

	StrBilanciam();

#ifdef DET
	fprintf(fp_det, "triangola output:\n\tM = ");
	_StampaRawMatr_i(ris);
#endif

	// return(M)
	return ris;
	// }
}

