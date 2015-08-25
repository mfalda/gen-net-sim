#include "HMM.h"

#define g_Sr globali.HMM.Sr
#define g_Sin globali.HMM.Sin
#define g_Sout globali.HMM.Sout
#define g_hubs globali.HMM.hubs
#define g_ind_Sc globali.HMM.ind_Sc
#define g_LG globali.HMM.LG
#define g_ind globali.HMM.ind
#define g_ind_s globali.HMM.ind_s
#define g_indm globali.HMM.indm
#define g_CCs globali.HMM.CCs
#define g_non_connessi globali.HMM.non_connessi
#define g_h globali.HMM.h
#define g_indp globali.HMM.indp
#define g_h_new globali.HMM.h_new
#define g_m globali.HMM.m
#define g_MC globali.HMM.MC
#define g_tmp1_i globali.HMM.tmp1_i
#define g_scalare_i globali.HMM.scalare_i
#define g_CC globali.HMM.CC
#define g_FB globali.HMM.FB
#define g_segno globali.HMM.segno
#define g_h1 globali.HMM.h1
#define g_h_new1 globali.HMM.h_new1
#define g_tmp2_i globali.HMM.tmp2_i
#define g_dim globali.HMM.dim
#define g_mod_type globali.HMM.mod_type
#define g_a1 globali.HMM.a1
#define g_a2 globali.HMM.a2
#define g_toll globali.HMM.toll
#define g_prior_p_subnet globali.HMM.prior_p_subnet
#define g_Freq_in globali.HMM.Freq_in
#define g_Freq_out globali.HMM.Freq_out
#define g_tmp1_d globali.HMM.tmp1_d
#define g_tmp2_d globali.HMM.tmp2_d
#define g_scalare_d globali.HMM.scalare_d
#define g_scalare1_d globali.HMM.scalare1_d
#define g_Prob globali.HMM.Prob
#define g_STin globali.HMM.STin
#define g_STout globali.HMM.STout
#define g_p globali.HMM.p
#define g_dist globali.HMM.dist
#define g_CCp globali.HMM.CCp
#define g_ppp globali.HMM.ppp
#define g_ris_cc globali.HMM.ris_cc
#define g_Sc_vett globali.HMM.Sc_vett
#define g_pm globali.HMM.pm
#define g_prob_mod globali.HMM.prob_mod
#define g_p_out globali.HMM.p_out
#define g_M globali.HMM.M
#define g_mod1 globali.HMM.mod1
#define g_mod2 globali.HMM.mod2
#define g_mod3 globali.HMM.mod3
#define g_tmpm_i globali.HMM.tmpm_i
#define g_aus globali.HMM.aus
#define g_Sc1 globali.HMM.Sc1
#define g_Sc2 globali.HMM.Sc2
#define g_Sc3 globali.HMM.Sc3
#define g_Sc globali.HMM.Sc
#define g_LC globali.HMM.LC
#define g_CCind globali.HMM.CCind

MATRICEi *HMM1(MATRICEi *ris, int N, double Cf_cl, double gamma, const VETTOREd *degree, const GString *INdegree, const Mod *modules, int L, const VETTOREd *prior_p_subnet1, int max_con, bool feedback, bool zero_nodes_with_indegree0, bool sepgraph, double r_tol, double a_tol, int iter) {
	bool autoreg, controllo, controllo2;
	int num, ri, it, Ng, ccc, i, Nrim, k, m1;
	double Cf_cl_int, mx, Mcc, Pm1, Pm2, Pm3, mn;
	struct RisProbM ris_pm1, ris_pm2, ris_pm3;
	GString *tmp_s = NULL, *etich = NULL;
	const char *OUT = "out", *FREE = "free";

	_Intestazione("\n*** HMM1 ***\n");

	CREAv_i(g_scalare_i, 1);
	CREAv_d(g_scalare_d, 1);
	CREAv_d(g_scalare1_d, 1);
	CREAstr(etich, "");
	CREAstr(tmp_s, "");
	g_prior_p_subnet = copia_v_d(g_prior_p_subnet, prior_p_subnet1, 1, LENGTHv_d(prior_p_subnet1));
//         L<-length(MODULES)
	// da parametro
//         g_CC<-g_FB<-AU<-g_MC<-DIM<-rep(0,L)
	// g_FB = rep_s_i(g_FB, 0, L);  // g_FB non viene mai modificato, quindi non mi serve una copia
	// AU = rep_s_i(AU, 0, L);  // AU non viene mai modificato, quindi non mi serve una copia
	CREAv_i(g_CC, L);
	CREAv_i(g_FB, L);
	CREAv_i(g_dim, L);
	g_MC = rep_s_i(g_MC, 0, L); // g_MC e` un valore derivato dal campo "rete"
	// DIM = rep_s_i(DIM, 0, L); // DIM non viene mai modificato, quindi non mi serve una copia
	autoreg = false;
//         for (i in (1:L))
//            g_MC[i]<-max(apply(MODULES[[i]]$net,1,sum))  #max in.degree
	for (i = 0; i < L; i++) {
		g_tmp1_i = somma_righe_i(g_tmp1_i, modules[i].rete);
		ASSEGNAv_i(g_CC, i + 1, modules[i].CC);
		ASSEGNAv_i(g_FB, i + 1, modules[i].feedback);
		ASSEGNAv_i(g_dim, i + 1, modules[i].dim_m);
		ASSEGNAv_i(g_MC, i + 1, max_v_i(g_tmp1_i));
		if (modules[i].autoreg)
			autoreg = true;
	}
	//         if (max.con<g_dim(MODULES[[L]]$net)[1])
	if (max_con < LENGTHm1_i(modules[L - 1].rete)) {
//           {prior.g_p.subnet[which(g_MC>max.con)]<-0
		g_tmp1_i = which_v_indxeq_i(g_tmp1_i, g_MC, max_con);
		assegna1_v_indx_d(g_prior_p_subnet, g_tmp1_i, 0.0);
//            prior.g_p.subnet<-prior.g_p.subnet/sum(prior.g_p.subnet)
		dividi1_vs_d(g_prior_p_subnet, somma_v_d(g_prior_p_subnet, false));
//           }
	}
	if (!feedback) {
//           {prior.g_p.subnet[which(g_FB==TRUE)]<-0
		g_tmp1_i = which_v_indxeq_i(g_tmp1_i, g_FB, true);
		assegna1_v_indx_d(g_prior_p_subnet, g_tmp1_i, 0.0);
//            prior.g_p.subnet<-prior.g_p.subnet/sum(prior.g_p.subnet)
		dividi1_vs_d(g_prior_p_subnet, somma_v_d(g_prior_p_subnet, false));
//           }
	}
//  g_LG<-rep(-1,N)
	g_LG = rep_s_i(g_LG, -1, N);
//  Mdiscr<-matrix(0,ncol=N,nrow=N)
	CREAm_i(ris, N, N);
	InitMatr_i(ris, 0);
//  if (max.con>N) max.con<-N
	if (max_con > N)
		max_con = N;
//  Freq.out<-Freq.in<-c(N,rep(0,N+1))
	CREAv_d(g_scalare_d, 1);
	ASSEGNAv_d(g_scalare_d, 1, (double) N);
	g_tmp1_d = rep_s_d(g_tmp1_d, 0.0, N + 1);
	g_Freq_out = vettore2v_d(g_Freq_out, g_scalare_d, g_tmp1_d);
	g_Freq_in = copia_v_d(g_Freq_in, g_Freq_out, 1, LENGTHv_d(g_Freq_out));
//  if (is.null(DEGREE))
	if (degree == NULL) {
//   {g_Prob<-c(NA,seq(1,N,1)^(-gamma),0)
		g_tmp1_d = seq_d(g_tmp1_d, 1.0, (double) N, 1.0);
		g_tmp2_d = exp_d(g_tmp2_d, g_tmp1_d, -gamma);
		ASSEGNAv_d(g_scalare1_d, 1, NA_REAL);
		ASSEGNAv_d(g_scalare_d, 1, 0.0);
		g_Prob = vettore3v_d(g_Prob, g_scalare1_d, g_tmp2_d, g_scalare_d);
//    g_Prob<-g_Prob/(sum(g_Prob,na.rm=TRUE))
		dividi1_vs_d(g_Prob, somma_v_d(g_Prob, true));
//   }
	}
//  else g_Prob<-c(DEGREE/(sum(DEGREE,na.rm=TRUE)),0)
	else {
		g_tmp1_d = dividi_vs_d(g_tmp1_d, degree, somma_v_d(degree, true));
		CREAv_d(g_Prob, LENGTHv_d(g_tmp1_d) + 1);
		// la funzione copia non modifica la lunghezza!
		g_Prob = copia_v_d(g_Prob, g_tmp1_d, 1, LENGTHv_d(g_tmp1_d));
		accoda1_vs_d(g_Prob, 0.0);
	}
//  g_STout<-g_Prob*N
	g_STout = moltiplica_vs_d(g_STout, g_Prob, N);
//  g_STin<-rep(0,N+2)
	g_STin = rep_s_d(g_STin, 0.0, N + 2);
//  if (INdegree=="out")  g_STin[1:(max.con+1)]<-N*g_Prob[1:(max.con+1)]/sum(g_Prob[1:(max.con+1)],na.rm=TRUE)
	if (!strncmp(INdegree->str, OUT, 3)) {
		g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 1, max_con + 1);
		g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
		// in somma cancNA = 1
		dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, true));
		assegna1_v_segmv_d(g_STin, 1, max_con + 1, g_tmp2_d);
	}
//  else   g_STin[1:(max.con+1)]<-NA
	else
		assegna1_v_segm_d(g_STin, 1, max_con + 1, NA_REAL);
//  if (zero.nodes.with.indegree0==TRUE)
	if (zero_nodes_with_indegree0) {
//    {if ((g_STin[1]!=0)&(!is.na(g_STin[1])))
		if (!Uguale(ACCEDIv_d(g_STin, 1), 0.0) && !ISNA(ACCEDIv_d(g_STin, 1))) {
//       {cat("\n Parameter 'zero.nodes.with.indegree0' was set to TRUE by the user; therefore, all nodes are required to have in.degree>=1")
			Rprintf("Parameter 'zero_nodes_with_indegree0' was set to TRUE by the user; therefore, all nodes are required to have in_degree >= 1");
//        g_STin[1]<-0
			ASSEGNAv_d(g_STin, 1, 0.0);
//        g_STin[2:(max.con+1)]<-N*g_Prob[2:(max.con+1)]/sum(g_Prob[2:(max.con+1)],na.rm=TRUE)   # (INdegree=="out")
			g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 2, max_con + 1);
			g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
			// in somma cancNA = 1
			dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, true));
			assegna1_v_segmv_d(g_STin, 2, max_con + 1, g_tmp2_d);
//       }
		}
//      g_STin[1]<-0 #if NA
		ASSEGNAv_d(g_STin, 1, 0.0); //if NA;
//    }
	}
//  g_aus<-cbind(g_STout*r.tol,rep(a.tol,length(g_STout)))
	g_tmp1_d = moltiplica_vs_d(g_tmp1_d, g_STout, r_tol);
	g_tmp2_d = rep_s_d(g_tmp2_d, a_tol, LENGTHv_d(g_STout));
	g_aus = cbind2v_d(g_aus, g_tmp1_d, g_tmp2_d);
//  g_toll<-apply(g_aus,1,max)
	g_toll = max_righe_d(g_toll, g_aus);
//  g_h<-seq(1,N,1)
	g_h = seq_i(g_h, 1, N, 1);
//  g_h.new<-vector()
	CREAv_i(g_h_new, 0);
//  Cf.cl.int<-Cf.cl
	Cf_cl_int = Cf_cl;
//  it<-0
	it = 0;
//  g_CCp<-rep(1,L)
	g_CCp = rep_s_d(g_CCp, 1.0, L);
	// anticipo la creazione di g_CCind
	g_CCind = mia_alloc(L, VETTOREi*);
	if (L > 0 && g_CCind == NULL) {
		Rprintf("Not enough memory (HMM # %d, g_CCind)", __LINE__ - 2);
		error("");
	}
	for (i = 0; i < L; i++)
		g_CCind[i] = NULL;
//  if (Cf.cl!=0)
	if (!Uguale(Cf_cl, 0.0)) {
//    {g_CCs<-sort(union(g_CC,g_CC))
		g_tmp1_i = elimina_doppi_i(g_tmp1_i, g_CC);
		g_CCs = ordina_i(g_CCs, g_tmp1_i, false);
//     g_CCind<-list()
//     Lc<-length(g_CCs)
		if (g_LC == 0)
			g_LC = LENGTHv_i(g_CCs) - 1;
		if (g_LC < LENGTHv_i(g_CCs)) {
			g_LC = max_s_i(g_LC, LENGTHv_i(g_CCs));
		}
//     g_p<-rep(1,Lc)
		g_p = rep_s_d(g_p, 1.0, g_LC);
//     for (i in (1:Lc)) g_CCind[[i]]<-which(g_CC==g_CCs[i])
		for (i = 1; i <= g_LC; i++) {
			g_tmp1_i = which_v_indxeq_i(g_tmp1_i, g_CC, ACCEDIv_i(g_CCs, i));
			g_CCind[i - 1] = copia_v_i(g_CCind[i - 1], g_tmp1_i, 1, LENGTHv_i(g_tmp1_i));
		}
//    }
	}
//  while (length(g_h)>1)
	while (LENGTHv_i(g_h) > 1) {
//    {Nrim<-length(g_h)
		Nrim = LENGTHv_i(g_h);
//     g_Sin<-apply(Mdiscr,1,sum)
		g_Sin = somma_righe_i(g_Sin, ris);
//     g_Sout<-apply(Mdiscr,2,sum)
		g_Sout = somma_colonne_i(g_Sout, ris);
//     #################
		//################;
//     Mcc<-cluster.coeff(Mdiscr)[[1]] #mean(g_CC[g_ind.possible])
		g_ris_cc = cluster_coeff2(g_ris_cc, ris, &Mcc); //mean(g_CC[ind_possible]);
//     if (Cf.cl>0)
		if (Cf_cl > 0) {
//        {if (Cf.cl.int!=Mcc)
			if (!Uguale(Cf_cl_int, Mcc)) {
//           {g_segno<-sign(g_CCs-Mcc)
				g_tmp1_i = somma_vs_i(g_tmp1_i, g_CCs, -Mcc);
				g_segno = segno_v_i(g_segno, g_tmp1_i);
//            g_indp<-which(g_segno>0)
				g_indp = which_v_indxgt_i(g_indp, g_segno, 0) ;
//            g_indm<-which(g_segno<=0)
				g_indm = which_v_indxle_i(g_indm, g_segno, 0) ;
//            if (Cf.cl.int>Mcc)
				if (Cf_cl_int > Mcc) {
//              {g_p[g_indp]<-(1-abs(Cf.cl.int-g_CCs[g_indp]))+g_p[g_indp]
					g_p = f_aux9_d(g_p, 1, g_CCs, g_indp, Cf_cl_int);
//               g_p[g_indm]<-(-abs(Cf.cl.int-g_CCs[g_indm]))+g_p[g_indm]
					g_p = f_aux9_d(g_p, 1, g_CCs, g_indm, Cf_cl_int);
//              }
				}
//            else
				else {
//              g_p[g_indm]<-(1-abs(Cf.cl.int-g_CCs[g_indm]))+g_p[g_indm]
					g_p = f_aux9_d(g_p, 1, g_CCs, g_indm, Cf_cl_int);
//               g_p[g_indp]<-(-abs(Cf.cl.int-g_CCs[g_indp]))+g_p[g_indp]
					g_p = f_aux9_d(g_p, 0, g_CCs, g_indp, Cf_cl_int);
//              }
				}
//            if (min(g_p)<=0) g_p[which(g_p<0)]<-0.0001 #if (min(g_p)<=0) g_p<-g_p-min(g_p)+0.001
				if (min_v_d(g_p) <= 0.0)
					assegna1_v_indxle_d(g_p, g_p, 0.0, 0.0001);
//            g_p<-g_p/sum(g_p)
				dividi1_vs_d(g_p, somma_v_d(g_p, false));
//            for (ccc in (1:Lc)) g_CCp[(g_CCind[[ccc]])]<-g_p[ccc]
				for (ccc = 1; ccc <= g_LC; ccc++)
					assegna1_v_indx_d(g_CCp, g_CCind[ccc - 1], ACCEDIv_d(g_p, ccc));
//            }
			}
//         else
			else {
//           g_p<-g_p+(1-abs(Cf.cl.int-g_CCs))
				g_tmp1_i = seq_i(g_tmp1_i, 1, LENGTHv_d(g_p), 1);
				g_p = f_aux9_d(g_p, 1, g_CCs, g_tmp1_i, Cf_cl_int);
//            g_p<-g_p/sum(g_p)
				dividi1_vs_d(g_p, somma_v_d(g_p, false));
//            for (ccc in (1:Lc)) g_CCp[(g_CCind[[ccc]])]<-g_p[ccc]
				for (ccc = 1; ccc <= g_LC; ccc++)
					assegna1_v_indx_d(g_CCp, g_CCind[ccc - 1], ACCEDIv_d(g_p, ccc));
//           }
			}
//        }
		}
//     g_ppp<-prior.g_p.subnet*CCp
		g_ppp = moltiplica_vv_d(g_ppp, g_prior_p_subnet, g_CCp);
//     if ((max(g_ppp)==0)|(is.na(max(g_ppp))))   stop("\n prior.g_p.subnet values incompatible with other parameter settings! \n")
		mx = max_v_d(g_ppp);
		if (Uguale(mx, 0.0) || ISNA(mx))
			error("g_prior_p_subnet values incompatible with other parameter settings! \n");
//     g_ppp<-g_ppp/sum(g_ppp,na.rm=TRUE)
		dividi1_vs_d(g_ppp, somma_v_d(g_ppp, true));
//     controllo<-1
		controllo = true;
//     controllo2<-1
		controllo2 = true;

//     if (length(g_h)<max(DIM))
		if (LENGTHv_i(g_h) < max_v_i(g_dim)) {
			// g_ppp[which(DIM>length(g_h))]<-0
			g_tmp1_i = which_v_indxgt_i(g_tmp1_i, g_dim, LENGTHv_i(g_h));
			assegna1_v_indx_d(g_ppp, g_tmp1_i, 0.0);
		}
//     if ((max(g_ppp)==0)|(is.na(max(g_ppp)))) {controllo<-0; controllo2<-0}
		mx = max_v_d(g_ppp);
		if (Uguale(mx, 0.0) || ISNA(mx)) {
			controllo = false;
			controllo2 = false;
		}
//     while (controllo==1)
		while (controllo) {
//     k<-min(3,length(which(g_ppp>0)))
			g_tmp1_i = which_v_indxgt_d(g_tmp1_i, g_ppp, 0.0);
			k = min_s_i(3, LENGTHv_i(g_tmp1_i));
//       g_m<-sampleB(1:L,k,prob=g_ppp)
			g_tmp1_i = seq_i(g_tmp1_i, 1, L, 1);
			g_m = sampleB_p(g_m, g_tmp1_i, k, false, g_ppp);
//       g_mod1<-MODULES[[(g_m[1])]]$net
			g_mod1 = copia_m_i(g_mod1, modules[ACCEDIv_i(g_m, 1) - 1].rete);
//       aus1<-probmod(g_M=g_mod1,g_h=g_h,g_Sin=g_Sin,g_Sout=g_Sout,g_STin=g_STin,g_STout=g_STout,Freq.in=Freq.in,Freq.out=Freq.out,g_toll=g_toll)
			g_Sc1 = probmod2(g_Sc1, g_mod1, g_h, g_Sin, g_Sout, g_STin, g_STout, g_Freq_in, g_Freq_out, g_toll, &ris_pm1);
//       Pm1<-aus1[[1]]
			Pm1 = ris_pm1.score;
//       if (k>1)
			if (k > 1) {
//         g_mod2<-MODULES[[(g_m[2])]]$net
				g_mod2 = copia_m_i(g_mod2, modules[ACCEDIv_i(g_m, 2) - 1].rete);
//          aus2<-probmod(g_M=g_mod2,g_h=g_h,g_Sin=g_Sin,g_Sout=g_Sout,g_STin=g_STin,g_STout=g_STout,Freq.in=Freq.in,Freq.out=Freq.out,g_toll=g_toll)
				g_Sc2 = probmod2(g_Sc2, g_mod2, g_h, g_Sin, g_Sout, g_STin, g_STout, g_Freq_in, g_Freq_out, g_toll, &ris_pm2);
//          Pm2<-aus2[[1]]
				Pm2 = ris_pm2.score;
//         }
			}
//       else Pm2<-NA
			else
				Pm2 = NA_REAL;
//       if (k>2)
			if (k > 2) {
//        g_mod3<-MODULES[[(g_m[3])]]$net

				g_mod3 = copia_m_i(g_mod3, modules[ACCEDIv_i(g_m, 3) - 1].rete);
//          aus3<-probmod(g_M=g_mod3,g_h=g_h,g_Sin=g_Sin,g_Sout=g_Sout,g_STin=g_STin,g_STout=g_STout,Freq.in=Freq.in,Freq.out=Freq.out,g_toll=g_toll)
				g_Sc3 = probmod2(g_Sc3, g_mod3, g_h, g_Sin, g_Sout, g_STin, g_STout, g_Freq_in, g_Freq_out, g_toll, &ris_pm3);
//          Pm3<-aus3[[1]]
				Pm3 = ris_pm3.score;
//         }
			}
//       else Pm3<-NA
			else
				Pm3 = NA_REAL;

//       g_pm<-prob.mod<-c(Pm1,Pm2,Pm3)
			g_pm = vettore3s_d(g_pm, Pm1, Pm2, Pm3);
			g_prob_mod = copia_v_d(g_prob_mod, g_pm, 1, LENGTHv_d(g_pm));
//       prob.mod[which(is.na(prob.mod))]<-0
			assegna1_v_indxNA_d(g_prob_mod, g_prob_mod, 0.0, false);
//       mn<-min(prob.mod)
			mn = min_v_d(g_prob_mod);
//       if (mn<=0) prob.mod<-prob.mod-mn+1/9
			if (mn <= 0.0)
				somma1_vs_d(g_prob_mod, 1.0 / 9.0 - mn);
//       prob.mod[which(is.na(g_pm))]<-0
			g_tmp1_i = which_v_indxNA_d(g_tmp1_i, g_pm, false);
			assegna1_v_indx_d(g_prob_mod, g_tmp1_i, 0.0);
// 	    g_ind<-which(prob.mod!=0)
			g_ind = which_v_indxne_d(g_ind, g_prob_mod, 0.0);
//       if (length(g_ind>0)) ERRORE???
			if (LENGTHv_i(g_ind) > 0) {
//         prob.mod<-prob.mod/sum(prob.mod)
				dividi1_vs_d(g_prob_mod, somma_v_d(g_prob_mod, false));
//          mod.type<-sampleB(c(1,2,3),1,prob=prob.mod)
				g_tmp1_i = vettore3s_i(g_tmp1_i, 1, 2, 3);
				g_mod_type = sampleB_p(g_mod_type, g_tmp1_i, 1, false, g_prob_mod);
//          if (mod.type==1) {g_M<-aus1[[2]]; g_Sc<-aus1[[3]];  if (INdegree=="free") g_hubs<-MODULES[[(g_m[1])]]$g_hubs else g_hubs<-MODULES[[(g_m[1])]]$hubsio}
				if (ACCEDIv_i(g_mod_type, 1) == 1) {
					// probmod NON modifica M
					g_M = copia_m_i(g_M, g_mod1);
					g_Sc = copia_m_d(g_Sc, g_Sc1);
					if (!strncmp(INdegree->str, FREE, 4))
						g_hubs = copia_v_i(g_hubs, modules[ACCEDIv_i(g_m, 1) - 1].hubs, 1, LENGTHv_i(modules[ACCEDIv_i(g_m, 1) - 1].hubs));
					else
						g_hubs = copia_v_i(g_hubs, modules[ACCEDIv_i(g_m, 1) - 1].hubsio, 1, LENGTHv_i(modules[ACCEDIv_i(g_m, 1) - 1].hubsio));
				}
// 	       if (mod.type==2) {g_M<-aus2[[2]]; g_Sc<-aus2[[3]];  if (INdegree=="free") g_hubs<-MODULES[[(g_m[2])]]$g_hubs else g_hubs<-MODULES[[(g_m[2])]]$hubsio}
				else if (ACCEDIv_i(g_mod_type, 1) == 2) {
					// probmod NON modifica M
					g_M = copia_m_i(g_M, g_mod2);
					g_Sc = copia_m_d(g_Sc, g_Sc2);
					if (!strncmp(INdegree->str, FREE, 4))
						g_hubs = copia_v_i(g_hubs, modules[ACCEDIv_i(g_m, 2) - 1].hubs, 1, LENGTHv_i(modules[ACCEDIv_i(g_m, 2) - 1].hubs));
					else
						g_hubs = copia_v_i(g_hubs, modules[ACCEDIv_i(g_m, 2) - 1].hubsio, 1, LENGTHv_i(modules[ACCEDIv_i(g_m, 2) - 1].hubsio));
				}
// 	       if (mod.type==3) {g_M<-aus3[[2]]; g_Sc<-aus3[[3]];  if (INdegree=="free") g_hubs<-MODULES[[(g_m[3])]]$g_hubs else g_hubs<-MODULES[[(g_m[3])]]$hubsio}
				else if (ACCEDIv_i(g_mod_type, 1) == 3) {
					// probmod NON modifica M
					g_M = copia_m_i(g_M, g_mod3);
					g_Sc = copia_m_d(g_Sc, g_Sc3);
					if (!strncmp(INdegree->str, FREE, 4))
						g_hubs = copia_v_i(g_hubs, modules[ACCEDIv_i(g_m, 3) - 1].hubs, 1, LENGTHv_i(modules[ACCEDIv_i(g_m, 3) - 1].hubs));
					else
						g_hubs = copia_v_i(g_hubs, modules[ACCEDIv_i(g_m, 3) - 1].hubsio, 1, LENGTHv_i(modules[ACCEDIv_i(g_m, 3) - 1].hubsio));
				}
//          controllo<-0
				controllo = false;
//          if (is.na(Pm1))   g_ppp[[(g_m[1])]]<-0
				if (ISNA(Pm1))
					ASSEGNAv_d(g_ppp, ACCEDIv_i(g_m, 1), 0.0);
//          if ((is.na(Pm2))&(k>1))   g_ppp[[(g_m[2])]]<-0
				if ((ISNA(Pm2)) && (k > 1))
					ASSEGNAv_d(g_ppp, ACCEDIv_i(g_m, 2), 0.0);
//          if ((is.na(Pm3))&(k>2))   g_ppp[[(g_m[3])]]<-0
				if ((ISNA(Pm3)) && (k > 2))
					ASSEGNAv_d(g_ppp, ACCEDIv_i(g_m, 3), 0.0);
//          g_ppp<-g_ppp/sum(g_ppp)
				dividi1_vs_d(g_ppp, somma_v_d(g_ppp, false));
//         }
			}
//       else
			else {
//         g_ppp[g_m]<-0
				assegna1_v_indx_d(g_ppp, g_m, 0.0);
//          g_ppp<-g_ppp/sum(g_ppp)
				dividi1_vs_d(g_ppp, somma_v_d(g_ppp, false));
//          if ((max(g_ppp)==0)|(is.na(max(g_ppp))))
				mx = max_v_d(g_ppp);
				if (Uguale(mx, 0.0) || ISNA(mx)) {
//           if (Cf.cl.int>0)
					if (Cf_cl_int > 0) {
//             Cf.cl.int<-max(0,Cf.cl.int-0.1)
						Cf_cl_int = max_s_d(0.0, Cf_cl_int - 0.1);
//               if (Cf.cl>0)
						if (Cf_cl > 0) {
//                  if (Cf.cl.int!=Mcc)
							if (!Uguale(Cf_cl_int, Mcc)) {
//                     g_segno<-sign(g_CCs-Mcc)
								g_tmp1_i = somma_vs_i(g_tmp1_i, g_CCs, -Mcc);
								g_segno = segno_v_i(g_segno, g_tmp1_i);
//                     g_indp<-which(g_segno>0)
								g_indp = which_v_indxlt_i(g_indp, g_segno, 0) ;
//                     g_indm<-which(g_segno<=0)
								g_indm = which_v_indxle_i(g_indp, g_segno, 0) ;
//                     if (Cf.cl.int>Mcc)
								if (Cf_cl_int > Mcc) {
//                       g_p[g_indp]<-(1-abs(Cf.cl.int-g_CCs[g_indp]))+g_p[g_indp]
									g_p = f_aux9_d(g_p, 1, g_CCs, g_indp, Cf_cl_int);
//                        g_p[g_indm]<-(-abs(Cf.cl.int-g_CCs[g_indm]))+g_p[g_indm]
									g_p = f_aux9_d(g_p, 1, g_CCs, g_indm, Cf_cl_int);
//                       }
								}
//                     else
								else {
//                       g_p[g_indm]<-(1-abs(Cf.cl.int-g_CCs[g_indm]))+g_p[g_indm]
									g_p = f_aux9_d(g_p, 1, g_CCs, g_indm, Cf_cl_int);
//                        g_p[g_indp]<-(-abs(Cf.cl.int-g_CCs[g_indp]))+g_p[g_indp]
									g_p = f_aux9_d(g_p, 1, g_CCs, g_indp, Cf_cl_int);
//                       }
								}
//                     if (min(g_p)<=0) g_p[which(g_p<0)]<-0.0001
								if (min_v_d(g_p) <= 0.0)
									assegna1_v_indxle_d(g_p, g_p, 0.0, 0.0001);
//                     g_p<-g_p/sum(g_p)
								dividi1_vs_d(g_p, somma_v_d(g_p, false));
//                     for (ccc in (1:Lc)) g_CCp[(g_CCind[[ccc]])]<-g_p[ccc]
								for (ccc = 1; ccc <= g_LC; ccc++)
									assegna1_v_indx_d(g_CCp, g_CCind[ccc - 1], ACCEDIv_d(g_p, ccc));
//                     }
							}
//                   else
							else {
//                   {g_p<-g_p+(1-abs(Cf.cl.int-g_CCs))
								g_tmp1_i = seq_i(g_tmp1_i, 1, LENGTHv_d(g_p), 1);
								g_p = f_aux9_d(g_p, 1, g_CCs, g_tmp1_i, Cf_cl_int);
//                   g_p<-g_p/sum(g_p)
								dividi1_vs_d(g_p, somma_v_d(g_p, false));
//                   for (ccc in (1:Lc)) g_CCp[(g_CCind[[ccc]])]<-g_p[ccc]
								for (ccc = 1; ccc <= g_LC; ccc++)
									assegna1_v_indx_d(g_CCp, g_CCind[ccc - 1], ACCEDIv_d(g_p, ccc));
//                   }
							}
//                 }
						}
//               g_ppp<-prior.g_p.subnet*CCp
						g_ppp = moltiplica_vv_d(g_ppp, g_prior_p_subnet, g_CCp);
//               if ((max(g_ppp)==0)|(is.na(max(g_ppp))))   stop("\n prior.g_p.subnet values incompatible with other parameter settings! \n")
						mx = max_v_d(g_ppp);
						if (Uguale(mx, 0) || ISNA(mx))
							error("g_prior_p_subnet values incompatible with other parameter settings! \n");
//               g_ppp<-g_ppp/sum(g_ppp,na.rm=TRUE)
						dividi1_vs_d(g_ppp, somma_v_d(g_ppp, true));
//              }
					}

//            else
					else {
//             if ((gamma==0)|(!is.null(DEGREE))) stop("computation failed! Try again or consider using different tolerance parameter settings or different topology parameter settings")
						if (Uguale(gamma, 0.0) || degree != NULL)
							error("computation failed! Try again or consider using different tolerance parameter settings or different topology parameter settings");
// 		         gamma<-max(gamma-0.2,0)
						gamma = max_s_d(gamma - 0.2, 0.0);
//              if (gamma!=0)
						if (!Uguale(gamma, 0.0)) {
// 		          g_Prob<-c(seq(1,N,1)^(-gamma),0)
							g_tmp1_d = seq_d(g_tmp1_d, 1.0, (double) N, 1.0);
							g_tmp2_d = exp_d(g_tmp2_d, g_tmp1_d, -gamma);
							ASSEGNAv_d(g_scalare_d, 1, 0.0);
							g_Prob = vettore2v_d(g_Prob, g_tmp2_d, g_scalare_d);
//  		           g_Prob<-g_Prob/(sum(g_Prob))
							dividi1_vs_d(g_Prob, somma_v_d(g_Prob, false));
//  		           g_STout<-c(NA,g_Prob*N)         # DEGREE==NULL
							g_tmp1_d = moltiplica_vs_d(g_tmp1_d, g_Prob, (double) N);
							ASSEGNAv_d(g_scalare_d, 1, NA_REAL);
							g_STout = vettore2v_d(g_STout, g_scalare_d, g_tmp1_d);
//                g_STin<-rep(0,N+2)
							g_STin = rep_s_d(g_STin, 0.0, N + 2);
//                if (INdegree=="out")
							if (!strncmp(INdegree->str, OUT, 3)) {
//                 if (zero.nodes.with.indegree0==TRUE)   g_STin[2:(max.con+1)]<-N*g_Prob[2:(max.con+1)]/sum(g_Prob[2:(max.con+1)],na.rm=TRUE)   #sicuramente siamo nel caso (INdegree=="out")
								if (zero_nodes_with_indegree0) {
									g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 2, max_con + 1);
									g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
									// in somma cancNA = 1
									dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, true));
									assegna1_v_segmv_d(g_STin, 2, max_con + 1, g_tmp2_d);
								}
//                  else  g_STin[1:(max.con+1)]<-N*g_Prob[1:(max.con+1)]/sum(g_Prob[1:(max.con+1)],na.rm=TRUE)
								else {
									g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 1, max_con + 1);
									g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
									// in somma cancNA = 1
									dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, true));
									assegna1_v_segmv_d(g_STin, 1, max_con + 1, g_tmp2_d);
								}
//                 }
							}
//                else   g_STin[1:(max.con+1)]<-NA
							else
								assegna1_v_segm_d(g_STin, 1, max_con + 1, NA_REAL);
//                if (zero.nodes.with.indegree0==TRUE)    g_STin[1]<-0
							if (zero_nodes_with_indegree0 == TRUE)
								ASSEGNAv_d(g_STin, 1, 0.0);
//                Cf.cl.int<-Cf.cl
							Cf_cl_int = Cf_cl;
// 		           if (Cf.cl>0)
							if (Cf_cl > 0) {
//                  if (Cf.cl.int!=Mcc)
								if (!Uguale(Cf_cl_int, Mcc)) {
//                     g_segno<-sign(g_CCs-Mcc)
									g_tmp1_i = somma_vs_i(g_tmp1_i, g_CCs, -Mcc);
									g_segno = segno_v_i(g_segno, g_tmp1_i);
//                     g_indp<-which(g_segno>0)
									g_indp = which_v_indxgt_i(g_indp, g_segno, 0) ;
//                     g_indm<-which(g_segno<=0)
									g_indm = which_v_indxle_i(g_indm, g_segno, 0);
//                     if (Cf.cl.int>Mcc)
									if (Cf_cl_int > Mcc) {
//                       g_p[g_indp]<-(1-abs(Cf.cl.int-g_CCs[g_indp]))+g_p[g_indp]
										g_p = f_aux9_d(g_p, 1, g_CCs, g_indp, Cf_cl_int);
//                        g_p[g_indm]<-(-abs(Cf.cl.int-g_CCs[g_indm]))+g_p[g_indm]
										g_p = f_aux9_d(g_p, 0, g_CCs, g_indm, Cf_cl_int);
//                       }
									}
//                     else
									else {
//                      g_p[g_indm]<-(1-abs(Cf.cl.int-g_CCs[g_indm]))+g_p[g_indm]
										g_p = f_aux9_d(g_p, 1, g_CCs, g_indm, Cf_cl_int);
//                        g_p[g_indp]<-(-abs(Cf.cl.int-g_CCs[g_indp]))+g_p[g_indp]
										g_p = f_aux9_d(g_p, 0, g_CCs, g_indp, Cf_cl_int);
//                       }
									}
//                     if (min(g_p)<=0) g_p[which(g_p<0)]<-0.0001
									if (min_v_d(g_p) <= 0)
										assegna1_v_indxle_d(g_p, g_p, 0.0, 0.0001);
//                     g_p<-g_p/sum(g_p)
									dividi1_vs_d(g_p, somma_v_d(g_p, false));
//                     for (ccc in (1:Lc)) g_CCp[(g_CCind[[ccc]])]<-g_p[ccc]
									for (ccc = 1; ccc <= g_LC; ccc++)
										assegna1_v_indx_d(g_CCp, g_CCind[ccc - 1], ACCEDIv_d(g_p, ccc));
//                     }
								}
//                   else
								else {
//                   g_p<-g_p+(1-abs(Cf.cl.int-g_CCs))
									g_tmp1_i = seq_i(g_tmp1_i, 1, LENGTHv_d(g_p), 1);
									g_p = f_aux9_d(g_p, 1, g_CCs, g_tmp1_i, Cf_cl_int);
									dividi1_vs_d(g_p, somma_v_d(g_p, false));
//                   for (ccc in (1:Lc)) g_CCp[(g_CCind[[ccc]])]<-g_p[ccc]
									for (ccc = 1; ccc <= g_LC; ccc++)
										assegna1_v_indx_d(g_CCp, g_CCind[ccc - 1], ACCEDIv_d(g_p, ccc));
//                   }
								}
//                 }
							}
//                 g_ppp<-prior.g_p.subnet*CCp
							g_ppp = moltiplica_vv_d(g_ppp, g_prior_p_subnet, g_CCp);
//                 if ((max(g_ppp)==0)|(is.na(max(g_ppp))))   stop("\n prior.g_p.subnet values incompatible with other parameter settings! \n")
							mx = max_v_d(g_ppp);
							if (Uguale(mx, 0.0) || ISNA(mx))
								error("g_prior_p_subnet values incompatible with other parameter settings! \n");
//                 g_ppp<-g_ppp/sum(g_ppp,na.rm=TRUE)
							dividi1_vs_d(g_ppp, somma_v_d(g_ppp, true));
//                 cat("\n WARNING: gamma set to",gamma,"in iteration",it+1,"connecting the remaining",length(g_h)+length(g_h.new),"g_hubs, because higher gamma is not compatible with current network structure and other parameters setting\n")
							g_string_sprintf(tmp_s, "gamma set to %.7g in iteration %d connecting the remaining %d g_hubs, because higher gamma is not compatible with current network structure and other parameters setting\n", gamma, it, LENGTHv_i(g_h) + LENGTHv_i(g_h_new));
							warning(tmp_s->str);
// 		           }
						}
// 		         else
						else {
// 		         g_Prob<-c(rep(1/N,N),0)
							g_Prob = rep_s_d(g_Prob, (double) 1.0 / N, N);
							g_Prob = accoda1_vs_d(g_Prob, 0.0);
//  		           g_STout<-c(NA,g_Prob*N)              # DEGREE==NULL
							g_tmp1_d = moltiplica_vs_d(g_tmp1_d, g_Prob, (double) N);
							ASSEGNAv_d(g_scalare_d, 1, NA_REAL);
							g_STout = vettore2v_d(g_STout, g_scalare_d, g_tmp1_d);
// 		           g_STin<-rep(0,N+2)
							g_STin = rep_s_d(g_STin, 0.0, N + 2);
//                if (INdegree=="out")
							if (!strncmp(INdegree->str, OUT, 3)) {
//                 if (zero.nodes.with.indegree0==TRUE)   g_STin[2:(max.con+1)]<-N*g_Prob[2:(max.con+1)]/sum(g_Prob[2:(max.con+1)],na.rm=TRUE)   #sicuramente siamo nel caso (INdegree=="out")
								if (zero_nodes_with_indegree0) {
									g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 2, max_con + 1);
									g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
									// in somma cancNA = 1
									dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, true));
									assegna1_v_segmv_d(g_STin, 2, max_con + 1, g_tmp2_d);
								}
//                  else  g_STin[1:(max.con+1)]<-N*g_Prob[1:(max.con+1)]/sum(g_Prob[1:(max.con+1)],na.rm=TRUE)
								else {
									g_tmp1_d = segmento_v_d(g_tmp1_d, g_Prob, 1, max_con + 1);
									g_tmp2_d = moltiplica_vs_d(g_tmp2_d, g_tmp1_d, (double) N);
									// in somma cancNA = 1
									dividi1_vs_d(g_tmp2_d, somma_v_d(g_tmp1_d, true));
									assegna1_v_segmv_d(g_STin, 1, max_con + 1, g_tmp2_d);
								}

//                 }
							}
//                else   g_STin[1:(max.con+1)]<-NA
							else
								assegna1_v_segm_d(g_STin, 1, max_con + 1, NA_REAL);
//                if (zero.nodes.with.indegree0==TRUE)    g_STin[1]<-0
							if (zero_nodes_with_indegree0)
								ASSEGNAv_d(g_STin, 1, 0.0);
//                Cf.cl.int<-Cf.cl
							Cf_cl_int = Cf_cl;
// 		           if (Cf.cl>0)
							if (Cf_cl > 0) {
//                  if (Cf.cl.int!=Mcc)
								if (!Uguale(Cf_cl_int, Mcc)) {
//                     g_segno<-sign(g_CCs-Mcc)
									g_tmp1_i = somma_vs_i(g_tmp1_i, g_CCs, -Mcc);
									g_segno = segno_v_i(g_segno, g_tmp1_i);
//                     g_indp<-which(g_segno>0)
									g_indp = which_v_indxgt_i(g_indp, g_segno, 0) ;
//                     g_indm<-which(g_segno<=0)
									g_indm = which_v_indxle_i(g_indm, g_segno, 0);
//                     if (Cf.cl.int>Mcc)
									if (Cf_cl_int > Mcc) {
//                       g_p[g_indp]<-(1-abs(Cf.cl.int-g_CCs[g_indp]))+g_p[g_indp]
										g_p = f_aux9_d(g_p, 1, g_CCs, g_indp, Cf_cl_int);
//                        g_p[g_indm]<-(-abs(Cf.cl.int-g_CCs[g_indm]))+g_p[g_indm] // ERRORE?
										g_p = f_aux9_d(g_p, 0, g_CCs, g_indm, Cf_cl_int);
//                       }
									}
//                     else
									else {
//                       g_p[g_indm]<-(1-abs(Cf.cl.int-g_CCs[g_indm]))+g_p[g_indm]
										g_p = f_aux9_d(g_p, 1, g_CCs, g_indm, Cf_cl_int);
//                        g_p[g_indp]<-(-abs(Cf.cl.int-g_CCs[g_indp]))+g_p[g_indp] // ERRORE?
										g_p = f_aux9_d(g_p, 0, g_CCs, g_indp, Cf_cl_int);
//                       }
									}
//                     if (min(g_p)<=0) g_p[which(g_p<0)]<-0.0001
									if (min_v_d(g_p) <= 0)
										assegna1_v_indxle_d(g_p, g_p, 0.0, 0.0001);
//                     g_p<-g_p/sum(g_p)
									dividi1_vs_d(g_p, somma_v_d(g_p, false));
//                     for (ccc in (1:Lc)) g_CCp[(g_CCind[[ccc]])]<-g_p[ccc]
									for (ccc = 1; ccc <= g_LC; ccc++)
										assegna1_v_indx_d(g_CCp, g_CCind[ccc - 1], ACCEDIv_d(g_p, ccc));
//                     }
								}
//                   else
								else {
//                   {g_p<-g_p+(1-abs(Cf.cl.int-g_CCs))
									g_tmp1_i = seq_i(g_tmp1_i, 1, LENGTHv_d(g_p), 1);
									g_p = f_aux9_d(g_p, 1, g_CCs, g_tmp1_i, Cf_cl_int);
//                   g_p<-g_p/sum(g_p)
									dividi1_vs_d(g_p, somma_v_d(g_p, false));
//                   for (ccc in (1:Lc)) g_CCp[(g_CCind[[ccc]])]<-g_p[ccc]
									for (ccc = 1; ccc <= g_LC; ccc++)
										assegna1_v_indx_d(g_CCp, g_CCind[ccc - 1], ACCEDIv_d(g_p, ccc));
//                   }
								}
//                 }
							}
//               g_ppp<-prior.g_p.subnet*CCp
							g_ppp = moltiplica_vv_d(g_ppp, g_prior_p_subnet, g_CCp);
//               if ((max(g_ppp)==0)|(is.na(max(g_ppp))))   stop("\n prior.g_p.subnet values incompatible with other parameter settings! \n")
							if (Uguale(max_v_d(g_ppp), 0.0) || !ISNA(max_v_d(g_ppp)))
								error("g_prior_p_subnet values incompatible with other parameter settings! \n");
//               g_ppp<-g_ppp/sum(g_ppp,na.rm=TRUE)
							dividi1_vs_d(g_ppp, somma_v_d(g_ppp, true));
//               cat("\n WARNING: power law distribution replaced by flat distribution in iteration",it+1,"connecting the remaining",length(g_h)+length(g_h.new),"g_hubs, because power-law distribution is not compatible with current network structure and other parameters setting\n")
							g_string_sprintf(tmp_s, "power law distribution replaced by flat distribution in iteration %d connecting the remaining %d g_hubs, because power-law distribution is not compatible with current network structure and other parameters setting\n", it + 1, LENGTHv_i(g_h) + LENGTHv_i(g_h_new));
							warning(tmp_s->str);
// 		          }
						}
//              }
					}
// 	         }#end if (max(g_ppp,na.rm=TRUE)==0)
				}//end if (max(g_ppp,na_rm=TRUE)==0)
//         }
			}

//      } #end controllo==1
		} //end controllo==1;

//    if (controllo2==1)
		if (controllo2) {
//   {Ng<-g_dim(g_M)[1]
			Ng = LENGTHm1_i(g_M);
//     if (length(g_hubs)==Ng) g_hubs<-sampleB(g_hubs,Ng%/%2)
			if (LENGTHv_i(g_hubs) == Ng)
				g_hubs = sampleB(g_hubs, g_hubs, floor(Ng / 2), false);
			// "g_h" verra` modificato, quindi faccio una copia
			g_h1 = copia_v_i(g_h1, g_h, 1, LENGTHv_i(g_h));
//     g_aus<-assign.nodes(g_M,Mdiscr,g_h,g_hubs=g_hubs,g_Sc=g_Sc,g_Sin,max.con)
			g_h_new1 = assign_nodes2(g_h_new1, g_M, ris, g_h1, g_hubs, g_Sc, g_Sin, max_con);
//     Mdiscr<-g_aus[[1]]
			// non serve, viene modificata per riferimento ,come pure la "g_h" seguente che corrisponde a g_aus[[2]]
//     ########livello gerarchico#############
			//#######livello gerarchico#############;
//     g_LG[setdiff(g_h,g_aus[[2]])]<-g_LG[setdiff(g_h,g_aus[[2]])]+1
			g_tmp1_i = setdiff_i(g_tmp1_i, g_h, g_h1);
			incr1_v_indx_i(g_LG, g_tmp1_i, 1);
//     #######################################
			//######################################;
//     g_h<-g_aus[[2]]
			g_h = copia_v_i(g_h, g_h1, 1, LENGTHv_i(g_h1));
//     g_h.new<-c(g_h.new,g_aus[[3]])
			g_h_new = accoda1_vv_i(g_h_new, g_h_new1);
//     g_Sout<-apply(Mdiscr,2,sum)
			g_Sout = somma_colonne_i(g_Sout, ris);
//     g_m<-max(g_Sout)+1
			m1 = max_v_i(g_Sout) + 1;
//     Freq.out[1:g_m]<-hist(g_Sout,breaks=seq(0,g_m,1),right=FALSE,plot=FALSE)$counts
			g_tmp2_i = seq_i(g_tmp2_i, 0, m1, 1);
			g_tmp1_i = hist1(g_tmp1_i, g_Sout, g_tmp2_i, 0, 1, 0);
			g_tmp1_d = promuovi_i(g_tmp1_d, g_tmp1_i);
			assegna1_v_segmv_d(g_Freq_out, 1, m1, g_tmp1_d);
//     g_Sin<-apply(Mdiscr,1,sum)
			g_Sin = somma_righe_i(g_Sin, ris);
//     g_m<-max(g_Sin)+1
			m1 = max_v_i(g_Sin) + 1;
//     Freq.in[1:g_m]<-hist(g_Sin,breaks=seq(0,g_m,1),right=FALSE,plot=FALSE)$counts
			g_tmp2_i = seq_i(g_tmp2_i, 0, m1, 1);
			g_tmp1_i = hist1(g_tmp1_i, g_Sin, g_tmp2_i, 0, 1, 0);
			g_tmp1_d = promuovi_i(g_tmp1_d, g_tmp1_i);
			assegna1_v_segmv_d(g_Freq_in, 1, m1, g_tmp1_d);
//     Nrim<-Nrim-Ng
			Nrim -= Ng;
//     if (Nrim<=1)
			if (Nrim <= 1) {
//      #C.current<-cluster.coeff(Mdiscr)[[1]]
				//C_current=cluster_coeff(Mdiscr)[[1]];
//       it<-it+1
				it++;
//       if ((it==1)&(Nrim==1)) g_h<-c(g_h,g_h.new)
				if ((it == 1) && (Nrim == 1))
					g_h = accoda1_vv_i(g_h, g_h_new);
//        else g_h<-g_h.new
				else
					g_h = copia_v_i(g_h, g_h_new, 1, LENGTHv_i(g_h_new));
//       g_h.new<-vector()
				CREAv_i(g_h_new, 0);
//       if (length(g_h)>0) g_h<-g_h[which(g_Sin[g_h]!=max.con)]
				if (LENGTHv_i(g_h) > 0) {
					g_tmp2_i = copia_v_indx_i(g_tmp2_i, g_Sin, g_h);
					g_tmp1_i = which_v_indxne_i(g_tmp1_i, g_tmp2_i, max_con);
					g_h = copia_v_indx_i(g_h, g_h, g_tmp1_i);
				}
//      }
			}
//    }
		}
//    else  #controllo2==0
		else  {//controllo2==0;
//     if (it>0)
			if (it > 0) {
//       if (length(g_h.new)==0) g_h<-vector()
				if (LENGTHv_i(g_h_new) == 0)
					CREAv_i(g_h, 0);
//        else
				else {
//          it<-it+1
					it++;
//           g_h<-g_h.new
					g_h = copia_v_i(g_h, g_h_new, 1, LENGTHv_i(g_h_new));
//           g_h.new<-vector()
					CREAv_i(g_h_new, 0);
//           if (length(g_h)>0) g_h<-g_h[which(g_Sin[g_h]!=max.con)]
					if (LENGTHv_i(g_h) > 0) {
						g_tmp2_i = copia_v_indx_i(g_tmp2_i, g_Sin, g_h);
						g_tmp1_i = which_v_indxne_i(g_tmp1_i, g_tmp2_i, max_con);
						g_h = copia_v_indx_i(g_h, g_h, g_tmp1_i);
					}
//          }
				}
//        }
			}
//       else #it==0
			else { //it==0;
//         if (length(g_h.new)==0) stop("\n prior.g_p.subnet values incompatible with other parameter settings! \n")
				if (LENGTHv_i(g_h_new) == 0)
					error("\n g_prior_p_subnet values incompatible with other parameter settings! \n");
//          else
				else {
//          {it<-it+1
					it++;
//           g_h<-c(g_h,g_h.new)
					g_h = copia_v_i(g_h, g_h_new, 1, LENGTHv_i(g_h_new));
//           g_h.new<-vector()
					CREAv_i(g_h_new, 0);
//           if (length(g_h)>0) g_h<-g_h[which(g_Sin[g_h]!=max.con)]
					if (LENGTHv_i(g_h) > 0) {
						g_tmp2_i = copia_v_indx_i(g_tmp2_i, g_Sin, g_h);
						g_tmp1_i = which_v_indxne_i(g_tmp1_i, g_tmp2_i, max_con);
						g_h = copia_v_indx_i(g_h, g_h, g_tmp1_i);
					}
//          }
				}
//        }
			}
//      }
		}
//    }# END WHILE  (length(g_h)>1)
	}// END WHILE  (length(g_h)>1);

// ###check graph connectivity and add links if it is not completely connected
//##check graph connectivity and add links if it is not completely connected;
// if (sepgraph==FALSE)
	if (!sepgraph) {
// g_dist<-check.conn(Mdiscr)
		g_dist = check_conn1(g_dist, ris);
//  non.connessi<-which(g_dist==Inf)
		g_non_connessi = which_v_indxeq_d(g_non_connessi, g_dist, R_PosInf);
//  L<-length(non.connessi)
		L = LENGTHv_i(g_non_connessi);
//  while (L>0)
		while (L > 0) {
//   Mdiscr<-connetti.scalefree(Mdiscr,g_STout,g_STin,g_dist,g_toll,max.con,und=FALSE)
			ris = connetti_scalefree(ris, g_STout, g_STin, g_dist, g_toll, max_con, false);
//    g_dist<-check.conn(Mdiscr)
			g_dist = check_conn1(g_dist, ris);
//    non.connessi<-which(g_dist==Inf)
			g_non_connessi = which_v_indxeq_d(g_non_connessi, g_dist, R_PosInf);
//    L<-length(non.connessi)
			L = LENGTHv_i(g_non_connessi);
//   }
		}
// }
	}
// if (autoreg==FALSE) diag(Mdiscr)<-0
	if (!autoreg)
		assegna1_s_diag_i(ris, 0);
// if (zero.nodes.with.indegree0==TRUE)
	if (zero_nodes_with_indegree0) {
// {###check that every gene has at least 1 regulator
		//##check that every gene has at least 1 regulator ;
//  g_Sr<-apply(Mdiscr,1,sum)
		g_Sr = somma_righe_i(g_Sr, ris);
//  g_ind<-which(g_Sr==0)
		g_ind = which_v_indxeq_i(g_ind, g_Sr, 0);
//  L<-length(g_ind)
		L = LENGTHv_i(g_ind);
//  if (L>0)
		if (L > 0) {
//   for (i in (1:L))
			for (i = 1; i <= L; i++) {
//     ri<-g_ind[i]
				ri = ACCEDIv_i(g_ind, i);
//      num<-1
				num = 1;
//      g_Sout<-apply(Mdiscr,2,sum)
				g_Sout = somma_colonne_i(g_Sout, ris);
//      g_Sc<-Score(S=g_Sout,ST=g_STout,Freq=Freq.out,n=1,g_toll)
				g_Sc_vett = score1(g_Sc_vett, g_Sout, g_STout, g_Freq_out, 1, g_toll);
//      g_ind.g_Sc<-which(g_Sc==(-Inf))
				g_ind_Sc = which_v_indxeq_d(g_ind_Sc, g_Sc_vett, R_NegInf);
//      g_Sc[g_ind.g_Sc]<-0
				assegna1_v_indx_d(g_Sc_vett, g_ind_Sc, 0.0);
//      if (length(which(g_Sc<=0))>0) g_Sc<-g_Sc-min(g_Sc)+1/(N^2)
				g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc_vett, 0.0);
				if (LENGTHv_i(g_tmp1_i) > 0)
					somma1_vs_d(g_Sc_vett, -(min_v_d(g_Sc_vett) + 1 / (N * N)));
//      g_Sc[g_ind.g_Sc]<-0
				assegna1_v_indx_d(g_Sc_vett, g_ind_Sc, 0.0);
//      if (autoreg==FALSE) g_Sc[ri]<-0
				if (!autoreg)
					ASSEGNAv_d(g_Sc_vett, ri, 0.0);
//      if (sum(g_Sc)==0)
				if (Uguale(somma_v_d(g_Sc_vett, false), 0.0)) {
// 	    {cat("\n WARNING: tolerance parameters were relaxed to allow each node to have at least 1 regulator")
					warning("tolerance parameters were relaxed to allow each node to have at least 1 regulator");
					g_tmp1_d = rep_s_d(g_tmp1_d, R_PosInf, N);
//        g_Sc<-Score(S=g_Sout,ST=g_STout,Freq=Freq.out,n=1,g_toll=rep(Inf,N))
					g_Sc_vett = score1(g_Sc_vett, g_Sout, g_STout, g_Freq_out, 1, g_tmp1_d);
//        g_ind.g_Sc<-which(g_Sc==(-Inf))
					g_ind_Sc = which_v_indxeq_d(g_ind_Sc, g_Sc_vett, R_NegInf);
// 	     g_Sc[g_ind.g_Sc]<-0
					assegna1_v_indx_d(g_Sc_vett, g_ind_Sc, 0.0);
//  	     if (length(which(g_Sc<=0))>0) g_Sc<-g_Sc-min(g_Sc)+1/(N^2)
					g_tmp1_i = which_v_indxle_d(g_tmp1_i, g_Sc_vett, 0.0);
					if (LENGTHv_i(g_tmp1_i) > 0)
						somma1_vs_d(g_Sc_vett, -(min_v_d(g_Sc_vett) + 1 / (N * N)));
//  	     g_Sc[g_ind.g_Sc]<-0
					assegna1_v_indx_d(g_Sc_vett, g_ind_Sc, 0.0);
//  	     if (autoreg==FALSE) g_Sc[ri]<-0
					if (!autoreg)
						ASSEGNAv_d(g_Sc_vett, ri, 0.0);
// 	    }
				}
//      g_p.out<-g_Sc/sum(g_Sc)
				g_p_out = dividi_vs_d(g_p_out, g_Sc_vett, somma_v_d(g_Sc_vett, false));
//      g_ind.s<-sampleB(seq(1,N,1),num,prob=g_p.out)
				g_tmp1_i = seq_i(g_tmp1_i, 1, N, 1);
				g_ind_s = sampleB_p(g_ind_s, g_tmp1_i, num, false, g_p_out);
//      Mdiscr[ri,g_ind.s]<-1
				assegna1_ms_rigaindx_i(ris, ri, g_ind_s, 1);
//      g_a1<-g_Sout[g_ind.s]
				g_a1 = copia_v_indx_i(g_a1, g_Sout, g_ind_s);
//      g_a2<-g_a1+1
				g_a2 = somma_vs_i(g_a2, g_a1, 1);
//      Freq.out[g_a1+1]<-Freq.out[g_a1+1]-1
				incr1_v_i(g_a1, 1);
				incr1_v_indx_d(g_Freq_out, g_a1, -1.0);
//      Freq.out[g_a2+1]<-Freq.out[g_a2+1]+1
				incr1_v_i(g_a2, 1);
				incr1_v_indx_d(g_Freq_out, g_a2, -1.0);
//     }
			}
//   }
		}
// }
	}
// etich<-paste("topology_",iter,".txt",sep="")
	g_string_sprintf(etich, "topology_%d.txt", iter);
// write.table(Mdiscr,etich,sep="\t",row.names=FALSE,col.names=FALSE)
	write_m_i(etich->str, ris);
// etich<-paste("Hierarch_level",iter,".txt",sep="")
	g_string_sprintf(etich, "Hierarch_level%d.txt", iter);
// write.table(cbind(1:N,g_LG),etich,sep="\t",row.names=FALSE,col.names=FALSE)
	g_tmp1_i = seq_i(g_tmp1_i, 1, N, 1);
	g_tmpm_i = cbind2v_i(g_tmpm_i, g_tmp1_i, g_LG);
	write_m_i(etich->str, g_tmpm_i);

	StrBilanciam();

	_Intestazione("\n*** Esco da HMM1 ***\n");

// return(Mdiscr)
	return ris;
// }
}

// HMM<-function(N=50, Cf.cl=0.3, gamma=2.2, DEGREE=NULL, INdegree=c("free", "out"), MODULES, prior.g_p.subnet=NULL, max.con=12, feedback=TRUE, zero.nodes.with.indegree0=TRUE, sepgraph=TRUE, r.tol=0.1, a.tol=1, iter=1)
SEXP HMM(SEXP N, SEXP Cf_cl, SEXP gamma, SEXP degree, SEXP INdegree, SEXP modules, SEXP prior_p_subnet, SEXP max_con, SEXP feedback, SEXP zero_nodes_with_indegree0, SEXP sepgraph, SEXP r_tol, SEXP a_tol, SEXP iter) {
	int i, nProtected = 0, l, N1, max_con1, iter1;
	double Cf_cl1, gamma1, r_tol1, a_tol1;
	bool feedback1, zero_nodes_with_indegree01, sepgraph1;
	VETTOREd *prior_p_subnet1, *degree1;
	MATRICEi *ris1 = NULL;
	GString *INdegree1;
	Mod *modules1;
	SEXP elem, ris;

	_InitDbg(false, false, false);

	_Intestazione("\n*** HMM ***\n");

	N1 = INTEGER_VALUE(N);
	Cf_cl1 = NUMERIC_VALUE(Cf_cl);
	gamma1 = NUMERIC_VALUE(gamma);
	degree1 = inVETTORE_d(degree, &nProtected);
	INdegree1 = inSTRINGA(INdegree, &nProtected, "INdegree");
	prior_p_subnet1 = inVETTORE_d(prior_p_subnet, &nProtected);
	max_con1 = INTEGER_VALUE(max_con);
	feedback1 = LOGICAL_VALUE(feedback);
	zero_nodes_with_indegree01 = LOGICAL_VALUE(zero_nodes_with_indegree0);
	sepgraph1 = LOGICAL_VALUE(sepgraph);
	r_tol1 = NUMERIC_VALUE(r_tol);
	a_tol1 = NUMERIC_VALUE(a_tol);
	iter1 = INTEGER_VALUE(iter);

	l = length(modules);
	modules1 = mia_alloc(l, Mod);
	if (l > 0 && modules1 == NULL) {
		Rprintf("Not enough memory (HMM # %d, modules1)", __LINE__ - 2);
		error("");
	}
	for (i = 0; i < l; i++) {
		elem = VECTOR_ELT(modules, i);
		modules1[i].codice = INTEGER_VALUE(VECTOR_ELT(elem, 0));
		modules1[i].rete = inMATRICE_i(VECTOR_ELT(elem, 1), &nProtected);
		modules1[i].hubs = inVETTORE_i(VECTOR_ELT(elem, 2), &nProtected);
		modules1[i].CC = INTEGER_VALUE(VECTOR_ELT(elem, 3));
		modules1[i].autoreg = LOGICAL_VALUE(VECTOR_ELT(elem, 4));
		modules1[i].feedback = LOGICAL_VALUE(VECTOR_ELT(elem, 5));
		modules1[i].hubsio = inVETTORE_i(VECTOR_ELT(elem, 6), &nProtected);
		modules1[i].simm = LOGICAL_VALUE(VECTOR_ELT(elem, 7));
		modules1[i].dim_m = INTEGER_VALUE(VECTOR_ELT(elem, 8));
	}

	InitGlobali();

	Rprintf("SBNT: ");
	for (i = 1; i <= iter1; i++) {
		ris1 = HMM1(ris1, N1, Cf_cl1, gamma1, degree1, INdegree1, modules1, l, prior_p_subnet1, max_con1, feedback1, zero_nodes_with_indegree01, sepgraph1, r_tol1, a_tol1, i);
		Rprintf(".");
	}
	Rprintf(" done\n");

	ris = daMATRICE_i(ris1, &nProtected);

	for (i = 0; i < l; i++) {
		CANCELLAm_i(modules1[i].rete);
		CANCELLAv_i(modules1[i].hubs);
		CANCELLAv_i(modules1[i].hubsio);
	}
	CANCELLAv_d(degree1);
	CANCELLAv_d(prior_p_subnet1);
	CANCELLAstr(INdegree1);

	CancGlobali();

	StrBilanciam();
	ControllaCanc();

	UNPROTECT(nProtected);

	libera(modules1);

	return ris;
}
