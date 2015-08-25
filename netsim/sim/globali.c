#include "globali.h"

void InitGlobali()
{
	globali.n_tipi = 0;
	globali.tipi = NULL;

	globali.sample.HL = NULL;
	globali.sample.q = NULL;
	globali.sample.x = NULL;
	globali.sample.y = NULL;

	globali.ret_vettore.v2 = NULL;
	globali.ret_vettore.v3 = NULL;

	globali.ret_lista.v2 = NULL;
	globali.ret_lista.m1 = NULL;

	globali.check_conn.Mt = NULL;
	globali.check_conn.Maus = NULL;
	globali.check_conn.tmp_i = NULL;
	globali.check_conn.colore = NULL;
	globali.check_conn.grigi = NULL;
	globali.check_conn.ind = NULL;
	globali.check_conn.adj = NULL;
	globali.check_conn.tmp1_d = NULL;
	globali.check_conn.scalare_d = NULL;

	globali.boole_result.valr = NULL;
	globali.boole_result.signr = NULL;
	globali.boole_result.tmp_r = NULL;

	globali.triangola.tmp1_i = NULL;
	globali.triangola.tmp2_i = NULL;
	globali.triangola.tmp3_d = NULL;
	globali.triangola.M_in = NULL;
	globali.triangola.sk = NULL;
	globali.triangola.coord_ind = NULL;
	globali.triangola.ind_aus = NULL;
	globali.triangola.ind = NULL;
	globali.triangola.coord = NULL;
	globali.triangola.tmp_coord = NULL;
	globali.triangola.tmp_coord1 = NULL;
	globali.triangola.Dmem = NULL;

	globali.cluster_coeff.tmpm_i = NULL;
	globali.cluster_coeff.ind = NULL;

	globali.assign_nodes.Sin_h = NULL;
	globali.assign_nodes.tmp_i = NULL;
	globali.assign_nodes.tmp_d = NULL;
	globali.assign_nodes.or_h = NULL;
	globali.assign_nodes.aus_h = NULL;
	globali.assign_nodes.M_in = NULL;
	globali.assign_nodes.Ord = NULL;
	globali.assign_nodes.p = NULL;
	globali.assign_nodes.ind_h = NULL;
	globali.assign_nodes.ri = NULL;
	globali.assign_nodes.co = NULL;
	globali.assign_nodes.ind = NULL;
	globali.assign_nodes.index = NULL;

	globali.target.aus = NULL;
	globali.target.tmp_d = NULL;
	globali.target.tmp1_d = NULL;
	globali.target.ind = NULL;
	globali.target.ind_reg = NULL;

	globali.scoremodular.tmpm_d = NULL;
	globali.scoremodular.tmp_i1 = NULL;
	globali.scoremodular.tmp_d1 = NULL;
	globali.scoremodular.T1 = NULL;
	globali.scoremodular.T2 = NULL;
	globali.scoremodular.old1 = NULL;
	globali.scoremodular.old2 = NULL;
	globali.scoremodular.new1 = NULL;
	globali.scoremodular.new2 = NULL;
	globali.scoremodular.toll1 = NULL;
	globali.scoremodular.toll2 = NULL;
	globali.scoremodular.a = NULL;
	globali.scoremodular.b = NULL;
	globali.scoremodular.ind = NULL;
	globali.scoremodular.m = NULL;
	globali.scoremodular.S1 = NULL;
	globali.scoremodular.S2 = NULL;
	globali.scoremodular.ind1 = NULL;
	globali.scoremodular.ind2 = NULL;

	globali.score_sf.tmpm_d = NULL;
	globali.score_sf.scalare_d = NULL;
	globali.score_sf.tmp_i1 = NULL;
	globali.score_sf.tmp_d1 = NULL;
	globali.score_sf.T1 = NULL;
	globali.score_sf.T2 = NULL;
	globali.score_sf.old1 = NULL;
	globali.score_sf.old2 = NULL;
	globali.score_sf.new1 = NULL;
	globali.score_sf.new2 = NULL;
	globali.score_sf.toll1 = NULL;
	globali.score_sf.toll2 = NULL;
	globali.score_sf.a = NULL;
	globali.score_sf.b = NULL;
	globali.score_sf.ind = NULL;
	globali.score_sf.m = NULL;
	globali.score_sf.S1 = NULL;
	globali.score_sf.S2 = NULL;
	globali.score_sf.ind1 = NULL;
	globali.score_sf.ind2 = NULL;
	globali.score_sf.indinf = NULL;
	globali.score_sf.tmp_i2 = NULL;
	globali.score_sf.indbad = NULL;
	globali.score_sf.ind0 = NULL;
	globali.score_sf.ind3 = NULL;

	globali.probmod.S_out = NULL;
	globali.probmod.S_in = NULL;
	globali.probmod.Sc = NULL;
	globali.probmod.checkIN = NULL;
	globali.probmod.checkOUT = NULL;
	globali.probmod.memory = NULL;
	globali.probmod.M_in = NULL;
	globali.probmod.M_out = NULL;
	globali.probmod.tmp_i1 = NULL;
	globali.probmod.tmp_i2 = NULL;
	globali.probmod.tmp_i3 = NULL;
	globali.probmod.indInf = NULL;
	globali.probmod.I = NULL;
	globali.probmod.ord_ind = NULL;
	globali.probmod.rs = NULL;
	globali.probmod.ind1 = NULL;
	globali.probmod.scalare_i = NULL;
	globali.probmod.I_add = NULL;
	globali.probmod.tmp_d1 = NULL;
	globali.probmod.score_matr1 = NULL;
	globali.probmod.score_matr2 = NULL;
	globali.probmod.score_matr3 = NULL;

	globali.dinamica.D = NULL;
	globali.dinamica.ind = NULL;
	globali.dinamica.y_prec = NULL;
	globali.dinamica.targetT = NULL;
	globali.dinamica.targ = NULL;
	globali.dinamica.incr = NULL;
	globali.dinamica.n = NULL;
	globali.dinamica.aus = NULL;
	globali.dinamica.tmp_ris = NULL;

	globali.createNEG.tmp_i = NULL;
	globali.createNEG.segno = NULL;
	globali.createNEG.ind = NULL;

	globali.lsoda.y = NULL;
	globali.lsoda.y1 = NULL;
	globali.lsoda.y_prec = NULL;
	globali.lsoda.targetT = NULL;
	globali.lsoda.targ = NULL;
	globali.lsoda.tmp_ris = NULL;

	globali.createRules.op = NULL;
	globali.createRules.tmp = NULL;

	globali.create_logicRule.x = NULL;
	globali.create_logicRule.s = NULL;
	globali.create_logicRule.scalare_d = NULL;
	globali.create_logicRule.nvect = NULL;
	globali.create_logicRule.e = NULL;
	globali.create_logicRule.prob = NULL;
	globali.create_logicRule.tmp_i = NULL;
	globali.create_logicRule.tmp1_i = NULL;
	globali.create_logicRule.pr_and = NULL;
	globali.create_logicRule.pr_or = NULL;
	globali.create_logicRule.scalare_i = NULL;
	globali.create_logicRule.blacklist = NULL;
	globali.create_logicRule.black_p = NULL;
	globali.create_logicRule.o = NULL;

	globali.mod1.ind = NULL;
	globali.mod1.tmp1_i = NULL;
	globali.mod1.x = NULL;
	globali.mod1.s = NULL;

	globali.mod3.M_out = NULL;
	globali.mod3.indok = NULL;
	globali.mod3.tmp3_i = NULL;
	globali.mod3.tmp1_d = NULL;
	globali.mod3.tmpSTin = NULL;
	globali.mod3.tmpSTout = NULL;
	globali.mod3.Freq_out = NULL;
	globali.mod3.tmp3_d = NULL;
	globali.mod3.p = NULL;
	globali.mod3.Sc = NULL;
	globali.mod3.tmp2_i = NULL;
	globali.mod3.ind_M = NULL;
	globali.mod3.ind1 = NULL;
	globali.mod3.indS = NULL;
	globali.mod3.indBS = NULL;
	globali.mod3.Sin = NULL;
	globali.mod3.ind = NULL;
	globali.mod3.indInf = NULL;
	globali.mod3.tmp1_i = NULL;
	globali.mod3.scalare_i = NULL;
	globali.mod3.scalare_d = NULL;
	globali.mod3.p_sc = NULL;

	globali.module1.scalare_i = NULL;
	globali.module1.Ng = NULL;
	globali.module1.tmp1_i = NULL;
	globali.module1.conn_matr = NULL;
	globali.module1.indices = NULL;

	globali.module2.tmp1_i = NULL;
	globali.module2.Ng = NULL;
	globali.module2.conn_matr = NULL;
	globali.module2.indices = NULL;

	globali.module3.tmp1_i = NULL;
	globali.module3.Ng = NULL;
	globali.module3.Ng_UP = NULL;
	globali.module3.conn_matr = NULL;
	globali.module3.indices = NULL;

	globali.connectivity_geometric.M = NULL;
	globali.connectivity_geometric.Mdiscr = NULL;
	globali.connectivity_geometric.x = NULL;
	globali.connectivity_geometric.y = NULL;
	globali.connectivity_geometric.xy = NULL;
	globali.connectivity_geometric.d = NULL;
	globali.connectivity_geometric.s = NULL;
	globali.connectivity_geometric.regulatedind = NULL;
	globali.connectivity_geometric.indL = NULL;
	globali.connectivity_geometric.Sr = NULL;
	globali.connectivity_geometric.aus = NULL;
	globali.connectivity_geometric.ind = NULL;
	globali.connectivity_geometric.ind1 = NULL;
	globali.connectivity_geometric.ind0 = NULL;
	globali.connectivity_geometric.tmp1_i = NULL;
	globali.connectivity_geometric.tmp2_i = NULL;
	globali.connectivity_geometric.tmp1_d = NULL;

	globali.connectivity_scalefree.M = NULL;
	globali.connectivity_scalefree.Mdiscr = NULL;
	globali.connectivity_scalefree.x = NULL;
	globali.connectivity_scalefree.y = NULL;
	globali.connectivity_scalefree.d = NULL;
	globali.connectivity_scalefree.s = NULL;
	globali.connectivity_scalefree.o = NULL;
	globali.connectivity_scalefree.regulatedind = NULL;
	globali.connectivity_scalefree.indL = NULL;
	globali.connectivity_scalefree.Sr = NULL;
	globali.connectivity_scalefree.inthenet = NULL;
	globali.connectivity_scalefree.not_inthenet = NULL;
	globali.connectivity_scalefree.not_regulated = NULL;
	globali.connectivity_scalefree.numposs = NULL;
	globali.connectivity_scalefree.give_outlink = NULL;
	globali.connectivity_scalefree.aus_give_outlink = NULL;
	globali.connectivity_scalefree.indInf = NULL;
	globali.connectivity_scalefree.ind_Sc = NULL;
	globali.connectivity_scalefree.a1 = NULL;
	globali.connectivity_scalefree.a2 = NULL;
	globali.connectivity_scalefree.primi = NULL;
	globali.connectivity_scalefree.indici = NULL;
	globali.connectivity_scalefree.num_v = NULL;
	globali.connectivity_scalefree.mem_o = NULL;
	globali.connectivity_scalefree.available = NULL;
	globali.connectivity_scalefree.campione = NULL;
	globali.connectivity_scalefree.linked = NULL;
	globali.connectivity_scalefree.ind_s = NULL;
	globali.connectivity_scalefree.Sout = NULL;
	globali.connectivity_scalefree.Sin = NULL;
	globali.connectivity_scalefree.ind = NULL;
	globali.connectivity_scalefree.indok = NULL;
	globali.connectivity_scalefree.ind1 = NULL;
	globali.connectivity_scalefree.ind0 = NULL;
	globali.connectivity_scalefree.tmp1_i = NULL;
	globali.connectivity_scalefree.tmp2_i = NULL;
	globali.connectivity_scalefree.tmp1_d = NULL;
	globali.connectivity_scalefree.tmp2_d = NULL;
	globali.connectivity_scalefree.scalare_i = NULL;
	globali.connectivity_scalefree.scalare_d = NULL;
	globali.connectivity_scalefree.Prob = NULL;
	globali.connectivity_scalefree.Freq_in = NULL;
	globali.connectivity_scalefree.Freq_out = NULL;
	globali.connectivity_scalefree.STin = NULL;
	globali.connectivity_scalefree.STout = NULL;
	globali.connectivity_scalefree.p = NULL;
	globali.connectivity_scalefree.toll1 = NULL;
	globali.connectivity_scalefree.toll_in = NULL;
	globali.connectivity_scalefree.toll_out = NULL;
	globali.connectivity_scalefree.Sc = NULL;
	globali.connectivity_scalefree.p_ind = NULL;
	globali.connectivity_scalefree.p_out = NULL;
	globali.connectivity_scalefree.aus = NULL;

	globali.connectivity_random.M = NULL;
	globali.connectivity_random.Mdiscr = NULL;
	globali.connectivity_random.aus = NULL;
	globali.connectivity_random.tmp_d = NULL;
	globali.connectivity_random.scalare_i = NULL;
	globali.connectivity_random.Pnum = NULL;
	globali.connectivity_random.ind = NULL;
	globali.connectivity_random.num = NULL;
	globali.connectivity_random.tmp_i = NULL;

	globali.connectivity_modular.M = NULL;
	globali.connectivity_modular.Mdiscr = NULL;
	globali.connectivity_modular.scalare_i = NULL;
	globali.connectivity_modular.tmp1_i = NULL;
	globali.connectivity_modular.tmp2_i = NULL;
	globali.connectivity_modular.tmp1_d = NULL;
	globali.connectivity_modular.tmp2_d = NULL;
	globali.connectivity_modular.Prob = NULL;
	globali.connectivity_modular.Freq_out = NULL;
	globali.connectivity_modular.Freq_in = NULL;
	globali.connectivity_modular.STout = NULL;
	globali.connectivity_modular.STin = NULL;
	globali.connectivity_modular.toll = NULL;
	globali.connectivity_modular.p_out = NULL;
	globali.connectivity_modular.scalare_d = NULL;
	globali.connectivity_modular.Sr = NULL;
	globali.connectivity_modular.p = NULL;
	globali.connectivity_modular.h = NULL;
	globali.connectivity_modular.prob_mod = NULL;
	globali.connectivity_modular.aus = NULL;
	globali.connectivity_modular.Sc_v = NULL;
	globali.connectivity_modular.Cg = NULL;
	globali.connectivity_modular.Sin = NULL;
	globali.connectivity_modular.Sout = NULL;
	globali.connectivity_modular.a1 = NULL;
	globali.connectivity_modular.a2 = NULL;
	globali.connectivity_modular.ind = NULL;
	globali.connectivity_modular.ind_Sc = NULL;
	globali.connectivity_modular.h_new = NULL;
	globali.connectivity_modular.mod_type = NULL;
	globali.connectivity_modular.ind_s = NULL;
	globali.connectivity_modular.aus0 = NULL;

	globali.simulateprofiles.ind = NULL;
	globali.simulateprofiles.tmpm1_i = NULL;
	globali.simulateprofiles.tmpm1_d = NULL;
	globali.simulateprofiles.reg = NULL;
	globali.simulateprofiles.M = NULL;
	globali.simulateprofiles.Mdiscr = NULL;
	globali.simulateprofiles.Mneg = NULL;
	globali.simulateprofiles.genenet = NULL;

	globali.simulatenet.M = NULL;
	globali.simulatenet.D = NULL;
	globali.simulatenet.Mneg = NULL;
	globali.simulatenet.Mdiscr = NULL;
	globali.simulatenet.reg = NULL;
	globali.simulatenet.aus = NULL;
	globali.simulatenet.genenet = NULL;
	globali.simulatenet.nulla = NULL;
	globali.simulatenet.R = NULL;
	globali.simulatenet.ris = NULL;
	globali.simulatenet.tmpm1_d = NULL;

}

void CancGlobali()
{
	libera(globali.tipi);

	CANCELLAv_i(globali.sample.HL);
	CANCELLAv_d(globali.sample.q);
	CANCELLAv_i(globali.sample.x);
	CANCELLAv_i(globali.sample.y);

	CANCELLAv_d(globali.ret_vettore.v2);
	CANCELLAv_d(globali.ret_vettore.v3);

	CANCELLAv_i(globali.ret_lista.v2);
	// questa e` stata cancellata da "daLISTA"
	//~ CANCELLAm_i(globali.ret_lista.m1);

	CANCELLAm_i(globali.check_conn.Mt);
	CANCELLAm_i(globali.check_conn.Maus);
	CANCELLAv_i(globali.check_conn.tmp_i);
	CANCELLAv_i(globali.check_conn.colore);
	CANCELLAv_i(globali.check_conn.grigi);
	CANCELLAv_i(globali.check_conn.ind);
	CANCELLAv_i(globali.check_conn.adj);
	CANCELLAv_d(globali.check_conn.tmp1_d);
	CANCELLAv_d(globali.check_conn.scalare_d);

	CANCELLAv_d(globali.boole_result.valr);
	CANCELLAv_i(globali.boole_result.signr);
	CANCELLAv_i(globali.boole_result.tmp_r);

	CANCELLAv_i(globali.triangola.tmp1_i);
	CANCELLAv_i(globali.triangola.tmp2_i);
	CANCELLAv_d(globali.triangola.tmp3_d);
	CANCELLAv_i(globali.triangola.M_in);
	CANCELLAv_i(globali.triangola.sk);
	CANCELLAv_i(globali.triangola.coord_ind);
	CANCELLAv_i(globali.triangola.ind_aus);
	CANCELLAv_i(globali.triangola.ind);
	CANCELLAm_i(globali.triangola.coord);
	CANCELLAm_i(globali.triangola.tmp_coord);
	CANCELLAm_i(globali.triangola.tmp_coord1);
	CANCELLAv_i(globali.triangola.Dmem);

	CANCELLAm_i(globali.cluster_coeff.tmpm_i);
	CANCELLAv_i(globali.cluster_coeff.ind);

	CANCELLAv_i(globali.assign_nodes.Sin_h);
	CANCELLAv_i(globali.assign_nodes.tmp_i);
	CANCELLAv_d(globali.assign_nodes.tmp_d);
	CANCELLAv_i(globali.assign_nodes.or_h);
	CANCELLAv_i(globali.assign_nodes.aus_h);
	CANCELLAv_i(globali.assign_nodes.M_in);
	CANCELLAv_i(globali.assign_nodes.Ord);
	CANCELLAv_d(globali.assign_nodes.p);
	CANCELLAv_i(globali.assign_nodes.ind_h);
	CANCELLAv_i(globali.assign_nodes.ri);
	CANCELLAv_i(globali.assign_nodes.co);
	CANCELLAv_i(globali.assign_nodes.ind);
	CANCELLAv_i(globali.assign_nodes.index);

	CANCELLAv_d(globali.target.aus);
	CANCELLAm_d(globali.target.tmp_d);
	CANCELLAv_d(globali.target.tmp1_d);
	CANCELLAv_i(globali.target.ind);
	CANCELLAv_i(globali.target.ind_reg);

	CANCELLAm_d(globali.scoremodular.tmpm_d);
	CANCELLAv_i(globali.scoremodular.tmp_i1);
	CANCELLAv_d(globali.scoremodular.tmp_d1);
	CANCELLAv_d(globali.scoremodular.T1);
	CANCELLAv_d(globali.scoremodular.T2);
	CANCELLAv_d(globali.scoremodular.old1);
	CANCELLAv_d(globali.scoremodular.old2);
	CANCELLAv_d(globali.scoremodular.new1);
	CANCELLAv_d(globali.scoremodular.new2);
	CANCELLAv_d(globali.scoremodular.toll1);
	CANCELLAv_d(globali.scoremodular.toll2);
	CANCELLAv_d(globali.scoremodular.a);
	CANCELLAv_d(globali.scoremodular.b);
	CANCELLAv_i(globali.scoremodular.ind);
	CANCELLAv_d(globali.scoremodular.m);
	CANCELLAv_d(globali.scoremodular.S1);
	CANCELLAv_d(globali.scoremodular.S2);
	CANCELLAv_i(globali.scoremodular.ind1);
	CANCELLAv_i(globali.scoremodular.ind2);

	CANCELLAm_d(globali.score_sf.tmpm_d);
	CANCELLAv_d(globali.score_sf.scalare_d);
	CANCELLAv_i(globali.score_sf.tmp_i1);
	CANCELLAv_d(globali.score_sf.tmp_d1);
	CANCELLAv_d(globali.score_sf.T1);
	CANCELLAv_d(globali.score_sf.T2);
	CANCELLAv_d(globali.score_sf.old1);
	CANCELLAv_d(globali.score_sf.old2);
	CANCELLAv_d(globali.score_sf.new1);
	CANCELLAv_d(globali.score_sf.new2);
	CANCELLAv_d(globali.score_sf.toll1);
	CANCELLAv_d(globali.score_sf.toll2);
	CANCELLAv_d(globali.score_sf.a);
	CANCELLAv_d(globali.score_sf.b);
	CANCELLAv_i(globali.score_sf.ind);
	CANCELLAv_d(globali.score_sf.m);
	CANCELLAv_d(globali.score_sf.S1);
	CANCELLAv_d(globali.score_sf.S2);
	CANCELLAv_i(globali.score_sf.ind1);
	CANCELLAv_i(globali.score_sf.ind2);
	CANCELLAv_i(globali.score_sf.indinf);
	CANCELLAv_i(globali.score_sf.tmp_i2);
	CANCELLAv_i(globali.score_sf.indbad);
	CANCELLAv_i(globali.score_sf.ind0);
	CANCELLAv_i(globali.score_sf.ind3);

	CANCELLAv_d(globali.probmod.S_out);
	CANCELLAv_d(globali.probmod.S_in);
	CANCELLAm_d(globali.probmod.score_matr1);
	CANCELLAm_d(globali.probmod.score_matr2);
	CANCELLAm_d(globali.probmod.score_matr3);
	CANCELLAm_i(globali.probmod.checkIN);
	CANCELLAm_i(globali.probmod.checkOUT);
	CANCELLAm_i(globali.probmod.memory);
	CANCELLAv_i(globali.probmod.M_in);
	CANCELLAv_i(globali.probmod.M_out);
	CANCELLAv_i(globali.probmod.tmp_i1);
	CANCELLAv_i(globali.probmod.tmp_i2);
	CANCELLAv_i(globali.probmod.tmp_i3);
	CANCELLAv_i(globali.probmod.indInf);
	CANCELLAv_i(globali.probmod.I);
	CANCELLAv_i(globali.probmod.ord_ind);
	CANCELLAv_i(globali.probmod.rs);
	CANCELLAv_i(globali.probmod.ind1);
	CANCELLAv_i(globali.probmod.scalare_i);
	CANCELLAv_i(globali.probmod.I_add);
	CANCELLAv_d(globali.probmod.tmp_d1);

	CANCELLAm_d(globali.dinamica.D);
	CANCELLAv_i(globali.dinamica.ind);
	CANCELLAv_d(globali.dinamica.y_prec);
	CANCELLAv_d(globali.dinamica.targetT);
	CANCELLAv_d(globali.dinamica.targ);
	CANCELLAv_d(globali.dinamica.incr);
	CANCELLAv_d(globali.dinamica.n);
	CANCELLAv_d(globali.dinamica.aus);
	CANCELLAm_d(globali.dinamica.tmp_ris);

	CANCELLAv_i(globali.createNEG.tmp_i);
	CANCELLAv_i(globali.createNEG.segno);
	CANCELLAv_i(globali.createNEG.ind);

	CANCELLAv_d(globali.lsoda.y);
	CANCELLAv_d(globali.lsoda.y1);
	CANCELLAv_d(globali.lsoda.y_prec);
	CANCELLAv_d(globali.lsoda.targetT);
	CANCELLAv_d(globali.lsoda.targ);
	CANCELLAm_d(globali.lsoda.tmp_ris);

	CANCELLAv_i(globali.createRules.op);
	CANCELLAv_i(globali.createRules.tmp);

	CANCELLAv_d(globali.create_logicRule.x);
	CANCELLAv_i(globali.create_logicRule.s);
	CANCELLAv_d(globali.create_logicRule.scalare_d);
	CANCELLAv_i(globali.create_logicRule.nvect);
	CANCELLAv_i(globali.create_logicRule.e);
	CANCELLAv_d(globali.create_logicRule.prob);
	CANCELLAv_i(globali.create_logicRule.tmp_i);
	CANCELLAv_i(globali.create_logicRule.tmp1_i);
	CANCELLAv_d(globali.create_logicRule.pr_and);
	CANCELLAv_d(globali.create_logicRule.pr_or);
	CANCELLAv_i(globali.create_logicRule.scalare_i);
	CANCELLAv_i(globali.create_logicRule.blacklist);
	CANCELLAv_i(globali.create_logicRule.black_p);
	CANCELLAv_i(globali.create_logicRule.o);

	CANCELLAm_i(globali.mod1.ind);
	CANCELLAv_i(globali.mod1.tmp1_i);
	CANCELLAv_i(globali.mod1.x);
	CANCELLAv_i(globali.mod1.s);

	CANCELLAv_i(globali.mod3.M_out);
	CANCELLAv_i(globali.mod3.indok);
	CANCELLAv_i(globali.mod3.tmp3_i);
	CANCELLAv_d(globali.mod3.tmp1_d);
	CANCELLAv_d(globali.mod3.tmpSTin);
	CANCELLAv_d(globali.mod3.tmpSTout);
	CANCELLAv_d(globali.mod3.Freq_in);
	CANCELLAv_d(globali.mod3.Freq_out);
	CANCELLAv_d(globali.mod3.tmp3_d);
	CANCELLAv_d(globali.mod3.p);
	CANCELLAv_d(globali.mod3.Sc);
	CANCELLAv_i(globali.mod3.tmp2_i);
	CANCELLAv_i(globali.mod3.ind_M);
	CANCELLAv_i(globali.mod3.ind1);
	CANCELLAv_i(globali.mod3.indS);
	CANCELLAv_i(globali.mod3.indBS);
	CANCELLAv_i(globali.mod3.Sin);
	CANCELLAv_i(globali.mod3.ind);
	CANCELLAv_i(globali.mod3.indInf);
	CANCELLAv_i(globali.mod3.tmp1_i);
	CANCELLAv_i(globali.mod3.scalare_i);
	CANCELLAv_d(globali.mod3.scalare_d);
	CANCELLAv_d(globali.mod3.p_sc);

	CANCELLAv_i(globali.module1.scalare_i);
	CANCELLAv_i(globali.module1.Ng);
	CANCELLAm_i(globali.module1.conn_matr);
	CANCELLAv_i(globali.module1.indices);
	CANCELLAv_i(globali.module1.tmp1_i);

	CANCELLAm_i(globali.module2.conn_matr);
	CANCELLAv_i(globali.module2.indices);
	CANCELLAv_i(globali.module2.tmp1_i);
	CANCELLAv_i(globali.module2.Ng);

	CANCELLAm_i(globali.module3.conn_matr);
	CANCELLAv_i(globali.module3.indices);
	CANCELLAv_i(globali.module3.tmp1_i);
	CANCELLAv_i(globali.module3.Ng);
	CANCELLAv_i(globali.module3.Ng_UP);

	CANCELLAm_d(globali.connectivity_geometric.M);
	CANCELLAm_i(globali.connectivity_geometric.Mdiscr);
	CANCELLAv_d(globali.connectivity_geometric.x);
	CANCELLAv_d(globali.connectivity_geometric.y);
	CANCELLAv_d(globali.connectivity_geometric.xy);
	CANCELLAv_d(globali.connectivity_geometric.d);
	CANCELLAv_i(globali.connectivity_geometric.s);
	CANCELLAv_i(globali.connectivity_geometric.regulatedind);
	CANCELLAv_i(globali.connectivity_geometric.indL);
	CANCELLAv_i(globali.connectivity_geometric.Sr);
	CANCELLAv_d(globali.connectivity_geometric.aus);
	CANCELLAv_i(globali.connectivity_geometric.ind);
	CANCELLAv_i(globali.connectivity_geometric.ind1);
	CANCELLAv_i(globali.connectivity_geometric.ind0);
	CANCELLAv_i(globali.connectivity_geometric.tmp1_i);
	CANCELLAv_i(globali.connectivity_geometric.tmp2_i);
	CANCELLAv_d(globali.connectivity_geometric.tmp1_d);

	CANCELLAm_d(globali.connectivity_scalefree.M);
	CANCELLAm_i(globali.connectivity_scalefree.Mdiscr);
	CANCELLAv_d(globali.connectivity_scalefree.x);
	CANCELLAv_d(globali.connectivity_scalefree.y);
	CANCELLAv_d(globali.connectivity_scalefree.d);
	CANCELLAv_i(globali.connectivity_scalefree.s);
	CANCELLAv_i(globali.connectivity_scalefree.o);
	CANCELLAv_i(globali.connectivity_scalefree.regulatedind);
	CANCELLAv_i(globali.connectivity_scalefree.indL);
	CANCELLAv_i(globali.connectivity_scalefree.Sr);
	CANCELLAv_i(globali.connectivity_scalefree.inthenet);
	CANCELLAv_i(globali.connectivity_scalefree.not_inthenet);
	CANCELLAv_i(globali.connectivity_scalefree.not_regulated);
	CANCELLAv_i(globali.connectivity_scalefree.numposs);
	CANCELLAv_i(globali.connectivity_scalefree.give_outlink);
	CANCELLAv_i(globali.connectivity_scalefree.aus_give_outlink);
	CANCELLAv_i(globali.connectivity_scalefree.indInf);
	CANCELLAv_i(globali.connectivity_scalefree.ind_Sc);
	CANCELLAv_i(globali.connectivity_scalefree.a1);
	CANCELLAv_i(globali.connectivity_scalefree.a2);
	CANCELLAv_i(globali.connectivity_scalefree.primi);
	CANCELLAv_i(globali.connectivity_scalefree.indici);
	CANCELLAv_i(globali.connectivity_scalefree.num_v);
	CANCELLAv_i(globali.connectivity_scalefree.mem_o);
	CANCELLAv_i(globali.connectivity_scalefree.available);
	CANCELLAv_i(globali.connectivity_scalefree.campione);
	CANCELLAv_i(globali.connectivity_scalefree.linked);
	CANCELLAv_i(globali.connectivity_scalefree.ind_s);
	CANCELLAv_i(globali.connectivity_scalefree.Sout);
	CANCELLAv_i(globali.connectivity_scalefree.Sin);
	CANCELLAv_i(globali.connectivity_scalefree.ind);
	CANCELLAv_i(globali.connectivity_scalefree.indok);
	CANCELLAv_i(globali.connectivity_scalefree.ind1);
	CANCELLAv_i(globali.connectivity_scalefree.ind0);
	CANCELLAv_i(globali.connectivity_scalefree.tmp1_i);
	CANCELLAv_i(globali.connectivity_scalefree.tmp2_i);
	CANCELLAv_d(globali.connectivity_scalefree.tmp1_d);
	CANCELLAv_d(globali.connectivity_scalefree.tmp2_d);
	CANCELLAv_i(globali.connectivity_scalefree.scalare_i);
	CANCELLAv_d(globali.connectivity_scalefree.scalare_d);
	CANCELLAv_d(globali.connectivity_scalefree.Prob);
	CANCELLAv_d(globali.connectivity_scalefree.Freq_in);
	CANCELLAv_d(globali.connectivity_scalefree.Freq_out);
	CANCELLAv_d(globali.connectivity_scalefree.STin);
	CANCELLAv_d(globali.connectivity_scalefree.STout);
	CANCELLAv_d(globali.connectivity_scalefree.p);
	CANCELLAv_i(globali.connectivity_scalefree.toll1);
	CANCELLAv_d(globali.connectivity_scalefree.toll_in);
	CANCELLAv_d(globali.connectivity_scalefree.toll_out);
	CANCELLAv_d(globali.connectivity_scalefree.Sc);
	CANCELLAv_d(globali.connectivity_scalefree.p_ind);
	CANCELLAv_d(globali.connectivity_scalefree.p_out);
	CANCELLAm_d(globali.connectivity_scalefree.aus);

	CANCELLAm_d(globali.connectivity_random.M);
	CANCELLAm_i(globali.connectivity_random.Mdiscr);
	CANCELLAv_d(globali.connectivity_random.aus);
	CANCELLAv_d(globali.connectivity_random.tmp_d);
	CANCELLAv_i(globali.connectivity_random.scalare_i);
	CANCELLAv_d(globali.connectivity_random.Pnum);
	CANCELLAv_i(globali.connectivity_random.ind);
	CANCELLAv_i(globali.connectivity_random.num);
	CANCELLAv_i(globali.connectivity_random.tmp_i);

	CANCELLAm_d(globali.connectivity_modular.M);
	CANCELLAm_i(globali.connectivity_modular.Mdiscr);
	CANCELLAv_i(globali.connectivity_modular.scalare_i);
	CANCELLAv_i(globali.connectivity_modular.tmp1_i);
	CANCELLAv_i(globali.connectivity_modular.tmp2_i);
	CANCELLAv_d(globali.connectivity_modular.tmp1_d);
	CANCELLAv_d(globali.connectivity_modular.tmp2_d);
	CANCELLAv_d(globali.connectivity_modular.Prob);
	CANCELLAv_d(globali.connectivity_modular.Freq_out);
	CANCELLAv_d(globali.connectivity_modular.Freq_in);
	CANCELLAv_d(globali.connectivity_modular.STout);
	CANCELLAv_d(globali.connectivity_modular.STin);
	CANCELLAv_d(globali.connectivity_modular.toll);
	CANCELLAv_d(globali.connectivity_modular.p_out);
	CANCELLAv_d(globali.connectivity_modular.scalare_d);
	CANCELLAv_i(globali.connectivity_modular.Sr);
	CANCELLAv_d(globali.connectivity_modular.p);
	CANCELLAv_i(globali.connectivity_modular.h);
	CANCELLAv_d(globali.connectivity_modular.prob_mod);
	CANCELLAm_d(globali.connectivity_modular.aus);
	CANCELLAv_d(globali.connectivity_modular.Sc_v);
	CANCELLAv_d(globali.connectivity_modular.Cg);
	CANCELLAv_i(globali.connectivity_modular.Sin);
	CANCELLAv_i(globali.connectivity_modular.Sout);
	CANCELLAv_i(globali.connectivity_modular.a1);
	CANCELLAv_i(globali.connectivity_modular.a2);
	CANCELLAv_i(globali.connectivity_modular.ind);
	CANCELLAv_i(globali.connectivity_modular.ind_Sc);
	CANCELLAv_i(globali.connectivity_modular.h_new);
	CANCELLAv_i(globali.connectivity_modular.mod_type);
	CANCELLAv_i(globali.connectivity_modular.ind_s);
	CANCELLAv_i(globali.connectivity_modular.aus0);

	CANCELLAv_i(globali.simulateprofiles.ind);
	CANCELLAm_i(globali.simulateprofiles.tmpm1_i);
	CANCELLAm_d(globali.simulateprofiles.tmpm1_d);
	CANCELLAm_i(globali.simulateprofiles.reg);
	// in questo caso M != conn*.M
	CANCELLAm_d(globali.simulateprofiles.M);
	// in questo caso Mdiscr != conn*.Mdiscr
	CANCELLAm_i(globali.simulateprofiles.Mdiscr);
	CANCELLAm_i(globali.simulateprofiles.Mneg);

	// M == conn*.M
	//~ CANCELLAm_d(globali.simulatenet.M);
	// questa e` stata cancellata da "daLISTA">
	//~ CANCELLAm_d(globali.simulatenet.tmpm1_d);
	// anche questa e` stata cancellata da "daLISTA"
	//~ CANCELLAm_d(globali.simulatenet.D);
	CANCELLAm_i(globali.simulatenet.Mneg);
	// Mdiscr == conn*.Mdiscr
	//~ CANCELLAm_i(globali.simulatenet.Mdiscr);
	CANCELLAm_i(globali.simulatenet.reg);

}
