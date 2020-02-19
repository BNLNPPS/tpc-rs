/* g2t_track.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
 Id: g2t_track.idl,v 1.22 20190930 14:13:56 jwebb Exp
     Log: g2t_track.idl,v
     Revision 1.22  20190930 14:13:56  jwebb
     Integrate HITS for forward tracking and forward calorimeter.

     n.b. deprecates the legacy HcalGeo RnD detector.

     Revision 1.21  20171002 15:29:39  jwebb
     Integration of ETOF into simulation

     Revision 1.20  20160918 22:37:52  fisyak
     pack no. of real TPC hits

     Revision 1.19  20160918 22:35:56  fisyak

     Revision 1.18  20151012 20:46:58  jwebb
     Hit definition and starsim to root interface for FTS.

     Revision 1.17  20140514 20:01:23  jwebb
     More support for HCAL.  Also note... last checkin of g2t_volume_id was to support FMS preshower.

     Revision 1.16  20120124 03:36:25  perev
     Add Etr

     Revision 1.15  20110803 20:11:57  jwebb
     Add MTD to the g2t hit tables.

     Revision 1.14  20110720 20:44:44  perev
     Fsc added

     Revision 1.13  20110126 19:21:17  perev
     FPD > STAR Soft

     Revision 1.12  20060922 20:16:42  potekhin
     Added HPD (RD)

     Revision 1.11  20060626 19:58:35  potekhin
     Include the gem barrel detector

     Revision 1.10  20050630 16:52:16
COMMENTS TRUNCATED */
#ifndef G2T_TRACK_H
#define G2T_TRACK_H
struct g2t_track_st {
  int eg_label; /* generator track label (0 if GEANT track) */
  int eg_pid; /* event generator particle id */
  int ge_pid; /* GEANT particle id */
  int hit_ctb_p; /* Id of first ctb hit on track linked list */
  int hit_eem_p; /* Id of first eem hit on track linked list */
  int hit_emc_p; /* Id of first emc hit on track linked list */
  int hit_esm_p; /* Id of first esm hit on track linked list */
  int hit_ftp_p; /* Id of first ftp hit on track linked list */
  int hit_gem_p; /* Id of first gem hit on track linked list */
  int hit_hpd_p; /* Id of first hpd hit on track linked list */
  int hit_ist_p; /* Id of first ist hit on track linked list */
  int hit_igt_p; /* Id of first igt hit on track linked list */
  int hit_fst_p; /* Id of first fst hit on track linked list */
  int hit_fgt_p; /* Id of first fgt hit on track linked list */
  int hit_fpd_p; /* Id of first fpd hit on track linked list */
  int hit_fsc_p; /* Id of first fsc hit on track linked list */
  int hit_mtd_p; /* Id of first mtd hit on track linked list */
  int hit_mwc_p; /* Id of first mwc hit on track linked list */
  int hit_pgc_p; /* Id of first pgc hit on track linked list */
  int hit_pmd_p; /* Id of first psc hit on track linked list */
  int hit_smd_p; /* pointer to first SHM linked list hit */
  int hit_ssd_p; /* Id of first ssd hit on track linked list */
  int hit_svt_p; /* Id of first svt hit on track linked list */
  int hit_pix_p; /* Id of first pix hit on track linked list */
  int hit_tof_p; /* Id of first tof hit on track linked list */
  int hit_tpc_p; /* Id of first tpc hit on track linked list */
  int hit_vpd_p; /* Id of first vpd hit on track linked list */
  int hit_etr_p; /* Id of first etr hit on track linked list */
  int hit_hca_p; /* Id of first hca hit on track linked list */
  int hit_fts_p; /* Id of first fts hit on track linked list */
  int hit_eto_p; /* Id of first etof hit on track linked list */
  int id; /* primary key */
  int is_shower; /* 1 if shower track, 0 if not */
  int itrmd_vertex_p; /* First intermediate vertex */
  int n_ctb_hit; /* Nhits in ctb */
  int n_eem_hit; /* Nhits in eem (endcap em cal) */
  int n_emc_hit; /* Nhits in emc */
  int n_esm_hit; /* Nhits in esm (endcap shower max) */
  int n_ftp_hit; /* Nhits in forward tpc */
  int n_gem_hit; /* Nhits in gem barrel */
  int n_hpd_hit; /* Nhits in hpd */
  int n_ist_hit; /* Nhits in ist */
  int n_igt_hit; /* Nhits in igt */
  int n_fst_hit; /* Nhits in fst [new Forward Silicon Tracker] */
  int n_fgt_hit; /* Nhits in fgt */
  int n_fpd_hit; /* Nhits in fpd */
  int n_fsc_hit; /* Nhits in fsc */
  int n_mtd_hit; /* Nhits in mtd */
  int n_mwc_hit; /* Nhits in mwc */
  int n_pgc_hit; /* Nhits in pgc  ???  */
  int n_pmd_hit; /* Nhits in pmd (PMD) */
  int n_smd_hit; /* number of hits in shower max */
  int n_ssd_hit; /* Nhits in ssd */
  int n_svt_hit; /* Nhits in svt */
  int n_pix_hit; /* Nhits in pix */
  int n_tof_hit; /* Nhits in tof */
  int n_tpc_hit; /* Nhits in tpc + (no. of real hits after TpcRS) << 8 */
  int n_vpd_hit; /* Nhits in vpd */
  int n_etr_hit; /* Nhits on etr */
  int n_hca_hit; /* Nhits on hca [hadronic calorimeter]*/
  int n_fts_hit; /* Nhits on fts */
  int n_eto_hit; /* Nhits on etof */
  int n_stg_hit; /* Nhits on sTGC [new Forward small thin gap chambers]*/
  int n_wca_hit; /* Nhits on wca [West EM calorimeter]*/
  int next_parent_p; /* Id of next parent track */
  int next_vtx_trk_p; /* Next daughter track of start vertex */
  int start_vertex_p; /* Id of start vertex of track */
  int stop_vertex_p; /* Id of stop vertex of this track */
  float charge; /* Charge */
  float e; /* Energy */
  float eta; /* Pseudorapidity */
  float p[3]; /* Momentum */
  float pt; /* Transverse momentum */
  float ptot; /* Total momentum */
  float rapidity; /* Rapidity */
};
#endif /* G2T_TRACK_H */
