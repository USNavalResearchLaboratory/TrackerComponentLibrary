/*

-Header_File npf_cspice.h ( Mice routine name pointer mapping structure )

-Abstract

   List the Mice call name to CSPICE call mappings.

-Disclaimer

   THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE
   CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.
   GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE
   ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE
   PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS"
   TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY
   WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A
   PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC
   SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE
   SOFTWARE AND RELATED MATERIALS, HOWEVER USED.

   IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA
   BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT
   LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,
   INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS,
   REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE
   REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.

   RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF
   THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY
   CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE
   ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.

-Required_Reading

   MICE.REQ

-Keywords

-Brief_I/O

   None.

-Detailed_Input

   None.

-Detailed_Output

   None.

-Parameters
*/

#ifndef npf_cspice_h
#define npf_cspice_h

/*
-Exceptions

   None.

-Files

   None.

-Particulars

   Below, the "npf" structure and the assignments within. This
   structure defines mapping between the routine name field "name"
   to the pointer corresponding to that routine's interface.

   A "name" with format name_c maps to a cspice_name call, a "name"
   name_s maps to a mice_name and cspice_name call from MATLAB where
   "name" represents the call base name.

   Update this list to add a new interface call to Mice.

-Examples

   None.

-Restrictions

   None.

-Literature_References

   None.

-Author_and_Institution

   E.D. Wright        (JPL)
   G. Chinn           (JPL)

-Version

   -Mice Mice 1.3.0 31-OCT-2012 (EDW)(SCK)

      Added interfaces:

         cspice_cidfrm (mice_cidfrm)
         cspice_cnmfrm (mice_cnmfrm)
         cspice_dafac
         cspice_dafbbs
         cspice_dafbfs
         cspice_dafcls
         cspice_dafcs
         cspice_dafdc
         cspice_dafec
         cspice_daffna
         cspice_daffpa
         cspice_dafgda
         cspice_dafgn
         cspice_dafgs
         cspice_dafopr
         cspice_dafopw
         cspice_dafus
         cspice_dcyldr
         cspice_dgeodr
         cspice_dlatdr
         cspice_dpgrdr
         cspice_drdcyl
         cspice_drdgeo
         cspice_drdlat
         cspice_drdpgr
         cspice_drdsph
         cspice_dsphdr
         cspice_dvpool
         cspice_edlimb (mice_edlimb)
         cspice_edterm
         cspice_fovray
         cspice_fovtrg
         cspice_frame
         cspice_frinfo (mice_frinfo)
         cspice_frmnam
         cspice_gfilum
         cspice_gfpa
         cspice_gfstol
         cspice_inedpl
         cspice_inelpl
         cspice_inrypl
         cspice_invort
         cspice_namfrm
         cspice_npedln (mice_npedln)
         cspice_npelpt (mice_npelpt)
         cspice_nplnpt (mice_nplnpt)
         cspice_occult
         cspice_phaseq
         cspice_pjelpl (mice_pjelpl)
         cspice_pl2nvc
         cspice_pl2nvp
         cspice_pl2psv
         cspice_psv2pl
         cspice_pxfrm2
         cspice_spkcls
         cspice_spkcpo (mice_spkcpo)
         cspice_spkcpt (mice_spkcpt)
         cspice_spkcvo (mice_spkcvo)
         cspice_spkcvt (mice_spkcvt)
         cspice_spkopn
         cspice_spkpvn
         cspice_spksfs
         cspice_spkw08
         cspice_surfpt (mice_surfpt)
         cspice_timdef_get
         cspice_timdef_set
         cspice_vprjp
         cspice_vprjpi
         cspice_vproj
         cspice_xfmsta

   -Mice 1.2.0 28-APR-2010 (EDW)

      Added interfaces:

         cspice_bodc2s (mice_bodc2s)
         cspice_ducrss
         cspice_dvcrss
         cspice_dvdot
         cspice_dvhat
         cspice_dvnorm
         cspice_dvsep
         cspice_ekfind
         cspice_ekgc
         cspice_ekgd
         cspice_ekgi
         cspice_eknelt
         cspice_gfrr
         cspice_vperp
         cspice_vnorm
         cspice_unitim

   -Mice 1.1.0 23-FEB-2009 (EDW)

      Added interfaces:

         cspice_el2cgv
         cspice_gfdist
         cspice_gfposc
         cspice_gfrfov
         cspice_gfsep
         cspice_gfsntc
         cspice_gfsubc
         cspice_gftfov
         cspice_lmpool
         cspice_nvc2pl
         cspice_nvp2pl
         cspice_saelgv
         cspice_vrotv
         cspice_wncard
         cspice_wnsumd (mice_wnsumd)
         cspice_wnvald

   -Mice 1.0.0 12-FEB-2008 (EDW)

-Index_Entries

*/

static struct npf
   {
   char  * name;
   void    (*pfunc)(int, mxArray*[], int, const mxArray*[]);
   }
   NPF[] =
     {
     { "axisar_c", &cspice_axisar },
     { "b1900_c",  &cspice_b1900  },
     { "b1950_c",  &cspice_b1950  },
     { "bodc2n_s", &mice_bodc2n   },
     { "bodc2s_s", &mice_bodc2s   },
     { "boddef_c", &cspice_boddef },
     { "bodn2c_s", &mice_bodn2c   },
     { "bods2c_s", &mice_bods2c   },
     { "bodvcd_c", &cspice_bodvcd },
     { "bodvrd_c", &cspice_bodvrd },
     { "cgv2el_c", &cspice_cgv2el },
     { "cidfrm_s", &mice_cidfrm   },
     { "ckcls_c",  &cspice_ckcls  },
     { "ckcov_c",  &cspice_ckcov  },
     { "ckgp_c",   &cspice_ckgp   },
     { "ckgpav_c", &cspice_ckgpav },
     { "ckobj_c",  &cspice_ckobj  },
     { "ckopn_c",  &cspice_ckopn  },
     { "ckw01_c",  &cspice_ckw01  },
     { "ckw02_c",  &cspice_ckw02  },
     { "ckw03_c",  &cspice_ckw03  },
     { "clight_c", &cspice_clight },
     { "clpool_c", &cspice_clpool },
     { "cnmfrm_s", &mice_cnmfrm   },
     { "conics_c", &cspice_conics },
     { "convrt_c", &cspice_convrt },
     { "cspice_mice", &cspice_mice},
     { "cyllat_c", &cspice_cyllat },
     { "cylrec_c", &cspice_cylrec },
     { "cylsph_c", &cspice_cylsph },
     { "dafac_c",  &cspice_dafac  },
     { "dafbbs_c", &cspice_dafbbs },
     { "dafbfs_c", &cspice_dafbfs },
     { "dafcls_c", &cspice_dafcls },
     { "dafcs_c",  &cspice_dafcs  },
     { "dafdc_c",  &cspice_dafdc  },
     { "dafec_c",  &cspice_dafec  },
     { "daffna_c", &cspice_daffna },
     { "daffpa_c", &cspice_daffpa },
     { "dafgda_c", &cspice_dafgda },
     { "dafgn_c",  &cspice_dafgn  },
     { "dafgs_c",  &cspice_dafgs  },
     { "dafopr_c", &cspice_dafopr },
     { "dafopw_c", &cspice_dafopw },
     { "dafus_c",  &cspice_dafus  },
     { "dcyldr_c", &cspice_dcyldr },
     { "deltet_c", &cspice_deltet },
     { "dgeodr_c", &cspice_dgeodr },
     { "dlatdr_c", &cspice_dlatdr },
     { "dpgrdr_c", &cspice_dpgrdr },
     { "dpr_c",    &cspice_dpr    },
     { "drdcyl_c", &cspice_drdcyl },
     { "drdgeo_c", &cspice_drdgeo },
     { "drdlat_c", &cspice_drdlat },
     { "drdpgr_c", &cspice_drdpgr },
     { "drdsph_c", &cspice_drdsph },
     { "dsphdr_c", &cspice_dsphdr },
     { "dtpool_s", &mice_dtpool   },
     { "ducrss_c", &cspice_ducrss },
     { "dvcrss_c", &cspice_dvcrss },
     { "dvdot_c",  &cspice_dvdot  },
     { "dvhat_c",  &cspice_dvhat  },
     { "dvnorm_c", &cspice_dvnorm },
     { "dvpool_c", &cspice_dvpool },
     { "dvsep_c",  &cspice_dvsep  },
     { "edlimb_s", &mice_edlimb   },
     { "edterm_c", &cspice_edterm },
     { "ekfind_c", &cspice_ekfind },
     { "ekgc_c",   &cspice_ekgc   },
     { "ekgd_c",   &cspice_ekgd   },
     { "ekgi_c",   &cspice_ekgi   },
     { "eknelt_c", &cspice_eknelt },
     { "el2cgv_c", &cspice_el2cgv },
     { "et2lst_c", &cspice_et2lst },
     { "et2utc_c", &cspice_et2utc },
     { "etcal_c",  &cspice_etcal  },
     { "eul2m_c",  &cspice_eul2m  },
     { "eul2xf_c", &cspice_eul2xf },
     { "fovray_c", &cspice_fovray },
     { "fovtrg_c", &cspice_fovtrg },
     { "frame_c",  &cspice_frame  },
     { "frinfo_s", &mice_frinfo   },
     { "frmnam_c", &cspice_frmnam },
     { "furnsh_c", &cspice_furnsh },
     { "gcpool_c", &cspice_gcpool },
     { "gdpool_c", &cspice_gdpool },
     { "georec_c", &cspice_georec },
     { "getfov_c", &cspice_getfov },
     { "gfdist_c", &cspice_gfdist },
     { "gfilum_c", &cspice_gfilum },
     { "gfoclt_c", &cspice_gfoclt },
     { "gfpa_c",   &cspice_gfpa   },
     { "gfposc_c", &cspice_gfposc },
     { "gfrfov_c", &cspice_gfrfov },
     { "gfrr_c",   &cspice_gfrr   },
     { "gfsep_c",  &cspice_gfsep  },
     { "gfsntc_c", &cspice_gfsntc },
     { "gfstol_c", &cspice_gfstol },
     { "gfsubc_c", &cspice_gfsubc },
     { "gftfov_c", &cspice_gftfov },
     { "gipool_c", &cspice_gipool },
     { "gnpool_c", &cspice_gnpool },
     { "halfpi_c", &cspice_halfpi },
     { "illum_c",  &cspice_illum  },
     { "ilumin_s", &mice_ilumin   },
     { "inedpl_c", &cspice_inedpl },
     { "inelpl_c", &cspice_inelpl },
     { "inrypl_c", &cspice_inrypl },
     { "invort_c", &cspice_invort },
     { "j1900_c",  &cspice_j1900  },
     { "j1950_c",  &cspice_j1950  },
     { "j2000_c",  &cspice_j2000  },
     { "j2100_c",  &cspice_j2100  },
     { "jyear_c",  &cspice_jyear  },
     { "kclear_c", &cspice_kclear },
     { "kdata_c",  &cspice_kdata  },
     { "kinfo_c",  &cspice_kinfo  },
     { "ktotal_c", &cspice_ktotal },
     { "latcyl_c", &cspice_latcyl },
     { "latrec_c", &cspice_latrec },
     { "latsph_c", &cspice_latsph },
     { "lmpool_c", &cspice_lmpool },
     { "ltime_c",  &cspice_ltime  },
     { "m2eul_c",  &cspice_m2eul  },
     { "m2q_c",    &cspice_m2q    },
     { "namfrm_c", &cspice_namfrm },
     { "nearpt_s", &mice_nearpt   },
     { "npedln_s", &mice_npedln   },
     { "npelpt_s", &mice_npelpt   },
     { "nplnpt_s", &mice_nplnpt   },
     { "nvc2pl_c", &cspice_nvc2pl },
     { "nvp2pl_c", &cspice_nvp2pl },
     { "occult_c", &cspice_occult },
     { "oscelt_c", &cspice_oscelt },
     { "pcpool_c", &cspice_pcpool },
     { "pdpool_c", &cspice_pdpool },
     { "pgrrec_c", &cspice_pgrrec },
     { "phaseq_c", &cspice_phaseq },
     { "pi_c",     &cspice_pi     },
     { "pipool_c", &cspice_pipool },
     { "pjelpl_s", &mice_pjelpl   },
     { "pl2nvc_c", &cspice_pl2nvc },
     { "pl2nvp_c", &cspice_pl2nvp },
     { "pl2psv_c", &cspice_pl2psv },
     { "psv2pl_c", &cspice_psv2pl },
     { "pxform_c", &cspice_pxform },
     { "pxfrm2_c", &cspice_pxfrm2 },
     { "q2m_c",    &cspice_q2m    },
     { "radrec_c", &cspice_radrec },
     { "rav2xf_c", &cspice_rav2xf },
     { "raxisa_c", &cspice_raxisa },
     { "reccyl_c", &cspice_reccyl },
     { "recgeo_c", &cspice_recgeo },
     { "reclat_c", &cspice_reclat },
     { "recpgr_c", &cspice_recpgr },
     { "recrad_c", &cspice_recrad },
     { "recsph_c", &cspice_recsph },
     { "rotate_c", &cspice_rotate },
     { "rotmat_c", &cspice_rotmat },
     { "rpd_c",    &cspice_rpd    },
     { "saelgv_c", &cspice_saelgv },
     { "scdecd_c", &cspice_scdecd },
     { "sce2c_c",  &cspice_sce2c  },
     { "sce2s_c",  &cspice_sce2s  },
     { "scencd_c", &cspice_scencd },
     { "scs2e_c",  &cspice_scs2e  },
     { "sct2e_c",  &cspice_sct2e  },
     { "sctiks_c", &cspice_sctiks },
     { "sincpt_s", &mice_sincpt   },
     { "spd_c",    &cspice_spd    },
     { "sphcyl_c", &cspice_sphcyl },
     { "sphlat_c", &cspice_sphlat },
     { "sphrec_c", &cspice_sphrec },
     { "spkcls_c", &cspice_spkcls },
     { "spkcov_c", &cspice_spkcov },
     { "spkcpo_s", &mice_spkcpo   },
     { "spkcpt_s", &mice_spkcpt   },
     { "spkcvo_s", &mice_spkcvo   },
     { "spkcvt_s", &mice_spkcvt   },
     { "spkezr_s", &mice_spkezr   },
     { "spkobj_c", &cspice_spkobj },
     { "spkopn_c", &cspice_spkopn },
     { "spkpos_s", &mice_spkpos   },
     { "spkpvn_s", &mice_spkpvn   },
     { "spksfs_s", &mice_spksfs   },
     { "spkw08_c", &cspice_spkw08 },
     { "srfxpt_s", &mice_srfxpt   },
     { "stpool_c", &cspice_stpool },
     { "str2et_c", &cspice_str2et },
     { "subpnt_s", &mice_subpnt   },
     { "subpt_s",  &mice_subpt    },
     { "subslr_s", &mice_subslr   },
     { "subsol_c", &cspice_subsol },
     { "surfnm_c", &cspice_surfnm },
     { "surfpt_s", &mice_surfpt   },
     { "sxform_c", &cspice_sxform },
     { "timdef_get_c", &cspice_timdef_get },
     { "timdef_set_c", &cspice_timdef_set },
     { "timout_c", &cspice_timout },
     { "tkvrsn_c", &cspice_tkvrsn },
     { "tpictr_c", &cspice_tpictr },
     { "tsetyr_c", &cspice_tsetyr },
     { "twopi_c",  &cspice_twopi  },
     { "twovec_c", &cspice_twovec },
     { "tyear_c",  &cspice_tyear  },
     { "unitim_c", &cspice_unitim },
     { "unload_c", &cspice_unload },
     { "unorm_c",  &cspice_unorm  },
     { "vdist_c",  &cspice_vdist  },
     { "vhat_c",   &cspice_vhat   },
     { "vnorm_c",  &cspice_vnorm  },
     { "vperp_c",  &cspice_vperp  },
     { "vprjp_c",  &cspice_vprjp  },
     { "vprjpi_c", &cspice_vprjpi },
     { "vproj_c",  &cspice_vproj  },
     { "vrotv_c",  &cspice_vrotv  },
     { "vsep_c",   &cspice_vsep   },
     { "wncard_c", &cspice_wncard },
     { "wncomd_c", &cspice_wncomd },
     { "wndifd_c", &cspice_wndifd },
     { "wnelmd_c", &cspice_wnelmd },
     { "wnexpd_c", &cspice_wnexpd },
     { "wnextd_c", &cspice_wnextd },
     { "wnfetd_c", &cspice_wnfetd },
     { "wnfild_c", &cspice_wnfild },
     { "wnfltd_c", &cspice_wnfltd },
     { "wnincd_c", &cspice_wnincd },
     { "wninsd_c", &cspice_wninsd },
     { "wnintd_c", &cspice_wnintd },
     { "wnreld_c", &cspice_wnreld },
     { "wnsumd_s", &mice_wnsumd   },
     { "wnunid_c", &cspice_wnunid },
     { "wnvald_c", &cspice_wnvald },
     { "xf2eul_c", &cspice_xf2eul },
     { "xf2rav_c", &cspice_xf2rav },
     { "xfmsta_c", &cspice_xfmsta },
     };

#endif


