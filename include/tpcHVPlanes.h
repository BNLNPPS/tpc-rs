/* tpcHVPlanes.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcHVPlanes.idl

  Table: tpcHVPlanes

       description: orientation of TPC HV planes:
                    central membrane
                    east gated grid
                    west gated grid

 */
#ifndef TPCHVPLANES_H
#define TPCHVPLANES_H
struct tpcHVPlanes_st {
	float CM_shift_z; /* physical z shift of the CM plane                      */
	float CM_tilt_x; /* x component of the CM plane's normal unit vector      */
	float CM_tilt_y; /* y component of the CM plane's normal unit vector      */
	float GGE_shift_z; /* physical z shift of the GG East plane                 */
	float GGE_tilt_x; /* x component of the GG East plane's normal unit vector */
	float GGE_tilt_y; /* y component of the GG East plane's normal unit vector */
	float GGW_shift_z; /* physical z shift of the GG West plane                 */
	float GGW_tilt_x; /* x component of the GG West plane's normal unit vector */
	float GGW_tilt_y; /* y component of the GG West plane's normal unit vector */
};
#endif /* TPCHVPLANES_H */
