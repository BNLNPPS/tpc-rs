#ifndef St_tss_tssparC_h
#define St_tss_tssparC_h

#include "tpcrs/config_structs.h"
#include "tss_tsspar.h"
struct St_tss_tssparC : tpcrs::ConfigStruct<tpcrs::IConfigStruct, St_tss_tssparC, tss_tsspar_st>
{
  char* 	fileout(int i = 0) 	{return Struct(i)->fileout;}
  int 	dynam(int i = 0) 	{return Struct(i)->dynam;}
  int 	format(int i = 0) 	{return Struct(i)->format;}
  int 	max_itime(int i = 0) 	{return Struct(i)->max_itime;}
  int 	max_pads(int i = 0) 	{return Struct(i)->max_pads;}
  int 	max_row(int i = 0) 	{return Struct(i)->max_row;}
  int 	max_sect(int i = 0) 	{return Struct(i)->max_sect;}
  int 	min_itime(int i = 0) 	{return Struct(i)->min_itime;}
  int 	min_pads(int i = 0) 	{return Struct(i)->min_pads;}
  int 	min_row(int i = 0) 	{return Struct(i)->min_row;}
  int 	min_sect(int i = 0) 	{return Struct(i)->min_sect;}
  int 	mode(int i = 0) 	{return Struct(i)->mode;}
  int 	nele_laser(int i = 0) {return Struct(i)->nele_laser;}
  int 	ngain(int i = 0) 	{return Struct(i)->ngain;}
  int 	nseg(int i = 0) 	{return Struct(i)->nseg;}
  int 	ntime(int i = 0) 	{return Struct(i)->ntime;}
  int 	printout(int i = 0) 	{return Struct(i)->printout;}
  int 	tpc_half(int i = 0) 	{return Struct(i)->tpc_half;}
  int 	reset(int i = 0) 	{return Struct(i)->reset;}
  float 	ave_ion_pot(int i = 0) {return Struct(i)->ave_ion_pot;}
  float 	bfield(int i = 0) 	{return Struct(i)->bfield;}
  float 	c_test(int i = 0) 	{return Struct(i)->c_test;}
  float 	diff_long(int i = 0) 	{return Struct(i)->diff_long;}
  float 	diff_trans(int i = 0) {return Struct(i)->diff_trans;}
  float 	gain_in(int i = 0)    {return Struct(i)->gain_in;}
  float 	gain_in(int sec, int row) {return gain(sec, row);}
  float 	gain_out(int i = 0)   {return Struct(i)->gain_out;}
  float 	gain_out(int sec, int row)  {return gain(sec, row);}
  float 	gain(int sec, int row);
  float 	prf_in(int i = 0) 	{return Struct(i)->prf_in;}
  float 	prf_out(int i = 0) 	{return Struct(i)->prf_out;}
  float 	sca_rms(int i = 0) 	{return Struct(i)->sca_rms;}
  float 	scale(int i = 0) 	{return Struct(i)->scale;}
  float 	step_size(int i = 0) 	{return Struct(i)->step_size;}
  float 	tau(int i = 0) 	{return Struct(i)->tau;}
  float 	threshold(int i = 0) 	{return Struct(i)->threshold;}
  float 	time_offset(int i = 0) {return Struct(i)->time_offset;}
  float 	v_test(int i = 0) 	{return Struct(i)->v_test;}
  float 	white_rms(int i = 0) 	{return Struct(i)->white_rms;}
  float 	wire_coupling_in(int i = 0) 	{return Struct(i)->wire_coupling_in;}
  float 	wire_coupling_out(int i = 0) 	{return Struct(i)->wire_coupling_out;}
  float 	x_laser(int i = 0) 	{return Struct(i)->x_laser;}
  float 	y_laser(int i = 0) 	{return Struct(i)->y_laser;}
  float 	z_laser(int i = 0) 	{return Struct(i)->z_laser;}
};
#endif
