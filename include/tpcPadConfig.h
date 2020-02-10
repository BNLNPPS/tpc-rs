/* tpcPadConfig.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcPadConfig.idl

  Table: tpcPadConfig

       description: TPC padrow configurations (tpcitpc) for 24 sectors

 sector status 0 => tpc, 1 => itpc
 */
#ifndef TPCPADCONFIG_H
#define TPCPADCONFIG_H
#define TPCPADCONFIG_SPEC \
"struct tpcPadConfig { \
	octet itpc[24]; \
};"
typedef struct tpcPadConfig_st {
	unsigned char itpc[24]; 
} TPCPADCONFIG_ST;
#endif /* TPCPADCONFIG_H */
