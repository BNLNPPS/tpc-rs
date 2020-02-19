/* tpcChargeEvent.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    tpcChargeEvent.idl

  Table: tpcChargeEvent

       description: Table of events depositing significant charge into the TPC;
                    bunchCrossing is stored as 32bit pairs because no 64bit integer
                    works with the database, i.e. the first bunch crissing is actually
                    (eventBunchCrossingsHigh[0] << 32)  eventBunchCrossingsLow[0]
                    with appropriate 32to64 conversion before bitshifting

 */
#ifndef TPCCHARGEEVENT_H
#define TPCCHARGEEVENT_H
struct tpcChargeEvent_st {
	int nChargeEvents; /* number of charge events in this record */
	unsigned int eventBunchCrossingsLow[4096]; /* number of bunches into the run when charge event occurred (32 low bits) */
	unsigned int eventBunchCrossingsHigh[4096]; /* number of bunches into the run when charge event occurred (32 high bits) */
	float eventCharges[4096]; /* metric of magnitude of charge deposited in TPC */
	int badBunch; /* collider bunch where most of the charge events occurred */
};
#endif /* TPCCHARGEEVENT_H */
