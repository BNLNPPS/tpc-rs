/* richvoltages.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    richvoltages.idl

  Table: richvoltages

       description:

 */
#ifndef RICHVOLTAGES_H
#define RICHVOLTAGES_H
#define RICHVOLTAGES_SPEC \
"struct richvoltages { \
	unsigned long runNumber; \
	unsigned long startStatusTime; \
	unsigned long endStatusTime; \
	unsigned long status; \
};"
typedef struct richvoltages_st {
  unsigned int runNumber; /*   */
  unsigned int startStatusTime; /*   */
  unsigned int endStatusTime; /*   */
  unsigned int status; /*   */
} RICHVOLTAGES_ST;
#endif /* RICHVOLTAGES_H */
