/* trgTimeOffset.h */
/* This file was made by the idl compiler "stic". Do not edit.
** This was generated for version '(unspecified)'
** Instead, edit the source idl file and then re-run the compiler.
** For help, type contact Craig Tull or Herb Ward. */
/* COMMENTS FROM IDL FILE:
    trgTimeOffset.idl

  Table: trgTimeOffset

       description:   trigger offset common to all detector-level clock modules  trigger offset common to all detector-level clock modules


 */
#ifndef TRGTIMEOFFSET_H
#define TRGTIMEOFFSET_H
struct trgTimeOffset_st {
  float offset; /*   standard trigger offset in micro-seconds  */
  float laserOffset; /*   laser trigger offset in micro-seconds  */
  float laserOffsetW; /*   laser extra trigger offset for West laser */
};
#endif /* TRGTIMEOFFSET_H */
