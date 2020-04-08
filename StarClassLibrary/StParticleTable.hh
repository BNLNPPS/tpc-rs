/***************************************************************************
 *
 * $Id: StParticleTable.hh,v 1.5 2005/09/22 20:09:20 fisyak Exp $
 *
 * Author: Thomas Ullrich, May 99 (based on Geant4 code, see below)
 ***************************************************************************
 *
 * The design of the StParticleDefinition class and all concrete
 * classes derived from it is largely based on the design of the
 * G4ParticleDefinition class from Geant4 (RD44).
 * Although the code is in large parts different (modified or rewritten)
 * and adapted to the STAR framework the basic idea stays the same.
 *
 ***************************************************************************
 *
 * $Log: StParticleTable.hh,v $
 * Revision 1.5  2005/09/22 20:09:20  fisyak
 * Make StLorentzVector persistent
 *
 * Revision 1.4  2003/09/02 17:59:35  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.3  2000/01/04 14:57:56  ullrich
 * Added friend declaration to avoid warning messages
 * under Linux.
 *
 * Revision 1.2  1999/12/21 15:14:26  ullrich
 * Modified to cope with new compiler version on Sun (CC5.0).
 *
 * Revision 1.1  1999/05/14 18:50:00  ullrich
 * Initial Revision
 *
 **************************************************************************/
#ifndef StParticleTable_hh
#define StParticleTable_hh
#include "Rtypes.h"

#include <iostream>
#include <string>
#include <map>
#include <vector>

#include "tpcrs/logger.h"

class StParticleDefinition;
using std::vector;
using std::map;
using std::string;

typedef vector<StParticleDefinition*> StVecPtrParticleDefinition;

class StParticleTable
{
 public:
  virtual ~StParticleTable();

  static StParticleTable* particleTable();
  static StParticleTable* instance();

  unsigned int entries() const;
  unsigned int size() const;

  bool contains(const string &) const;              // by name
  bool contains(int) const;                         // by PDG encoding
  bool containsGeantId(int) const;                  // by Geant3 id

  StParticleDefinition* findParticle(const string &)  const; // by name
  StParticleDefinition* findParticle(int)  const;           // by PDG encoding
  StParticleDefinition* findParticleByGeantId(int) const;   // by Geant3 id

  void insert(StParticleDefinition*);
  void erase(StParticleDefinition*);

  void dump(std::ostream & = LOG_INFO);

  StVecPtrParticleDefinition allParticles() const;

  friend class nobody;

 private:
  StParticleTable();
  StParticleTable(const StParticleTable &right);

  static StParticleTable* mParticleTable;

  typedef map<int, int>                      mGeantPdgMapType;
  typedef map<int, StParticleDefinition*>    mPdgMapType;
  typedef map<string, StParticleDefinition*> mNameMapType;

  mGeantPdgMapType   mGeantPdgMap;     // Geant3 IDs only
  mPdgMapType        mPdgMap;          // PDG IDs only
  mNameMapType       mNameMap;         // complete list
};
#endif






