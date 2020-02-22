#ifndef StMaker_h
#define StMaker_h

#include <string>

#include "TDatime.h"
#include "TTable.h"


class StChain
{
 public:

  TTable* GetDataBase(std::string path, const TDatime *td=0);

  /// Returns the current date. Should it mimick the date of the currently
  /// processed run?
  int GetDate() { return GetDateTime().GetDate(); }

  const TDatime GetDateTime() const { return TDatime(20151215, 0); }
};


class StMaker
{
 public:
  static StChain* gStChain;

  static StChain* GetChain() { return gStChain; }
};

#endif
