#ifndef St_db_Maker_h
#define St_db_Maker_h

class TDatime;
class TTable;

class St_db_Maker
{
 public:

  static Int_t GetValidity(const TTable* tb, TDatime* const val);
};

#endif
