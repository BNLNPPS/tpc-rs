TDataSet *CreateTable() { 
// -----------------------------------------------------------------
// bfc/.make/db/.const/StarDb/Geometry/tpc/.iTPCSurvey/iTPCSurvey Allocated rows: 24  Used rows: 24  Row size: 56 bytes
//  Table: iTPCSurvey_st[0]--> iTPCSurvey_st[23]
// ====================================================================
// ------  Test whether this table share library was loaded ------
  if (!TClass::GetClass("St_iTPCSurvey")) return 0;
iTPCSurvey_st row;
St_iTPCSurvey *tableSet = new St_iTPCSurvey("iTPCSurvey",24);
//
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =          1; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\x09P\x01\x00\x08P\x01\x00\x07P\x01\x00\x06P\x01\x00\x05P\x01\x00\x04P\x01\x00\x03P\x01\x00\x02P\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =          2; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xc6O\x01\x00\xc7O\x01\x00\xc8O\x01\x00\xc9O\x01\x00\xcaO\x01\x00\xcbO\x01\x00\xccO\x01\x00\xcdO\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =          3; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xd4O\x01\x00\xd5O\x01\x00\xd6O\x01\x00\xd7O\x01\x00\xd8O\x01\x00\xd9O\x01\x00\xdaO\x01\x00\xdbO\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =          4; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xe2O\x01\x00\xe3O\x01\x00\xe4O\x01\x00\xe5O\x01\x00\xe6O\x01\x00\xe7O\x01\x00\xe8O\x01\x00\xe9O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =          5; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xf0O\x01\x00\xf1O\x01\x00\xf2O\x01\x00\xf3O\x01\x00\xf4O\x01\x00\xf5O\x01\x00\xf6O\x01\x00\xf7O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =          6; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xfeO\x01\x00\xffO\x01\x00\x00P\x01\x00\xc0O\x01\x00\xbfO\x01\x00\xbeO\x01\x00\xbdO\x01\x00\xbcO\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =          7; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xb5O\x01\x00\xb4O\x01\x00\xb3O\x01\x00\xb2O\x01\x00\xb1O\x01\x00\xb0O\x01\x00\xafO\x01\x00\xaeO\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =          8; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xa7O\x01\x00\xa6O\x01\x00\xa5O\x01\x00\xa4O\x01\x00\xa3O\x01\x00\xa2O\x01\x00\xa1O\x01\x00\xa0O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =          9; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\x99O\x01\x00\x98O\x01\x00\x97O\x01\x00\x96O\x01\x00\x95O\x01\x00\x94O\x01\x00\x93O\x01\x00\x92O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         10; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\x00\x00\x00\x00\x00\x00\x00\x00\x89O\x01\x00\x88O\x01\x00\x87O\x01\x00\x86O\x01\x00\x85O\x01\x00\x84O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         11; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00DO\x01\x00EO\x01\x00FO\x01\x00GO\x01\x00HO\x01\x00IO\x01\x00JO\x01\x00KO\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         12; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00RO\x01\x00SO\x01\x00TO\x01\x00UO\x01\x00VO\x01\x00WO\x01\x00XO\x01\x00YO\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         13; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\x60O\x01\x00aO\x01\x00bO\x01\x00cO\x01\x00dO\x01\x00eO\x01\x00fO\x01\x00gO\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         14; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00nO\x01\x00oO\x01\x00pO\x01\x00qO\x01\x00rO\x01\x00sO\x01\x00tO\x01\x00uO\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         15; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00|O\x01\x00}O\x01\x00~O\x01\x00\x7fO\x01\x00\x80O\x01\x00@O\x01\x00?O\x01\x00>O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         16; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x007O\x01\x006O\x01\x005O\x01\x004O\x01\x003O\x01\x002O\x01\x001O\x01\x000O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         17; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00)O\x01\x00(O\x01\x00\x27O\x01\x00&O\x01\x00%O\x01\x00$O\x01\x00#O\x01\x00\x22O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         18; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\x1bO\x01\x00\x1aO\x01\x00\x19O\x01\x00\x18O\x01\x00\x17O\x01\x00\x16O\x01\x00\x15O\x01\x00\x14O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         19; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\x0dO\x01\x00\x0cO\x01\x00\x0bO\x01\x00\x0aO\x01\x00\x09O\x01\x00\x08O\x01\x00\x07O\x01\x00\x06O\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         20; // ;
    row.Angle	 =  -4.72e-05; // ;
    row.dx	 =    -0.0001; // ;
    row.dy	 =     0.0024; // ;
    row.ScaleX	 =   0.000326; // ;
    row.ScaleY	 =   0.000109; // ;
 memcpy(&row.comment,"SN006\x00N\x01\x00\xc4N\x01\x00\xc5N\x01\x00\xc6N\x01\x00\xc7N\x01\x00\xc8N\x01\x00\xc9N\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         21; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xd0N\x01\x00\xd1N\x01\x00\xd2N\x01\x00\xd3N\x01\x00\xd4N\x01\x00\xd5N\x01\x00\xd6N\x01\x00\xd7N\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         22; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xdeN\x01\x00\xdfN\x01\x00\xe0N\x01\x00\xe1N\x01\x00\xe2N\x01\x00\xe3N\x01\x00\xe4N\x01\x00\xe5N\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         23; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xecN\x01\x00\xedN\x01\x00\xeeN\x01\x00\xefN\x01\x00\xf0N\x01\x00\xf1N\x01\x00\xf2N\x01\x00\xf3N\x01",32);// 
tableSet->AddAt(&row);
memset(&row,0,tableSet->GetRowSize());
    row.Id	 =         24; // ;
    row.Angle	 =          0; // ;
    row.dx	 =          0; // ;
    row.dy	 =          0; // ;
    row.ScaleX	 =          0; // ;
    row.ScaleY	 =          0; // ;
 memcpy(&row.comment,"\x00\xfaN\x01\x00\xfbN\x01\x00\xfcN\x01\x00\xfdN\x01\x00\xfeN\x01\x00\xffN\x01\x00\x00O\x01\x00\xc0N\x01",32);// 
tableSet->AddAt(&row);
// ----------------- end of code ---------------
 return (TDataSet *)tableSet;
}
