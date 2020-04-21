#include "Riostream.h"
#include "StBichsel/Bichsel.h"
#include <assert.h>
TString   Bichsel::m_Tags[kTotal] = {"P10", "Bi", "PAI"};
dEdxParameterization* Bichsel::m_dEdxParameterizations[kTotal] = {0, 0, 0};
Bichsel* Bichsel::fgBichsel = 0;


Bichsel::Bichsel(const Char_t* tag, Int_t keep3D) : m_Type(-1), m_Tag(tag), m_dEdxParameterization(0)
{

  for (Int_t k = 0; k < kTotal; k++) if (m_Tag.Contains(m_Tags[k].Data(), TString::kIgnoreCase)) {m_Type = k; break;}

  assert(m_Type >= 0);

  if (! m_dEdxParameterizations[m_Type])
    m_dEdxParameterizations[m_Type] = new dEdxParameterization(m_Tag.Data(), keep3D);

  m_dEdxParameterization = m_dEdxParameterizations[m_Type];
  fgBichsel = this;
}


Bichsel* Bichsel::Instance(const Char_t* tag, Int_t keep3D)
{
  if (!fgBichsel) new Bichsel(tag, keep3D);

  return fgBichsel;
}


void Bichsel::Clean()
{
  for (Int_t k = 0; k < kTotal; k++)
    if (m_dEdxParameterizations[k]) {
      delete m_dEdxParameterizations[k];
      m_dEdxParameterizations[k] = 0;
    }
}


void Bichsel::Print()
{
  LOG_INFO << "Bichsel:: " << m_Tag << '\n';

  if (m_dEdxParameterization) m_dEdxParameterization->Print();
}
// $Id: Bichsel.cxx,v 1.15 2015/12/24 00:16:25 fisyak Exp $
// $Log: Bichsel.cxx,v $
// Revision 1.15  2015/12/24 00:16:25  fisyak
// Add TpcRS model and macros
//
