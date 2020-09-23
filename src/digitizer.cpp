#include "tpcrs/detail/digitizer.h"

#include "altro.h"

namespace tpcrs { namespace detail {


void Digitizer::SimulateAltro(std::vector<short>::iterator first, std::vector<short>::iterator last, bool cancel_tail) const
{
  Altro altro_sim(last - first, &*first);

  if (cancel_tail) {
    //        ConfigAltro(ONBaselineCorrection1, ONTailcancellation, ONBaselineCorrection2, ONClipping, ONZerosuppression)
    altro_sim.ConfigAltro(                    0,                  1,                     0,          1,                 1);
    altro_sim.ConfigTailCancellationFilter(cfg_.S<tpcAltroParams>().Altro_K1,
                                           cfg_.S<tpcAltroParams>().Altro_K2,
                                           cfg_.S<tpcAltroParams>().Altro_K3,
                                           cfg_.S<tpcAltroParams>().Altro_L1,
                                           cfg_.S<tpcAltroParams>().Altro_L2,
                                           cfg_.S<tpcAltroParams>().Altro_L3);
  }
  else {
    altro_sim.ConfigAltro(0, 0, 0, 1, 1);
  }

  altro_sim.ConfigZerosuppression(cfg_.S<tpcAltroParams>().Altro_thr, cfg_.S<tpcAltroParams>().Altro_seq, 0, 0);
  altro_sim.RunEmulation();

  int i = 0;
  for (auto it = first; it != last; ++it, ++i) {
    if (*it && !altro_sim.ADCkeep[i]) { *it = 0;}
  }
}


void Digitizer::SimulateAsic(std::vector<short>& ADC) const
{
  int t1 = 0;
  int nSeqLo = 0;
  int nSeqHi = 0;

  for (unsigned int tb = 0; tb != ADC.size(); ++tb) {
    if (ADC[tb] <= cfg_.S<asic_thresholds>().thresh_lo) {
      if (! t1) ADC[tb] = 0;
      else {
        if (nSeqLo <= cfg_.S<asic_thresholds>().n_seq_lo ||
            nSeqHi <= cfg_.S<asic_thresholds>().n_seq_hi)
        {
          for (unsigned int t = t1; t <= tb; t++) ADC[t] = 0;
        }
      }

      t1 = nSeqLo = nSeqHi = 0;
    }

    nSeqLo++;

    if (!t1) t1 = tb;
    if (ADC[tb] > cfg_.S<asic_thresholds>().thresh_hi) {nSeqHi++;}
  }
}

} }
