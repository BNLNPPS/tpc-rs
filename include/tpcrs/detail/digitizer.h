#pragma once

#include "TRandom.h"

#include "tpcrs/configurator.h"
#include "tpcrs/tpcrs_core.h"


namespace tpcrs { namespace detail {

using ChargeContainer = std::vector<tpcrs::SimulatedCharge>;

class Digitizer
{
 public:

  Digitizer(const tpcrs::Configurator& cfg) :
    cfg_(cfg),
    digi_(cfg)
  {}

  template<typename OutputIt>
  void Digitize(unsigned int sector, const ChargeContainer& binned_charge, OutputIt digitized) const;

 private:

  void SimulateAltro(std::vector<short>::iterator first, std::vector<short>::iterator last, bool cancel_tail) const;
  void SimulateAsic(std::vector<short>& ADC) const;

  const tpcrs::Configurator& cfg_;
  tpcrs::DigiChannelMap digi_;
};


template<typename OutputIt>
void Digitizer::Digitize(unsigned int sector, const ChargeContainer& binned_charge, OutputIt digitized) const
{
  double pedRMS = cfg_.S<TpcResponseSimulator>().AveragePedestalRMSX;
  double ped = cfg_.S<TpcResponseSimulator>().AveragePedestal;

  std::vector<short> ADCs_(binned_charge.size(), 0);

  auto bc = binned_charge.begin();
  auto adcs_iter = ADCs_.begin();

  for (auto ch = digi_.channels.begin(); ch != digi_.channels.end(); ch += digi_.n_timebins)
  {
    double gain = cfg_.S<tpcPadGainT0>().Gain[sector-1][ch->row-1][ch->pad-1];

    if (gain <= 0) {
      bc        += digi_.n_timebins;
      adcs_iter += digi_.n_timebins;
      continue;
    }

    for (int i=0; i != digi_.n_timebins; ++i, ++bc, ++adcs_iter)
    {
      int adc = int(bc->charge / gain + gRandom->Gaus(ped, pedRMS) - ped);
      // Zero negative values
      adc = adc & ~(adc >> 31);
      // Select minimum between adc and 1023, i.e. overflow at 1023
      adc = adc - !(((adc - 1023) >> 31) & 0x1) * (adc - 1023);

      *adcs_iter = adc;
    }
  }

  for (auto adcs_iter = ADCs_.begin(); adcs_iter != ADCs_.end(); adcs_iter += digi_.n_timebins)
  {
    SimulateAltro(adcs_iter, adcs_iter + digi_.n_timebins, true);
  }

  auto ch = digi_.channels.begin();
  adcs_iter = ADCs_.begin();

  for (auto bc = binned_charge.begin(); bc != binned_charge.end(); ++bc, ++ch, ++adcs_iter)
  {
    if (*adcs_iter == 0) continue;
    *digitized = tpcrs::DigiHit{sector, ch->row, ch->pad, ch->timebin, *adcs_iter, bc->track_id};
  }
}


} }
