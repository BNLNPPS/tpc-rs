#pragma once

#include "TRandom.h"

#include "tpcrs/configurator.h"
#include "tpcrs/tpcrs_core.h"


namespace tpcrs { namespace detail {


class Digitizer
{
 public:

  Digitizer(const tpcrs::Configurator& cfg) :
    cfg_(cfg),
    digi_(cfg)
  {}

  template<typename InputIt, typename OutputIt>
  OutputIt Digitize(InputIt first_channel, InputIt last_channel, OutputIt digitized) const;

  template<typename InputIt, typename OutputIt>
  OutputIt Digitize(unsigned int sector, InputIt first_ch, InputIt last_ch, OutputIt digitized) const;

 private:

  void SimulateAltro(std::vector<short>::iterator first, std::vector<short>::iterator last, bool cancel_tail) const;
  void SimulateAsic(std::vector<short>& ADC) const;

  const tpcrs::Configurator& cfg_;
  tpcrs::DigiChannelMap digi_;
};


template<typename InputIt, typename OutputIt>
OutputIt Digitizer::Digitize(InputIt first_channel, InputIt last_channel, OutputIt digitized) const
{
  static DigiChannelMap digi(cfg_, 0);

  double pedRMS = cfg_.S<TpcResponseSimulator>().AveragePedestalRMSX;
  double ped = cfg_.S<TpcResponseSimulator>().AveragePedestal;

  auto ch_charge = first_channel;
  std::vector<short> ADCs_(digi.n_timebins, 0);
  std::vector<short> IDTs_(digi.n_timebins, 0);

  auto digitize = [ped, pedRMS](float charge, double gain) -> int
  {
    int adc = int(charge / gain + gRandom->Gaus(ped, pedRMS) - ped);
    // Zero negative values
    adc = adc & ~(adc >> 31);
    // Select minimum between adc and 1023, i.e. overflow at 1023
    adc = adc - !(((adc - 1023) >> 31) & 0x1) * (adc - 1023);
    return adc;
  };

  for (auto ch = digi.first(); !(digi.last() < ch); )
  {
    double gain = cfg_.S<tpcPadGainT0>().Gain[ch.sector-1][ch.row-1][ch.pad-1];

    // Go to the next pad
    if (gain <= 0) {
      digi.next_pad(ch);
      continue;
    }

    if (ch < ch_charge->channel || ch_charge == last_channel)
    { // digitize zero signal and continue
      ADCs_[ch.timebin-1] = digitize(0, gain);
      IDTs_[ch.timebin-1] = ch_charge->track_id;
    }
    else if (ch_charge->channel < ch)
    {
      // Skip channels with charges due to zero gain
      while(ch_charge->channel < ch)
        ++ch_charge;
    }
    else // equal channels
    {
      // digitize non-zero signal from ch_charge
      ADCs_[ch.timebin-1] = digitize(ch_charge->charge, gain);
      IDTs_[ch.timebin-1] = ch_charge->track_id;
      ++ch_charge;
    }

    // Pad boundary
    if (ch.timebin == digi.n_timebins)
    {
      SimulateAltro(std::begin(ADCs_), std::end(ADCs_), true);

      for (unsigned int tb = 1; tb != digi.n_timebins; ++tb)
      {
        if (ADCs_[tb-1] == 0) continue;
        *digitized = tpcrs::DigiHit{ch.sector, ch.row, ch.pad, tb, ADCs_[tb-1], IDTs_[tb-1]};
      }
    }

    digi.next(ch);
  }

  return digitized;
}


template<typename InputIt, typename OutputIt>
OutputIt Digitizer::Digitize(unsigned int sector, InputIt first_ch, InputIt last_ch, OutputIt digitized) const
{
  double pedRMS = cfg_.S<TpcResponseSimulator>().AveragePedestalRMSX;
  double ped = cfg_.S<TpcResponseSimulator>().AveragePedestal;

  std::vector<short> ADCs_(digi_.total_timebins(), 0);

  auto bc = first_ch;
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

  for (auto bc = first_ch; bc != last_ch; ++bc, ++ch, ++adcs_iter)
  {
    if (*adcs_iter == 0) continue;
    *digitized = tpcrs::DigiHit{sector, ch->row, ch->pad, ch->timebin, *adcs_iter, bc->track_id};
  }
}


} }
