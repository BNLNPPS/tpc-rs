#ifndef tpcrs_DigiData_h
#define tpcrs_DigiData_h

#include <vector>


namespace tpcrs {

/**
 * Associates a packed hardware index with the data.
 *
 *                sector   row    pad   timebin
 * projected max      24    72    182       512
 * bits                5 +   7 +   10 +      10 = 32 = 4 bytes = unsigned int
 * max value          32   128   1024      1024 
 */
struct DigiChannel
{
  unsigned int sector : 5, row : 7, pad : 10, timebin : 10;
  short adc;
  short idt;
};


class DigiData
{
 public:

  const std::vector<DigiChannel>& channels() const { return channels_; }

  void Add(unsigned int sector, unsigned int row, unsigned int pad, short* ADCs, short* IDTs, int n_timebins)
  {
    bool in_cluster = false;

    for (unsigned int tb = 0; tb < n_timebins; ++tb)
    {
      if (!ADCs[tb])
        in_cluster = false;

      if (ADCs[tb] && !in_cluster)
        in_cluster = true;

      if (in_cluster)
        channels_.push_back(DigiChannel{sector, row, pad, tb, ADCs[tb], IDTs[tb]});
    }
  }

 private:

  std::vector<DigiChannel> channels_;
};

}

#endif
