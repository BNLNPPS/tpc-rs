## Currently available tests for tpc-rs

Generated Interaction | Base Name          | Total Events | Reference
---                   | ---                | ---          | ---
dAu, 200 GeV          | `starY16_dAu200`   | 10           | STAR, Run 16, rcf16000_1_100evts.fzd
pp, 200 GeV           | `starY15_pp200a`   | 10           | STAR, Run 15, rcf15010_1_100evts.fzd
pp, 200 GeV           | `starY15_pp200b`   | 10           | STAR, Run 15, rcf15001_0_100evts.fzd
AuAu, 200 GeV         | `starY14_AuAu200a` | 10           | STAR, Run 14, rcf15000_1_100evts.fzd
AuAu, 200 GeV         | `starY14_AuAu200b` | 10           | STAR, Run 14, rcf14000_1_100evts.geant.root
He3Au, 200 GeV        | `starY14_He3Au200` | 10           | STAR, Run 14, rcf14010_1_100evts.fzd
CuAu, 200 GeV         | `starY12_CuAu200`  | 10           | STAR, Run 12, rcf12003_1_100evts.fzd
UU, 200 GeV           | `starY12_UU200`    | 10           | STAR, Run 12, rcf12002_1_100evts.fzd
pp, 500 Gev           | `starY12_pp500`    | 10           | STAR, Run 12, rcf12001_1_1000evts.fzd
pp, 500 GeV           | `starY11_pp500`    | 10           | STAR, Run 11, rcf10100_90_4000evts_minb.fzd
AuAu, 62 GeV          | `starY10_AuAu62`   | 10           | STAR, Run 10, rcf10033_1_100evts.fzd
AuAu, 11 GeV          | `starY10_AuAu11`   | 10           | STAR, Run 10, rcf10031_1_100evts.fzd

p = proton
d = deuteron
Au = gold
Cu = copper
He3 = helium-3
Reference: See https://www.star.bnl.gov/devcgi/weekDEVjobStatus.pl

## Data files

The archives with data files (ROOT and YAML) containing input and the respective
output values for the tests can be downloaded from a google drive.

Build Type | -mXX | Notes      | tar.gz Archive
Release    | 32   | no output? | [17N9nJ5iZxX-O2uVdGDpRHSkvwu32G4jh](https://drive.google.com/open?id=17N9nJ5iZxX-O2uVdGDpRHSkvwu32G4jh)
Debug      | 32   |            | [1XESdyqg6kXhu5gPrgS3drq3UIBu0f1l5](https://drive.google.com/open?id=1XESdyqg6kXhu5gPrgS3drq3UIBu0f1l5)
Debug      | 32   | <= v0.0.6  | [1Y2qPEN-sypcbWdrzFe4K7gGNCay5jpB8](https://drive.google.com/open?id=1Y2qPEN-sypcbWdrzFe4K7gGNCay5jpB8)
