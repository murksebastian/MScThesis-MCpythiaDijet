# MScThesis-MCpythiaDijet
MC event generation for dijet and trijet event topologies produced at the LHC

##### dijet01.cc #####

- the main function takes two arguments (via user input): a random number seed and the number of events to generate; in order to run N events with RndmNSeed j, for instance, one could type

### make dijet01
### ./dijet01 j N > dijet01_j

in which case the PYTHIA readout (run statistics, etc.) will be saved to a file called "dijet01_j" (otherwise it will be printed in the terminal)
- running order of the strong coupling constant is set to 1
- by default, PYTHIA SlowJet algorithm is set to anti-kT algorithm with pTmin = 50 GeV
- kinematic phase space cuts: PhaseSpace:mHatMin = 2400, PhaseSpace:pTHatMin = 380
- kinematic selection cut on leading (subleading) jet is 440 (50) GeV
- angular selection cuts are |y^*|<1.7, |yB|<1.1
- the following jet parameters are stored: invariant mass, pT, y, y^*, yB, chi, phi (only used for sanity checks)
- events are grouped into seven mass intervals according to the dijet invariant mass (range is from 2.5 TeV < m < 5.4 TeV + m > 5.4 TeV)
- additional run statistics are printed; it is possible to check what fraction of events passes certain selection cuts, ends up in a particular mass interval, etc.
- histograms are output in the form of .txt files (two columns: 1) x value and 2) number of entries)

##### trijet01R10m3.cc #####

- works analogous to dijet01.cc, with an additional associated jet (ardQCD:all = off,HardQCD:3parton = on)
- additional kinematic phase space cut for the subsubleading jet: PhaseSpace:pTHat5Min = 114

##### Plotting macros; included in all of these is the option to correct the LO dijet topology up to NLO using k-factors derived from NLOJET++; they are extracted from the file "kfacOriginal.root"

##### diOkinVar.py #####

- allows to plot histograms of kinematic dijet variables for different running orders of the strong coupling constant using the .txt files produced with "dijet01.cc"

##### ditrikinVar.py #####

- allows to compare the kinematic variables of dijet and trijet events in a single histogram

##### chiHistRatioRvsW.py #####

- used to compared normalised angular distributions of the form 1/N dN/dchi vs. chi for weighted and unweighted trijet topologies together with the LO and NLO corrected dijet results

##### chiHistRatioWvsW.py #####

- allows to plot normalised angular distributions of the form 1/N dN/dchi vs. chi for two or more weighted trijet topologies generated with different minimum angular separation thresholds together with the LO and NLO corrected dijet results
