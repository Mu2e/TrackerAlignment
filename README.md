# Tracker Alignment with Millepede-II

Ryunosuke O'Neil
([@ryuwd](https://github.com/ryuwd))

A Tracker Alignment utility for performing Track-based alignment with the [Millepede-II software package](https://www.desy.de/~kleinwrt/MP2/doc/html/index.html). Only no-field cosmic tracks are currently supported.

My thesis on this work: [FERMILAB-MASTERS-2020-07](https://inspirehep.net/literature/1849484).

### Alignment set-up
![Alignment flow(2)](https://user-images.githubusercontent.com/56410978/82936768-fa2e6500-9f86-11ea-81fe-b9f0bf20e842.png)
- 'Digis' refers to art files containing 'digi' data products. Current configuration is for cosmics in the extracted position (no field).
- Track Reco uses Richie's Time Fit ( [Mu2e-doc-33162](https://mu2e-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=33162) )
- 'MILLE' means to write input files for the 'PEDE' executable. The `MILLE` stage is implemented by `MilleDataWriter`.
- 'Align Tracker' refers to the stage carried out by `AlignedTrackerMaker` (TrackerConditions). It uses alignment constants provided to the Proditions interface to change tracker straw positions accordingly.
- 'Track Collection' is where the `AlignTrackCollector` module selects tracks, calculates the needed quantities for alignment, and writes the data to file using `MilleDataWriter`.
- 'PEDE' refers to the millepede executable that performs the alignment fit given the output file(s) produced by the `MilleDataWriter` class during the `AlignTrackCollector` job(s).


### MC Alignment test
1. Set up your backing build
```
mkdir BackingBuild
cd BackingBuild
muse backing HEAD
```
2. Build TrackerAlignment
```
git clone git@github.com:mu2e/TrackerAlignment
muse setup
muse build -j 4
```
3. Set up a TrackerAlignment helper environment.
```
source TrackerAlignment/setup.sh
```

2. Choose a working directory and setup your starting misalignment
```
cd /exp/mu2e/data/users/$USER/
mkdir alignment-test
cd alignment-test

mu2ealign new MDC2020_startup_planes MDC2020_best v1_3

```
This command generates job.fcl, alignconstants_in.txt, revision.txt, and sources.txt. This is effectively a working directory for one alignment iteration. This also saves the Mu2e Offline revision at the time of generation. The file env.sh is also created which when sourced will add the current directory to your MU2E_SEARCH_PATH so that alignment files placed here are found by Offline.

Here, the Tracker plane misalignment is overwritten with 'MDC2020_startup' configuration. The corresponding DbService text file in `${TRKALIGN_BASE}/test/misalignments/` is copied. Then the remaining conditions are set to use purpose 'MDC2020_best' version 'v1_3'

Note that by default the weak modes are constrained to zero. You may also choose to fix specific planes instead, or set an uncertainty on the constraint by changing
```
physics.analyzers.AlignTrackCollector.WeakConstraints : "None/Fix/Measurement"
physics.analyzers.AlignTrackCollector.FixPlane : [ 5, 30 ]

```

3. Run Track Collection + PEDE interactively

You can run one process to collect the track data like this:
```bash
# running one job only
mu2ealign run
```

Or, you may run multiple processes:
```bash
# running multiple jobs e.g. 4 jobs processing 8 input files, 
# or 8 / 4 = 2 input art files per job
mu2ealign_genparallel 4 8
mu2ealign run
```

Once the jobs finish, you can run the first alignment fit:
```bash
mu2ealign pede
```
Now you should make a new directory, and import the produced alignment constants:
```bash
mkdir iter1 && cd iter1
mu2ealign new ../alignconstants_out.txt
```
Then return to the beginning of this step.

Or, you can run multiple iterations automatically (easier, recommended)
```bash
# 
mu2ealign autorun 3 # for 3 alignment iterations and one nominal run
```

# Issues
- Only the track drift time measurement is included in the chi^2 which is why chi^2/Ndof tends to peak at 0.5, rather than 1 when you plot it. Also Millepede will report a global chi2/ndof of around 0.5.
   - In the actual cosmic fit, there are two terms added to the chi-squared for the track per straw hit (see the definition [here](https://github.com/Mu2e/Offline/blob/ff2d1d20467d56c67c3c035784de95d0df47f490/CosmicReco/src/PDFFit.cc#L317-L332))


# Alignment Derivatives

## Terms
### Global degrees of freedom
These are alignment parameters that currently describe a translation or rotation in x, y, z applied to one Plane, and (later) Panel. These are 'global' because their values are common across all input tracks.

### Local degrees of freedom
These are the track fit parameters for one track. One set of local parameters describes one track.

## Derivatives?

For each global and local degree of freedom, Millepede requires a partial derivative of the track residual with respect to that degree of freedom.

There can be as many as 100 or 1000 global degrees of freedom, however of those Millepede only needs to know what the non-zero derivatives are. In this case, only those DOFs corresponding to the Plane and Panel holding a crossed (or hit) Straw will be nonzero. This sparse storage scheme saves a lot of space.

For the cosmic track time fit, there is a numerical differentiation utility to calculate these quantities. 
