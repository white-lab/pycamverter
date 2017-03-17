# Change Log

## 0.2.0 (2017-03-17)

Features

  - Added support for using ProteomeDiscoverer search files, alongside MASCOT.
  - Started annotating export files with generator's pycamverter version.
  - Support csv scan lists.

Bugfixes

  - Fixed issue with peptideData being blank when ambiguous peptides were
    present.
  - Fixed issues with viewing peptides with similar types of modifications
  - Merge proteins from ambiguous peptide matches

## 0.1.2 (2017-03-10)

Features

  - Upped ProteoWizard version to 3.0.10577
  - Added C13 isotope labeling
  - Improved speed of spectra matching

## 0.1.1 (2017-03-07)

Features

  - Use python multiprocessing to cut data processing time in half.

## 0.1.0 (2017-03-06)

Features

  - Initial release, spun off from pyproteome.
