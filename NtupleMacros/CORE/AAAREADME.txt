
This directory contains only files that are _ABSOLUTELY_ needed by
_ALL_ users. Please limit content of the files to what is common for
main physics analyses. The main reason for that is supportability and
robustness of the common code. We have common code to avoid mistakes
in essential places. The more people use the code the higher is the
chance it will be properly debugged. Unfortunately if we drop here
everything that one can only think to be useful, we get complete mess,
which no one can really support.

Good examples:
  * Common selectors (no "SUSY" or "WW" selectors and stuff please)
    * electronSelections.* 
    * muonSelections.*
    * metSelectors.*
  * utilities
    * duplicate identification code

Examples of code that should be somewhere else:
  * beatifications like methods one-two lines long
  * trivial cuts, like Z-veto
  * analysis specific code - soft muon matching
  * stuff that in general seems to be fine just to be copied

Please keep number of files as small as possible.
