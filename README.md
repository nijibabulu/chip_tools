# chip_tools

This is a set of scripts used in the manuscript Analysis of SMAD1/5 target genes in a sea anemone reveals ZSWIM4-6 as a novel BMP signaling modulator. The steps are as follows:
1. run peakzilla on each chip/input pair (`peakzilla_qnorm_mapq_patched.py`). The version in this repository has a few patches from [the original](https://github.com/steinmann/peakzilla), but this does not affect default behavior of the script, excepting that fragment sizes from bedpe files are estimated from all the read pairs instead of only the first 100000. In the manuscript, the argument `--qnorm` was used.
1. join peaks that overlap between peakzilla runs into a single tsv file (`join_peaks.py`) using `--max-distance 100`.
1. associate the peaks with genes (`associate_peaks.py`).
1. filter the peaks (`final_filter.py`). 
