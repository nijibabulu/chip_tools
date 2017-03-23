# chip_tools

This is a set of tools that I use for a peakzilla-based ChIP analysis pipeline (in development). The steps are as follows:
1. run peakzilla on each chip/input pair (`peakzilla_qnorm_mapq_patched.py`). The version in this repository has a few patches, but this does not affect default behavior of the script, excepting that fragment sizes from bedpe files are estimated from all the read pairs instead of only the first 100000.
1. join peaks that overlap between peakzilla runs into a single tsv file (`join_peaks.py`)
1. associate the peaks with genes (`associate_peaks.py`).
1. filter the peaks (`filter_peaks.py`). This step is most specific to the the user, as peak quality differs quite a bit between ChIP datasets. it's recommended to come up with your own filter and perhaps write an independent script (espceially given the patchwork code in this script).
