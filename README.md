
PrecisionFDA Brain Biomarker Challenge
======================================

MSKCC MedPhys Team Submission
-----------------------------
Joseph O. Deasy,
James Mathews,
Alexandra Miller,
Jung Hun Oh,
Maryam Pouryahya,
Allen Tannenbaum,
Maria Thor


  1. **[Biomarkers](#Biomarkers)**
  1. **[Usage](#Usage)**
  2. **[Build the utilities](#Build)**
  3. **[Installation](#Installation)**
  4. **[Documentation](#Documentation)**

Biomarkers <a name="Biomarkers"></a>
----------

The sub-challenge 1 biomarkers we found to be especially predictive of survival outcome are the 12 genes:

```
TIMP1
ZNF32
ADM
EIF3E
CAVIN1
PLAT
HAUS4
MYADM
MYL9
RRNAD1
EMP3
UBR5
```

The sub-challenge 2 biomarkers we found to be especially predictive of survival outcome are the 8 chromosome loci:

```
19q13.13
9p24.1
1p13.3
9p21.3
1p36.21
7p11.2
10q23.32
1p34.3
```

The sub-challenge 3 biomarkers we found to be especially predictive of survival outcome are the 7 genes:

```
ASB3
CRTAC1
PCMTD1
RAB34
SH3BGRL3
ULBP2
EIF3E
```

Usage <a name="Usage"></a>
-----
First build and install the package by following the instructions below. Or use `source("bbutilities/R/utilities.R")`. Then make sure you've downloaded the data files from the [challenge website](https://precision.fda.gov/challenges/8).

In an R session:

```
source('challenge_models.R')
```

Build the utilities <a name="Build"></a>
-------------------
At the Unix command line:

```
R CMD INSTALL bbutilities
R CMD build bbutilities
R CMD check bbutilities_1.0.0.tar.gz
R CMD Rd2pdf bbutilities
```

You can also run the script `build.sh`.

Installation <a name="Installation"></a>
------------
In an R session:

```
install.packages("bbutilities_1.0.0.tar.gz", repos = NULL, type="source")
```

Documentation <a name="Documentation"></a>
-------------

See `bbutilities.pdf`. After installation, you can also use `help(...)` on the following functions:

```
assess_subtype_with_anova
assess_subtype_with_mglm
assess_with_glm
assess_with_ttest
build_and_report_final_model
build_glm
calculate_stats
cross_validation_glm
load_data
merge_rankings
merge_specific_rankings
pipeline
rank_features
write_model_parameters_to_file
write_subset
``` 

