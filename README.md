
PrecisionFDA Brain Biomarker Challenge
======================================

MSKCC MedPhys Team Submission
=============================

  1. **[Usage](#Usage)**
  2. **[Build the utilities](#Build)**
  3. **[Installation](#Installation)**
  4. **[Documentation](#Documentation)**

Usage <a name="Usage"></a>
-----
Make sure you've downloaded the data files from the [challenge website](https://precision.fda.gov/challenges/8).

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
build_glm
calculate_stats
create_folds_subtype
cross_validation
load_data
merge_rankings
merge_specific_rankings
pipeline
rank_features
write_subset
``` 




