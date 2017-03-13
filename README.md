General
-------------------------
This directory contains project files for the following publication:

Tanger, P., Klassen, S., et al. 2017. Field-based high throughput phenotyping rapidly identifies genomic regions controlling yield components in rice. Scientific Reports. 7:42839

Comments and requests should be directed to:
paul.tanger14@alumni.colostate.edu

While all data and code is made freely available, we would appreciate citation when appropriate.


Files
-------------------------
There are 3 directories:

rilpopr_final - R scripts

input           - input files

output          - output files that scripts generate

Requirements
-------------------------
Obviously you need R installed (> 3.0.0 should do).

Initially set the working directory to rilpopr_final/

Based on wherever you unzipped the file from Datadryad.
```r
setwd("fill/in/your/path/scripts_tanger_etal_2017/")
```

If you downloaded from datadryad, you need to get the scripts.

These are hosted at bitbucket.org:

https://bitbucket.org/paultanger/scripts_tanger_etal_2017

You can either use git clone (if you have git installed),

or just download the files directly:
```r
# with git clone
setwd("../")
system("git clone https://paultanger@bitbucket.org/paultanger/scripts_tanger_etal_2017.git temp")
system("rm -r rilpopr_final")
system("mv temp rilpopr_final")
setwd("rilpopr_final/")
# then open the new readme.md

# or download directly
url = "https://bitbucket.org/paultanger/scripts_tanger_etal_2017/get/UPDATETHIS.zip"
thefile = tempfile(tmpdir=tempdir(), fileext=".zip")
download.file(url, thefile, method="curl")
unzip(thefile, junkpaths=T, exdir="../rilpopr_final", overwrite=T)
unlink(thefile)
```

If you downloaded from bitbucket, you need to get the data:
```r
url = "http://datadryad.org/bitstream/handle/10255/UPDATETHIS"
thefile = tempfile(tmpdir=tempdir(), fileext=".zip")
download.file(url, thefile, method="curl")
unzip(thefile, junkpaths=T, exdir="../input", overwrite=F)
unlink(thefile)
# and create output directory
system("mkdir ../output/")
```

The following packages are required:

This code only needs to be run once.

```r
install.packages("reshape2")
install.packages("Hmisc")
install.packages("lsmeans")
install.packages("lme4")
install.packages("plyr")
install.packages("ggplot2")
install.packages("grid")
install.packages("gridExtra")
install.packages("scales")
install.packages("corrplot")
install.packages("RColorBrewer")
install.pacakges("qtl")
install.packages("Parallel")
install.packages("doParallel")
```

Then, load these packages.  This needs to be done every time before running the scripts.

```r
library("reshape2")
library("Hmisc")
library("lsmeans")
library("lme4")
library("plyr")
library("ggplot2")
library("grid")
library("gridExtra")
library("scales")
library("corrplot")
library("RColorBrewer")
library("qtl")
library("parallel")
library("doParallel")
```

How to Use
-------------------------
The files used to process the data and generate figures, tables, and objects for further manipulation:

Feel free to explore each script individually to examine how the data was processed and modify as you see fit.

```r
source("RUNFIRST.R")
```

This loads libraries and a couple handy functions

```r
source("CombineDS2013Manual_HTP_Data.R")
```
This script takes the raw phenotyping data and creates two objects with plot means.  One object is for subsequent steps (DS2013AllMeansNewWithCTDFixHI.Robject) and one object is for running correlations later (AllmeansWideNew.Robject).

```r
source("2013FielDatafilter_and_format.R")
```
This script takes the raw field data for manually measured traits and creates a file for subsequent use (DS2013phenofinalFIXROWCOLfixfield_20160929_1943.tsv)

```r
source("RmodelingmarkdownDS2013alldataFinalLSmeansWithParents.R") 
```
This calculates LSmeans and variance components used in subsequent steps and Table S2.  Table S3 was obtained via JMP

```r
source("SetChrOrderToRef.R")
source("TableS4ALL_QTLs.R")
```
Produces Table S4.

```r
source("callfromIR64function.R")
source("geno_for_Rqtl_andJoinMap_v3FullPop.R")
```

Due to filesize, the raw GBS reads have been deposited in the NCBI Short Read Archive (SRA),
here: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA378634

The keyfile to process the files is stored in the metadata file associated with the reads, and also included in the datadryad data package.

tassel code final.txt contains the steps using tassel to process the reads.  The output (filtered.hmp.txt in the input directory) is used by these scripts in an iterative process importing and exporting from JoinMap to construct the final map (genotypematrixRestMapFinal_20150514_1239.csv).  Also used by these scripts is a function from another R package that is associated with a manuscript in preparation and hosted here: https://github.com/mckaylab/TSPmap.  The functions used are finddups() and removedups().  More information about how to install that package and how to call the functions can be found at the github site.

AfterGenMapFiltering
AfterGenMapProcessing
AfterGenMapProcessingAfterRemoveBadChr12
FinalRqtlMapFiltering
FinalFinalRqtlMapFiltering

These are useful functions used by scripts above:
fixOrderOfChrGenMap.R
FormatForRqtlFakecMWithXLineName.R
FormatForRqtlFakecM.R
FormatForRQTL.R
FormatForJoinMap.R
ExtractGenoFromCross.R
filtermarkersv1.R
ParallelDropOneMarker.R

```r
source("QTLv4server.R")
source("QTLscan2v1server.R")
```
This calculates the QTL and permutations

```r
source("QTLgetGenesFinalLSmeansPhenos.R")
```

This runs the model selection: FinalLSmeansstepwiseDFafterfixchr_20150624_1332.RObject
It also uses this: choose_best_model_function.R

```r
source("GraphDS2013ManualHTP_Data_FigS2.R")
```
Produces Figure S2

```r
source("FigS3code.R")
```
Produces Figure S3

```r
source("Fig2.R")
```
Produces Figure 2

```r
source("X.R")
```
Produces Figure 3.. produced with JMP based on Lovell's code

```r
source("NewFig4.R")
```
Produces Figure 4


License
-------------------------
Copyright (c) 2017, Paul Tanger, Brook Moyers, John Lovell

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.