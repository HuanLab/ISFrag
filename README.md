ISFrag R Package User Manual
================
Sam Shen, Jian Guo, Tao Huan
24/03/2021

-   [Part 1: Introduction and
    Installation](#part-1-introduction-and-installation)
-   [Part 2: MS1 Feature Extraction](#part-2-ms1-feature-extraction)
    -   [2.1 XCMS Feature Extraction](#21-xcms-feature-extraction)
    -   [2.2 Additional Featuretable
        Input](#22-additional-featuretable-input)
-   [Part 3: MS2 Annotation](#part-3-ms2-annotation)
-   [Part 4: In-source Fragments
    Identification](#part-4-in-source-fragments-identification)
-   [Part 5: Results Export](#part-5-results-export)
    -   [5.1 Export ISF Result
        Featuretable](#51-export-isf-result-featuretable)
    -   [5.2 Export ISF Relationship
        Tree](#52-export-isf-relationship-tree)
-   [Part 6: Additional Details and
    Notes](#part-6-additional-details-and-notes)

# Part 1: Introduction and Installation

`ISFrag` is an R package for identifying and annotating in-source
fragments in LCMS metabolite featuretable. The package is written in the
language R and its source code is publicly available at
<https://github.com/HuanLab/ISFrag.git>.

To install `ISFrag` package R version 4.0.0 or above is required, and we
recommend using RStudio to complete the installation and usage of
`ISFrag` by following the steps below:

``` r
# Install "devtools" package from CRAN if you do not already have it installed.
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}

# Load "devtools" package.
library(devtools)

# Install "ISFrag" from Github using "devtools".
if (!requireNamespace("ISFrag", quietly = TRUE)){
  install_github("HuanLab/ISFrag")
}

# Load "ISFrag" package.
library(ISFrag)
```

# Part 2: MS1 Feature Extraction

`ISFrag` supports multiple ways to generate an MS1 featuretable. Users
can choose to use `XCMS` to extract features from mzXML files, upload
their own featuretable in csv format, or combine both features extracted
by `XCMS` with their own featuretable. For the rest of the tutorial, to
view details and additional parameters of functions, type:
`help("<function name>")`.

## 2.1 XCMS Feature Extraction

One or multiple mzXML files from DDA, DIA, or fullscan analyses can be
analyzed at once using XCMS to extract MS1 features. All mzXML file(s)
need to be placed in a separate folder containing no other irrelevant
mzXML files. Note: for multi-sample analyses, peak alignment and filling
will be performed by XCMS. Additional details of XCMS is available at:
<https://rdrr.io/bioc/xcms/man/>.

``` r
# MS1directory specifies the full directory of the folder containing mzXML file(s).
MS1directory <- "X:/Users/Sam_Shen/ISFtest20210127/RP(-)/RP(-)1/fullscan"

# The generate.featuretable() function outputs a dataframe formatted featuretable as well as an MSnbase object.
featureTable <- generate.featuretable(MS1directory = MS1directory, type = "single", peakwidth = c(5,20))
head(featureTable)
```

    ##          mz      rt   rtmin   rtmax maxo
    ## F1 1123.664 862.081 860.071 864.590 2226
    ## F2 1122.684 862.081 860.071 864.087 3096
    ## F3 1121.678 862.081 859.568 864.590 3908
    ## F4 1107.725 974.010 972.502 975.516 1728
    ## F5 1107.662 862.081 860.574 864.087 2966
    ## F6 1105.740 974.511 971.999 976.521 2180

## 2.2 Additional Featuretable Input

To use a custom featuretable for `ISFrag` analysis or combine unique
features from both the `XCMS` and custom featuretable, use the
`add.features()` function. The `featuretable` parameter is not needed if
users wish to only use their own custom featuretable for analyses. In
order for `ISFrag` to succesfully read the provided csv file, it must
contain only columns in the following order: m/z, retention time, min
retention time, max retention time, followed by an additional column
containing the intensities of features detected in each sample. Note:
column 3 and column 4 are the retention time of the feature edges, and
all three columns containing retention time information should be in
seconds.

``` r
# ft_directory specifies the directory of the custom csv file.
ft_directory <- "X:/Users/Sam_Shen/ISFtest20210127/RP(-)"
# ft_name specifies the name of the custom csv file.
ft_name_single <- "NISTplasmaDDARP(-)1featuretable.csv"
ft_name_multi <- "NISTplasmaDDARP(-)Alignedfeaturetable.csv"

# Sample csv feature table for single file analysis.
setwd(ft_directory)
head(read.csv(ft_name_single, header = T, stringsAsFactors = F))
```

    ##          mz     rt     rtmin    rtmax Intensity
    ## 1 194.90540 25.444 19.459002 35.47100  5737.125
    ## 2  85.02976 29.384  7.153002 70.32402  1118.750
    ## 3  87.00925 29.384  7.531002 70.32402  2168.750
    ## 4 111.00890 29.384  7.153002 72.63300  4777.000
    ## 5 173.00910 29.384  4.920000 75.61200  1909.250
    ## 6 173.00920 29.384  5.726000 64.40298  1909.250

``` r
# Sample csv feature table for a 3-file analysis, not that there are 3 columns containing feature intensity from each sample.
head(read.csv(ft_name_multi, header = T, stringsAsFactors = F))
```

    ##         mz      rt   rtmin   rtmax NISTplasmaDDARP...1_P1.B.1_01_12331
    ## 1 44.99868 1272.12 1272.12 1272.12                               80584
    ## 2 44.99871 1442.40 1442.40 1442.40                               82195
    ## 3 44.99881 1681.08 1681.08 1681.08                               43457
    ## 4 56.99625 1688.40 1688.40 1688.40                                4393
    ## 5 56.99628   52.74   52.74   52.74                                5108
    ## 6 59.01415 1222.14 1222.14 1222.14                              821808
    ##   NISTplasmaDDARP...2_P1.B.1_01_12332 NISTplasmaDDARP...3_P1.B.1_01_12333
    ## 1                               76202                               77925
    ## 2                               65772                               82181
    ## 3                               39078                               49169
    ## 4                                4225                                4342
    ## 5                                5882                                5352
    ## 6                              838424                              850325

``` r
# The add.features() function outputs a dataframe formatted featuretable containing
featureTable <- add.features(featureTable = featureTable, ft_directory = ft_directory, ft_name = ft_name_single)
head(featureTable)
```

    ##          mz      rt   rtmin   rtmax maxo
    ## F1 1123.664 862.081 860.071 864.590 2226
    ## F2 1122.684 862.081 860.071 864.087 3096
    ## F3 1121.678 862.081 859.568 864.590 3908
    ## F4 1107.725 974.010 972.502 975.516 1728
    ## F5 1107.662 862.081 860.574 864.087 2966
    ## F6 1105.740 974.511 971.999 976.521 2180

# Part 3: MS2 Annotation

One or multiple mzXML files from DDA analyses are needed to assign MS2
spectrum to features and perform annotation. The number of mzXML file(s)
provided here does not need to correspond with the number of mzXML files
used in the feature extraction step earlier. All mzXML file(s) need to
be placed in a separate folder containing no other irrelevant mzXML
files. In addition, the standard library used to perform annotation must
be in msp format.

``` r
# MS2directory specifies the full directory of the folder containing DDA mzXML file(s).
MS2directory <- "X:/Users/Sam_Shen/ISFtest20210127/RP(-)/RP(-)1/DDA"

# The ms2.tofeaturetable() function assigns MS2 spectra from the provided DDA files to the MS1 featuretable. It returns a new featuretable with additional columns containing MS2 fragment information.
featureTable <- ms2.tofeaturetable(MS2directory = MS2directory, featureTable = featureTable)
head(featureTable)
```

    ##          mz      rt   rtmin   rtmax maxo MS2_match MS2mz MS2int PeaksCount
    ## F1 1123.664 862.081 860.071 864.590 2226     FALSE     0      0          0
    ## F2 1122.684 862.081 860.071 864.087 3096     FALSE     0      0          0
    ## F3 1121.678 862.081 859.568 864.590 3908     FALSE     0      0          0
    ## F4 1107.725 974.010 972.502 975.516 1728     FALSE     0      0          0
    ## F5 1107.662 862.081 860.574 864.087 2966     FALSE     0      0          0
    ## F6 1105.740 974.511 971.999 976.521 2180     FALSE     0      0          0
    ##    fromFile ISF_level
    ## F1        0         0
    ## F2        0         0
    ## F3        0         0
    ## F4        0         0
    ## F5        0         0
    ## F6        0         0

``` r
# Now, use the feature.annotation() function to annotate features in the featuretable against a standard database in msp format.This functions returns a featuretable containing additional columns with annotation information.
lib_directory <- "X:/Users/Sam_Shen/Library" # directory containing the library file
lib_name <- "MoNA-export-LC-MS-MS_Negative_Mode.msp" # name of the library file
featureTable <- feature.annotation(featureTable = featureTable, lib_directory = lib_directory, lib_name = lib_name, dp = 0.1)
```

    ## [1] "Using cores:"
    ## [1] 15

``` r
head(featureTable)
```

    ##          mz      rt   rtmin   rtmax maxo MS2_match MS2mz MS2int PeaksCount
    ## F1 1123.664 862.081 860.071 864.590 2226     FALSE     0      0          0
    ## F2 1122.684 862.081 860.071 864.087 3096     FALSE     0      0          0
    ## F3 1121.678 862.081 859.568 864.590 3908     FALSE     0      0          0
    ## F4 1107.725 974.010 972.502 975.516 1728     FALSE     0      0          0
    ## F5 1107.662 862.081 860.574 864.087 2966     FALSE     0      0          0
    ## F6 1105.740 974.511 971.999 976.521 2180     FALSE     0      0          0
    ##    fromFile ISF_level Annotation DPscore
    ## F1        0         0    unknown       0
    ## F2        0         0    unknown       0
    ## F3        0         0    unknown       0
    ## F4        0         0    unknown       0
    ## F5        0         0    unknown       0
    ## F6        0         0    unknown       0

# Part 4: In-source Fragments Identification

In-source fragments are identified from Level 3 to Level 1 fragments
through functions `find.leve3()`, `find.level2()`, `find.level1()`,
respectively. Each function returns a list of featuretable, where each
featuretable contains a parent feature in the first row, with remaining
rows containing candidate in-source fragment features. Note: these
functions must be used in order level 3, 2, 1.

The `find.level3()` function takes in the `MS1directory` string,
`MS1files` string vector outputted by the `generate.featuretable()`
function, and the analysis `type` (single or multi). The `find.level2()`
functions takes in the output of `find.level3()`, and similarly the
`find.level1()` function takes in the output of `find.level2()`.

``` r
# Identify level 3 in-source fragments.
level3 <- find.level3(MS1directory = MS1directory, MS1.files = MS1.files, featureTable = featureTable, type = "single")
```

    ## [1] "Using cores:"
    ## [1] 15

``` r
# Identify level 2 in-source fragments.
level2 <- find.level2(ISFtable = level3)
```

    ## [1] "Using cores:"
    ## [1] 15

``` r
# Identify level 1 in-source fragments.
level1 <- find.level1(ISF_putative = level2)
```

    ## [1] "Using cores:"
    ## [1] 15

# Part 5: Results Export

Once all level 3, 2, and 1 in-source fragments are identified,
`run get.ISFrag.results()` function to summarize the analysis results.
This step must be done in order to export either ISF relationship tree
or featuretable.

``` r
# Summarize ISFrag results after identifying all level 3, 2, and 1 in-source fragments.
results <- get.ISFrag.results(ISF_List = level1, featureTable = featureTable)
```

## 5.1 Export ISF Result Featuretable

Either the complete featuretable with ISF relationship annotated in
additional columns, or a detailed ISF parent-fragment relationship
dataframe for a single parent feature can be exported.

``` r
# Get complete featuretable with all features and ISF relationship annotations.
resultFT <- export.ISFrag.results(ISFresult = results)
head(resultFT)
```

    ##          mz      rt   rtmin   rtmax maxo MS2_match MS2mz MS2int PeaksCount
    ## F1 1123.664 862.081 860.071 864.590 2226     FALSE     0      0          0
    ## F2 1122.684 862.081 860.071 864.087 3096     FALSE     0      0          0
    ## F3 1121.678 862.081 859.568 864.590 3908     FALSE     0      0          0
    ## F4 1107.725 974.010 972.502 975.516 1728     FALSE     0      0          0
    ## F5 1107.662 862.081 860.574 864.087 2966     FALSE     0      0          0
    ## F6 1105.740 974.511 971.999 976.521 2180     FALSE     0      0          0
    ##    fromFile ISF_level Annotation DPscore Num_Level2 Num_Level1
    ## F1        0         0    unknown       0          0          0
    ## F2        0         0    unknown       0          0          0
    ## F3        0         0    unknown       0          0          0
    ## F4        0         0    unknown       0          0          0
    ## F5        0         0    unknown       0          0          0
    ## F6        0         0    unknown       0          0          0

``` r
# Get detailed featuretable for a specified parent feature (eg. feature F1).
detailedResults <- export.ISFrag.detailed(ISF_List = level1, featureID = "F1")
# Here the first row is the parent feature F1, and the remaining rows are its in-source fragment features.
head(detailedResults)
```

    ##              mz      rt   rtmin   rtmax maxo MS2_match MS2mz MS2int PeaksCount
    ## F1    1123.6641 862.081 860.071 864.590 2226     FALSE     0      0          0
    ## F10   1100.6641 862.081 860.071 864.087 2186     FALSE     0      0          0
    ## F18   1062.6617 862.081 860.574 863.586 1550     FALSE     0      0          0
    ## F591   614.3136 862.081 859.568 864.087 1912     FALSE     0      0          0
    ## F1380  460.2588 862.081 860.574 863.586 1972     FALSE     0      0          0
    ##       fromFile ISF_level Annotation DPscore     ppcor
    ## F1           0         0    unknown       0 0.0000000
    ## F10          0   Level_3    unknown       0 0.9837755
    ## F18          0   Level_3    unknown       0 0.9716118
    ## F591         0   Level_3    unknown       0 0.8958322
    ## F1380        0   Level_3    unknown       0 0.9549163

## 5.2 Export ISF Relationship Tree

ISF relationship trees show the hierarchical relationship of the parent
feature and its in-source fragment features. Tree diagrams for all
parent features can be exported at once, or the tree diagram for a
specified parent feature can be exported. When using
`plot.tree.single()` to draw a single tree diagram, the provided feature
ID must be that of a parent feature which contains either level 2 or 1
in-source fragment features.

``` r
# Specify the directory the tree diagrams should be plotted to.
output_dir <- "X:/Users/Sam_Shen/ISFtest20210127/RP(-)"

# Plot tree diagrams for all parent features.
plot.tree.all(ISFresult = results, directory = output_dir)

# Plot tree diagram for a single specified parent feature (eg. feature F9).
plot.tree.single(ISFresult = results, featureID = "F9", directory = output_dir)
```

<div id="htmlwidget-629922aebe872b0a26a2" style="width:672px;height:480px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-629922aebe872b0a26a2">{"x":{"diagram":"digraph {\n\ngraph [rankdir = \"TB\"]\n\n\n\n  \"1\" [label = \"F9: mz=1101.71, rt=908, ParentFeature\", style = \"filled,rounded\", shape = \"box\", tooltip = \"- name: F9: mz=1101.71, rt=908, ParentFeature\", fillcolor = \"#63B8FF\", fontcolor = \"#FFFFFF\"] \n  \"2\" [label = \"F1374: mz=461.27, rt=908, PPcor=0.97, level=Level_2\", style = \"filled,rounded\", shape = \"box\", tooltip = \"- name: F1374: mz=461.27, rt=908, PPcor=0.97, level=Level_2\", fillcolor = \"#FFA54F\", fontcolor = \"#FFFFFF\"] \n\"1\"->\"2\" [arrowhead = \"vee\", color = \"grey35\", penwidth = \"2\"] \n}","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script>

# Part 6: Additional Details and Notes
