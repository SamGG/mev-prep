# cytoPreparePercent // xtraPrepareDataset

cytoPreparePercent aims to transform, center and scale percents from cytometry analysis in order to use MeV (http://www.tm4.org/mev.html) for carrying out analysis.

##Introduction

Flow Cytometry and other medium or high throughput technics measure many variables at once in a single sample. An experiment consists in many samples representing a few conditions. In FC, variables are usually the percentage of each gated population in each sample. The percentage could be organised as a spreadsheet (aka data frame), in which lines represent populations, columns respresent samples, and the intersection is the observed percentage. Such a data set is multi-dimensionnal: each sample has a coordinate in the multi-dimensionnal space of populations. The challenge is to explore and analyse the experiment in a global fashion. Classical analysis based on making graph and statistics on a population by population basis is tedious. Heatmap is a synthetic representation that can be useful. An interactive implementation is available in MeV, allowing zooming, selection and clustering. To achieve a general representation, the data set needs to be transformed in order to retain the pertinent information of each population. The aim of this tool is ease the overall process that results in a data set ready for MeV for example.

A typical processing is to transform data via a logarithmic function and then to center. Data are then visualized and statistics (usually line by line) are computed. But visualisation don't represent the data the statistics use in reality. For example, a t-test computes a score in which the difference between two conditions is related to the pooled dispersion of the conditions. This score leads to computation of a p-value and the selection of the corresponding parameter if it is statistically significant. Sometimes a parameter is statistically selected although it does not show an interesting pattern in the visualisation. This sounds as the duality between fold change and p-value filtering when filtering genes out of a micro-array experiment. In fact the described processing leads to fold change whereas the t-test relates those fold changes to the dispersion within each group. So an interesting option is to apply a similar transformation that leads to a better agreement visualiation and statistics.


##Typical use

The percentages are organised as a spreadsheet such as MeV requirement.

The data are first transform in a logarithmic scale to cope the wide dynamic range.

Within each population, the data are then centered, allowing relative difference to be compared.
This leads to fold change data in logarithmic scale.

While fold changes are interesting per se, statistical computation and selection are based on the relative difference of FC to their dispersion within conditions. Scaling to the median within dispersion is proposed here. Resulting data could be compared between parameter because the resulting within dispersion is adjusted to 0.5. Others per-line scaling are also available: min-max, percentile and overall dispersion. Percentile scaling is a robust alternative to minmax, because it is less sensible to extreme values. 

If a reference conditions is available, a final centering allows


##Features

The processing is organized in successive steps that are described here.

###Transform

Percentage can be very different from population to population. Percentages cover a wide range of value, from 0.1% up to 50%. In classical boxplots the percentage are presented using a logarithimic scale, what allows compression while keeping relative differences. Usually, the interest is relative differences between conditions, not the percentage itself. Although logarithim transform is compressing the top of the range, it tends to expand the bottom of the range. This is especially important  because noise is usually more present at low percentage where percentages are obtained from a few hundreds of events (or may be less). To avoid to get high differences resulting from noise at low percentages, a constant could be added before. This features a soft low-end limit. An alternative transformation is the arcsinh, which behaves like a logarithm. The soft limit is obtained by dividing the percentage by a constant before applying the arcsinh.

###Center

###Scale

###Option post-centering to a reference

