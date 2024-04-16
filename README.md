# dengue-singapore

This repo contains code and data to support Finch _et al_ "Climate variation and serotype competition drive dengue outbreak dynamics in Singapore"(in press). This study quantifies the role of climate variation and serotype competition in shaping dengue risk in Sinapore using over 20 years of weekly case data, and integrates these findings into an early warning system framework able to predict dengue outbreaks up to 2 months ahead. Model fitting is performed in `INLA`. Climate and dengue case data are available in the `data` folder. Note that data on dengue serotype frequencies shown in Figure 1 are not included in this repo as they are not publicly available.

To reproduce the analysis, open the .Rproj file. The following scripts read in the data and run the analysis:
- [R/00_load-data.R](https://github.com/EmilieFinch/dengue-singapore/blob/main/R/00_load-data.R)
- [R/01_run-selected-models.R](https://github.com/EmilieFinch/dengue-singapore/blob/main/R/01_run-selected-models.R)
- [R/02_run-tscv-predictions](https://github.com/EmilieFinch/dengue-singapore/blob/main/R/02_run-tscv-predictions.R)

The following scripts contain functions needed for analysis. These do not need to be run seperately as they are sourced in the main workflow scripts.
- [R/create-lagged-data_fn.R](https://github.com/EmilieFinch/dengue-singapore/blob/main/R/create-lagged-data_fn.R)
- [R/fit-inla_fn.R](https://github.com/EmilieFinch/dengue-singapore/blob/main/R/fit-inla_fn.R)
- [R/tscv-prediction_fn.R](https://github.com/EmilieFinch/dengue-singapore/blob/main/R/tscv-prediction_fn.R)
- [R/utils_fn.R](https://github.com/EmilieFinch/dengue-singapore/blob/main/R/utils_fn.R)

Model outputs are included in the `output` folder and manuscript figures can be generated by knitting [reports/manuscript-figures.Rmd](https://github.com/EmilieFinch/dengue-singapore/blob/main/reports/manuscript-figures.Rmd).

This analysis was performed using R 4.1.1 and INLA_21.11.22.
