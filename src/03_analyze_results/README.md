# Analyze PrediXcan predictions in SAGE

Compile summary statistics comparing gene expression predictions versus RNA-Seq measurements in SAGE.

## Prerequisites
The analysis requires the following R packages:
-- `data.table`
-- `purrr`
-- `broom`
-- `ggplot2`
-- `dplyr`
-- `optparse`

Install these with any standard technique, such as
```R
install.packages(c("data.table","purrr","broom","ggplot2","dplyr","optparse")
```

## Running
```bash
./analyze_predictions.sh
```
