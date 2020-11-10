## Variant lists as of 10/11/2020

* A combination of WES+WGS when there were samples with cvd2/cvd3/inf/neu measurements.

To find out which variants were exclusively from WES, we do
```
gunzip -c INTERVAL-wes+wgs.variantlist.gz | awk '$6<2000' | wc -l
```
giving 636,046 out of 116,133,865.

## Variant lists as of 09/11/2020

* The data contains both WES and WGS genotypes
  -- the sample overlap was only five in the four Olink panels.

* It is expected that there will be very minor changes with additional work on sample/genotype QC.
