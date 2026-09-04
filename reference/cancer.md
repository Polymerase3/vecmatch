# Patients with Colorectal Cancer and Adenoma

This is a synthetically generated dataset containing metadata for
healthy individuals and patients diagnosed with colorectal cancer or
adenomas. The primary purpose of this dataset in the context of matching
is to balance the `status` groups across various covariates and achieve
optimal matching quality.

## Usage

``` r
data(cancer)
```

## Format

A data frame (`cancer`) with 1,224 rows and 5 columns:

- status:

  Patient's health status, which can be one of the following: `healthy`,
  `adenoma`, `crc_benign` (benign colorectal carcinoma), or
  `crc_malignant` (malignant colorectal carcinoma).

- sex:

  Patient's biological sex, recorded as either `M` (male) or `F`
  (female).

- age:

  Patient's age, represented as a continuous numeric variable.

- bmi:

  Patient's Body Mass Index (BMI), represented as a continuous numeric
  variable.

- smoker:

  Smoking status of the patient, recorded as `yes` or `no`.

## Source

Data generated artificially
