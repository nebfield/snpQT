#!/usr/bin/env bash

plink --bfile plink_8 --test-missing

awk '{ if ($5 < 10e-5) print $1, $2 }' plink.missing > fail_missingness.txt

plink --bfile plink_8 --remove fail_missingness.txt --make-bed --out plink_9


