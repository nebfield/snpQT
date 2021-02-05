#!/bin/awk -f

/Working directory/ { print FILENAME,$0 }
/loaded from .bim/ { loaded["variants"] = $1 }
/\) loaded from .fam/ { loaded["people"] = $1;
    loaded["males"] = substr($3, 2);
    loaded["females"] = $5;
    loaded["ambig"] = $7 }
/remaining phenotypes/ { pheno["cases"] = $4;
    pheno["controls"] = $8 }
/are controls.[[:blank:]]*\(/ { pheno["missing"] = $NF }
       
END {
    if (pheno["missing"] == "")
	pheno["missing"] = 0
    print loaded["variants"], "variants"
    print loaded["people"], "samples"
    print pheno["cases"] + pheno["controls"] + pheno["missing"], "phenotypes",
        "(" pheno["cases"], "cases,", pheno["controls"], "controls,",
	"and", pheno["missing"], "are missing)"
}
