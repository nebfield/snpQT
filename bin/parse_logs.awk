#!/bin/awk -f

/Working directory/ { loaded["wd"] = $3 }
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
    
    print FILENAME, loaded["variants"], loaded["people"],
	pheno["cases"] + pheno["controls"] + pheno["missing"],
	pheno["cases"], pheno["controls"], pheno["missing"],
	loaded["wd"]
}
