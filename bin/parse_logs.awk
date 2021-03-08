#!/bin/awk -f

BEGIN {
    plink1=0
    plink2=0
}

/PLINK v1.9/ { plink1=1 }
/PLINK v2.0/ { plink2=1 }

/Working directory/ { loaded["wd"] = $3 }

plink1 &&  /loaded from .bim/ { loaded["variants"] = $1 }
plink1 &&  /\) loaded from .fam/ { loaded["people"] = $1;
    loaded["males"] = substr($3, 2);
    loaded["females"] = $5;
    loaded["ambig"] = $7 }
plink1 &&  /remaining phenotypes/ { pheno["cases"] = $4;
    pheno["controls"] = $8 }
plink1 && /are controls.[[:blank:]]*\(/ { pheno["missing"] = $NF }

plink2 && /variants loaded from/ { loaded["variants"] = $1 }
plink2 && /samples \(/ { loaded["people"] = $1;
    loaded["males"] = $5;
    loaded["females"] = $3;
    loaded["ambig"] = $7 }
       
END {
    if (pheno["missing"] == "")
	pheno["missing"] = 0
    if (pheno["cases"] == "")
    {
	# if one is missing, they all are
	pheno["cases"] = 0
	pheno["controls"] = 0
	pheno["missing"] = 0
    }
    print FILENAME, loaded["variants"], loaded["people"],
	pheno["cases"] + pheno["controls"] + pheno["missing"],
	pheno["cases"], pheno["controls"], pheno["missing"],
	loaded["wd"]
}
