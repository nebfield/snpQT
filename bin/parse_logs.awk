#!/bin/awk -f

/Working directory/ { print FILENAME,$0 }
/loaded/ { print $1,$2 }
/remaining phenotypes/ { print $4,$6,$7,$8,$10 }
/are controls.[[:blank:]]*\(/ { print $NF, "phenotypes are missing)"
       
