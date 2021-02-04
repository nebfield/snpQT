#!/bin/awk -f

/loaded/ { print $1,$2 }
/Working directory/ { print FILENAME,$0 }
/remaining phenotypes/ { print FILENAME,$4,$6,$7,$8,$10 }
       
