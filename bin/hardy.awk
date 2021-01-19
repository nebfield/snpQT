awk '{ if ($3=="TEST" || $3=="UNAFF" && $9 <0.0000001) print $0 }' \
	plink.hwe > plinkzoomhwe.hwe
