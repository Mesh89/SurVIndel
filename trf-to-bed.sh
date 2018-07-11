awk '
BEGIN {
	seq = ""
}


{
	if ($0 ~ /Sequence:/) {
		seq = $2
	} else if (seq != "" && NF == 15) {
		print seq"\t"$1"\t"$2"\t"$3
	}
}'
