open INPUT, $ARGV[0];

while (<INPUT>) {
	if (m/"\d+";"([^;]*)";"(.*) # \d+";[^;]*;"(\d+)";"(\d+)";.*/ && $4 > 0) {
		print "In $2: $1 riallocata $4 volte/i (dim. max. $3)\n";
	}
}

close INPUT;