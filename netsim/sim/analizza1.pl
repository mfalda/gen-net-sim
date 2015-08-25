# funziona solo con i vettori, le matrici (poche) vanno analizzate a mano, per ora
#~ $ARGV[0] = "memoria_v_i.csv";
open INPUT, $ARGV[0];
open OUTPUT, ">$ARGV[1]";

while (<INPUT>) {
	if (m/"\d+";"([^;]*)";[^;]*;[^;]*;\"([^;]+)\";[^;]*;"(.*) # (\d+)";"(.*) # (\d+)"/) {
		push(@vars, [$1, $3, $4, $5, $6, $2]);
	}
}

print OUTPUT "\"File\";\"Nome1\";\"Dim.max.\";\"Da1\";\"A1\";\"Nome2\";\"Dim.max.\";\"Da2\";\"A2\"\n";

for ($i = 0; $i < $#vars; $i++) {
	for ($j = $i + 1; $j < $#vars; $j++) {
		if (($vars[$i][1] ne "" && $vars[$j][1] ne "") && ($vars[$i][1] eq $vars[$i][3]) && ($vars[$j][1] eq $vars[$j][3]) && ($vars[$i][1] eq $vars[$j][1]) && ($vars[$i][4] < $vars[$j][2] || $vars[$j][4] < $vars[$i][2])) {
			#~ print "In " . $vars[$i][1] . ":\n";
			#~ print "\tvariabile" . $vars[$i][0] . ": " . $vars[$i][2] . " -> " . $vars[$i][4] . "\n";
			#~ print "\tvariabile" . $vars[$j][0] . ": " . $vars[$j][2] . " -> " . $vars[$j][4] . "\n";
			print OUTPUT "\"$vars[$i][1]\";\"$vars[$i][0]\";\"$vars[$i][5]\";\"$vars[$i][2]\";\"$vars[$i][4]\";\"$vars[$j][0]\";\"$vars[$j][5]\";\"$vars[$j][2]\";\"$vars[$j][4]\"\n";
		}
	}
}

close OUTPUT;
close INPUT;