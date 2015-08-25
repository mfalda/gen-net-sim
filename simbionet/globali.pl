# a partire dagli oggetti cancellati e da il nome del file scrive le strutture per "globali.h", il construttore e il distruttore per "globali.c", le definizioni per i sinonimi da porre nel file (ALLA FINE), infine sostituisce tutte le variabili citate

$txt = << "END"
	CANCELLAv_i(S);
	CANCELLAv_i(P);
	CANCELLAv_i(ind);
	CANCELLAv_i(ind2);
	CANCELLAv_i(tmp1_i);
	CANCELLAv_i(Hio);
	CANCELLAv_i(H);
	CANCELLAm_i(pl);
END
;

$ns = "hubs";
local($/) = undef;
open(INPUT, $ns . ".c")
    or die "non riesco ad aprire il file '$filename' in lettura: $!\n";
$testo = <INPUT>;
close(INPUT);

@l = split(";", $txt);

print "\tstruct {\n";
foreach (@l) {
	if (m/.*CANCELLAv_(.)\((.*)\)/) {
		print "\t\tVETTORE$1 *$2;\n"
	}
	elsif (m/.*CANCELLAm_(.)\((.*)\)/) {
		print "\t\tMATRICE$1 *$2;\n"
	}
	elsif (m/.*CancellaLISTA\((.*), (.*)\)/) {
		print "\t\tLISTA *$1;\n"
	}
}
print "\t} $ns;\n\n";

foreach (@l) {
	if (m/.*CANCELLA._.\((.*)\)/ || m/.*CancellaLISTA\((.*), (.*)\)/) {
		print "\tglobali.$ns.$1 = NULL;\n";
	}
}
print "\n";

foreach (@l) {
	if (m/.*(CANCELLA._.)\((.*)\)/) {
		print "\t$1(globali.$ns.$2);\n";
	}
	elsif (m/.*(CancellaLISTA)\((.*), (.*)\)/) {
		print "\t$1(globali.$ns.$2, $3);\n";
	}
}
print "\n";

foreach (@l) {
	if (m/.*CANCELLA._.\((.*)\)/ || m/.*CancellaLISTA\((.*), .*\)/) {
		$nome = $1;
		print "#define g_$nome globali.$ns.$nome\n";
		$testo =~ s/\b$nome\b(.)/g_$nome$1/g;
	}
}

open OUTPUT, ">" . $ns . "1.c";
print OUTPUT $testo;
close(OUTPUT);