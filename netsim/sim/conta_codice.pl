# estrae le righe senza i commenti

use strict;
use warnings;

my ($line, $dr, $f, $c, $tot, $comm);

if (scalar(@ARGV) < 1) {
   print "conta_codice.pl [-c] dir\n";
   print "-c: mostra le righe dei commenti e quelle vuote\n";
   exit;
}
if ($ARGV[0] eq "-c") {
  $comm = 1;
  $dr = $ARGV[1];
}
else {
  $comm = 0;
  $dr = $ARGV[0];
}
opendir(DIR, $dr) || die "non riesco a leggere la cartella '$dr!'";
CICLO:foreach $f (readdir(DIR)) {
	next CICLO if (-d $f || $f !~ m/^.*\.[cpp|h|d]/);
	open INPUT, "$f";
	print "File $f:\n\n";
	$tot = 0;
	while ($line = <INPUT>) {
		if ($line =~ m#^\s*$#) {
			print "\tlinea vuota in $.\n" if ($comm == 1);
		}
		if ($line =~ m#/\*#) {
			$c = $.; # linea corrente
		}
		elsif ($line =~ m#^\s*//#) {
			print "\tcommento in $.\n" if ($comm == 1);
		}
		elsif ($line =~ m#\*/#) {
			print "\tcommento tra $c e $.\n" if ($comm == 1);
			$c = 0;
		}
		else {
			$tot++;
		}
	}
	print "\nTotale: $tot linee di codice.\n\n\n";
	close INPUT;
}

