use strict;
use warnings;

my ($i, $uno, $commento, $file, $paren, $paren1, @livelli, $livello);

$file = $ARGV[0];
#~ $file = "connectivity_modular.c";
open INPUT, $file;
open OUTPUT, ">$file.tsv";
$paren = 0;
$paren1 = 0;
$uno = 0;
$commento = 0;
$i = 0;
while (<INPUT>) {
	$i++;
  next if (m#^\s*//# || $commento == 1);
  if (m#.*\*/#) {
    $commento = 0;
    next;
  }
	if ($uno == 1) {
    $livello = pop(@livelli);
    print join(", ", @{$livello}), ", $i\n";
    print OUTPUT @{$livello}[0], "\t", @{$livello}[3], "\t", @{$livello}[2] . "\n";
    print OUTPUT "\t", $i, "\t", @{$livello}[2] . "\n\n";
		$paren--;
    $paren1--;
    $uno = 0;
	}
	if (m/^\s*(for|while|do)\b/) {
    $paren1++;
    push(@livelli, ["$1", $paren + 1, $paren1, $i]);
		if (m/{/) {
			$paren++;
			$uno = 0;
		}
		else {
			$uno = 1;
		}
	}
	elsif (m/{/) {
		$paren++;
	}
	elsif (m/}/) {
    $paren--;
    if ($#livelli >= 0 && $paren == $livelli[-1][1] - 1) {
      $livello = pop(@livelli);
      print join(", ", @{$livello}), ", $i\n";
      print OUTPUT @{$livello}[0], "\t", @{$livello}[3], "\t", @{$livello}[2] . "\n";
      print OUTPUT "\t", $i, "\t", @{$livello}[2] . "\n\n";
      $paren1--;
    }
	}
}

close OUTPUT;
close INPUT;
