#!/usr/bin/perl
print "Content-type: text/html\n\n";
open(HITS,"xtalcompcounter.txt");
$hits = <HITS>;
close(HITS);
++$hits;
open(HITS,">xtalcompcounter.txt");
print HITS $hits;
close(HITS);
print "$hits";
exit;
