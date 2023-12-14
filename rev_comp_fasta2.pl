#!/usr/bin/perl
use strict;
use warnings;

print "type in the path of the file\n";
my $file_name = <>;
chomp($file_name); 

open (FASTA, $file_name) or die "error #!"; 

$/ = ">";
<FASTA>;    
while (my $entry = <FASTA>){
    $entry = reverse $entry;
    $entry =~ tr/ACGTacgt/TGCAtgca/;
    print "$entry \n";
}

close(FASTA);
