#!/usr/bin/perl -w
use strict;

############################################################################
# script to remove position or column from a multi-Fasta file
# in function of a given character
############################################################################


my $inFile = $ARGV[0]; #'example_seq.fasta';
my $char = $ARGV[1]; #'N';
my @headers = ();
my @sequences = ();
my $index = 0;
my $outFile = 'results.fna';
open(IN,'<',$inFile) or die "Unable to read file $inFile: $!\n";
while( defined( my $line = <IN> ) ){
    chomp($line);
    if( $line =~ m/^>/ ){
        $headers[$index] = $line;
        $index++;
    }
    else{
        $sequences[$index-1] .= $line;
    }
}
close(IN);
my %lookup = ();
for(my $i=0;$i<=$#sequences;$i++){
    my $seq = $sequences[$i];
    my $len = length($seq);
        for(my $j=0;$j<$len;$j++){
        my $residue = substr($seq,$j,1);
        if( $residue eq $char ){
            $lookup{$j} = 1;
        }
    }
}
#print "# Skipped the following positions (zero indexed):\n";
#print "# ",join(", ", sort {$a <=> $b} keys (%lookup)), "\n";
#print "# Cleaned sequences:\n";
#open(OUT,'>',$outFile) or die "Unable to write file $outFile: $!\n";
for(my $i=0;$i<=$#headers;$i++){
    my $head = $headers[$i];
    my $seq = $sequences[$i];
    my $len = length($seq);
    my $out = '';
    for(my $j=0;$j<$len;$j++){
        my $residue = substr($seq,$j,1);
        $out .= $residue unless exists $lookup{$j};
    }
    print $head, "\n", $out, "\n";
    #print OUT $head, "\n", $out, "\n";
}
#close(OUT);
#print "\n";
#print "End of program! Your result is written in file $outFile\n";
