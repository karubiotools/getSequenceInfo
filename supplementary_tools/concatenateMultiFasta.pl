#!/usr/bin/perl -w
use strict;

################################################################
# Script allowing to concatenate multiFasta file,
# generating an output file containing a single sequence
################################################################

# example of use: perl concatenateMultiFasta.pl multiFasta_file.fasta
# other example:  perl concatenateMultiFasta.pl *.fasta

my @listFastaFiles = @ARGV;

foreach my $multiFasta ( @listFastaFiles ) {
        my $outFasta = 'concatenated_'.$multiFasta ;
        open(FILE,"<$multiFasta") || die ("Error opening $multiFasta $!");
        open(OUT, '>', $outFasta) or die $!;
        print OUT ">concatenation $multiFasta\n";
        while (my $row = <FILE>) {
                chomp $row;
                if ($row=~m/^>/){
                }
                else{
                    print OUT "$row\n";
                }
        }
}

close(FILE);
close(OUT);
