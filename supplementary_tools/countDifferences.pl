#!/usr/bin/perl -w
use strict;
use warnings;
#use DateTime;
#use Date::Calc qw(:all);
use Bio::SeqIO;
#use String::Similarity;
#use StringOO::Similarity;
#use lib qw(./StringOO);
#use Similarity::similarity;

#use local::lib "$FindBin::Bin/..";  ### points to ~/mydir1 and local::lib finds lib
#use lib "$FindBin::StringOO/../lib";
#package Similarity.pm;
#use String::Approx 'adistr';
#package String::Similarity;
#package String::Similarity qw(similarity);
#use String::Similarity;

##################################################################
## Perl program to count variants from aligned multi-FASTA file
##################################################################

my $file = $ARGV[0];
my %hSeq = ();
my @arrayID = ();
my $resultFile = $ARGV[1]; #"resultVariants.xls";
my $resultSimilarity = $ARGV[2]; #"resultSimilarity.xls";
my $chain = "";

        my $seqIO = Bio::SeqIO->new(-format=>'Fasta', -file=>$file);

        while (my $seq = $seqIO->next_seq()) {
                my $seqID = $seq->id;
                my $seqNuc = $seq->seq;
                push @arrayID, $seqID;
                $hSeq{$seqID} = $seqNuc;
                #my @seqArray = split //, $seqNuc;
        }

        open(RES, ">", $resultFile) or die "error open file $!:";
        open(SIM, ">", $resultSimilarity) or die "error open file $!:";

        print RES "\t";
        print SIM "\t";
        foreach my $id (@arrayID) {
                print RES "$id\t";
                print SIM "$id\t";
        }
        print RES "\n";
        print SIM "\n";

        for(my $i = 0; $i <= $#arrayID; $i++){
                print RES "$arrayID[$i]\t";
                print SIM "$arrayID[$i]\t";
                for(my $j = 0; $j <= $#arrayID; $j++){  #for(my $j = $i+1; $j <= $#arrayID; $j++){
                        my $snpVariants = compare($hSeq{$arrayID[$i]}, $hSeq{$arrayID[$j]});
                        #my $similarity = adistr($hSeq{$arrayID[$i]}, $hSeq{$arrayID[$j]});
                        my $similarity = ident($hSeq{$arrayID[$i]}, $hSeq{$arrayID[$j]}); # replace similarity by ident
                        my $roundSim = sprintf('%.3f', $similarity);
                        print RES "$snpVariants\t";
                        print SIM "$roundSim\t";
                        #print RES ("SNPs between  $arrayID[$i] and $arrayID[$j] = $snpVariants\t");
                }
                print RES "\n";
                print SIM "\n";
        }

close(RES) or die "error when closing file $!:";
close(SIM) or die "error when closing file $!:";

        #foreach my $id (@arrayID) {
        #       my @seqArray = split //, $hSeq{$id};
        #
        #}


#------------------------------------------------------------------------------
# compare two aligned sequences
sub compare {
        my ($seq1,$seq2) = @_;
        #$seq1 = uc $seq1;
        #$seq2 = uc $seq2;
        my @seqArray1 = split //, $seq1;
        my @seqArray2 = split //, $seq2;
        my $variants = 0;

        for(my $i = 0; $i <= $#seqArray1; $i++){
                if( $seqArray1[$i] ne $seqArray2[$i] ){
                        $variants++;
                }
        }

        return $variants;
}

# calculate percent identity between two aligned sequences
sub ident {
        my ($seq1,$seq2) = @_;
        #$seq1 = uc $seq1;
        #$seq2 = uc $seq2;
        my @seqArray1 = split //, $seq1;
        my @seqArray2 = split //, $seq2;
        my $sim = 0;

        for(my $i = 0; $i <= $#seqArray1; $i++){
                if( $seqArray1[$i] eq $seqArray2[$i] ){
                        $sim++;
                }
        }

        return ($sim / ($#seqArray1 + 1));
}

#------------------------------------------------------------------------------
