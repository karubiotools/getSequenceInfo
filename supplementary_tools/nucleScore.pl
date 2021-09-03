#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
#use Shannon::Entropy qw/entropy/;
use File::Basename;
#use Bio::Species;

#use FindBin;
#use lib "$FindBin::RealBin/../perl5";

#my $input = $ARGV[0];
#chercher comment faire une liste perl pour input
#my @liste = split(/,/, $input);
#my $recap_total_seq = $ARGV[1];

#my ($input, $recap_total_seq) = @ARGV;

#my $start = time();

#my $file = ""; #= $ARGV[0];
#my $recap_total_seq = "nucleScore_result.xls";

#open (RECAP,'>', $recap_total_seq) or die "could not open $!";
  print "File\tA percent\tT percent\tC percent\tG percent\tGC percent\tAT/GC ratio\tNucleScore\tATG\tTGA\tTAG\tTAA\tGenome size (bp)\n";
#close(RECAP);


#FASTA files
#if(@ARGV){

  #for (my $i = 0; $i <= $#ARGV; $i++) {
    #if ($ARGV[$i]=~/-output/i or $ARGV[$i]=~/-o/i) {
    #           $recap_total_seq = $ARGV[$i+1];
    #}
  #}


 #open (RECAP,'>>', $recap_total_seq) or die "could not open $!";

#refaire le for pour la liste input
for my $arg (@ARGV){
#for my $arg (@liste){
#    if ($arg =~ m/.fasta/ or $arg =~ m/.fna/ or $arg =~ m/.fa/){

        #print "Traitement du fichier de sequence: $arg\n";
        #print "Traitement du fichier de sequence: $arg\n";
        #my $file = $arg;
        my $file = $arg;


        my $seqIO = Bio::SeqIO->new(-format=>'Fasta', -file=>$file);
        my $globalSeq = "";
        while (my $seq = $seqIO->next_seq()) {
                my $seqID = $seq->id;
                my $seqNuc = $seq->seq;
                $globalSeq .= $seqNuc;
                #push @arrayID, $seqID;
                #$hSeq{$seqID} = $seqNuc;
                #my @seqArray = split //, $seqNuc;
        }

        my $gcpercent = gc_percent($globalSeq);
        my ($ade, $thy, $gua, $cyt, $n, $length) = number_nuc_length_seq($file);
        my ($aPercent, $tPercent, $gPercent, $cPercent, $nPercent) = nucleotid_percent($ade, $thy, $gua, $cyt, $n, $length);

        my $atgcRatio =  atgc_ratio($ade, $thy, $gua, $cyt);

        my @percentList = ($aPercent, $tPercent, $gPercent, $cPercent, $nPercent);

        my $variance = shift_data_variance(@percentList);
        my $nucleScore = nucle_score($variance, $gcpercent, $atgcRatio, $length);
        #my $entropy = entropy($globalSeq);

        #print "The sequence length for $file is: $length\n";
        #print "A percent: $aPercent\n";
        #print "T percent: $tPercent\n";
        #print "G percent: $gPercent\n";
        #print "C percent: $cPercent\n";
        #print "N percent: $nPercent\n";

        #print "GC percent: $gcpercent\n";

        #print "AT/GC ratio: $atgcRatio\n";

        #print "NucleScore: $nucleScore\n";

        #print "Shannon Entropy: $entropy\n\n";

        #print "3 digits:\n";
        my @trinucs=($globalSeq=~/(?=(.{3}))/g);
        my %tri_count=();
        $tri_count{$_}++ for @trinucs;
        #print $_,":",$tri_count{$_},"\n" for sort keys(%tri_count);
        #print "\n2 digits:\n";
        my @trinucs2=($globalSeq=~/(?=(.{2}))/g);
        my %tri_count2=();
        $tri_count2{$_}++ for @trinucs2;
        #print $_,":",$tri_count2{$_},"\n" for sort keys(%tri_count2);

        my $atg = $tri_count{'ATG'};
        my $tga = $tri_count{'TGA'};
        my $tag = $tri_count{'TAG'};
        my $taa = $tri_count{'TAA'};

        #print "--------------------------------------\n\n";

        my $label = basename($file);


        #Summary file
        #print RECAP "$file\t$aPercent\t$tPercent\t$cPercent\t$gPercent\t$gcpercent\t$atgcRatio\t$nucleScore\t$entropy\t$aaa\t$aat\n";
        print "$label\t$aPercent\t$tPercent\t$cPercent\t$gPercent\t$gcpercent\t$atgcRatio\t$nucleScore\t$atg\t$tga\t$tag\t$taa\t$length\n";
        #}
  }
  #close (RECAP) or die "close file error : $!";
#}

#my $end = time();

#my $total = $end - $start;

#print "***** Total time (in seconds) is: $total *****\n";

#------------------------------------------------------------------------------
# number nucleotid and length
sub number_nuc_length_seq {
        my ($fastaFile) = @_;
        my $ade = 0;
        my $thy = 0;
        my $gua = 0;
        my $cyt = 0;
        my $n = 0;
        my $length = 0;

        open (FASTA, "<", $fastaFile) or die "Could not open $!";
        while (<FASTA>) {
                chomp;
                if ($_ !~ />/) {
                        my @seq = split //, $_;

                        for my $nuc (@seq) {
                                $length +=1 ;
                                if ($nuc =~ /a/i) {$ade+=1;}
                                elsif ($nuc =~ /t/i) {$thy+=1;}
                                elsif ($nuc =~ /g/i) {$gua+=1;}
                                elsif ($nuc =~ /c/i) {$cyt+=1;}
                                elsif ($nuc =~ /n/i) {$n+=1;}
                        }
                }
        }
        close(FASTA) or die "Error close file :$!";
        return ($ade, $thy, $gua, $cyt, $n, $length);

}

#------------------------------------------------------------------------------
# compute percentage of nucleotid
sub nucleotid_percent {
        my($ade, $thy, $gua, $cyt, $n, $length) = @_;

        my $adePercent = $ade / $length * 100;
        my $thyPercent = $thy / $length * 100;
        my $guaPercent = $gua / $length * 100;
        my $cytPercent = $cyt / $length * 100;
        my $nPercent = $n / $length * 100;

        return ($adePercent, $thyPercent, $guaPercent, $cytPercent, $nPercent);

}

#------------------------------------------------------------------------------
# compute GC pourcent
sub gc_percent {
        my ($seq) = @_;

        my @charSeq = split(//, uc($seq));
        my %hashFlank = ();

        foreach my $v (@charSeq) {
                $hashFlank{$v} += 1;
        }

        if (! $hashFlank{'G'}) { $hashFlank{'G'} = 0;}
        if (! $hashFlank{'C'}) { $hashFlank{'C'} = 0;}

        if(length($seq) == 0) {
                return 0;
        }
        else {
                return (($hashFlank{'G'} + $hashFlank{'C'}) / (length($seq))) * 100;
        }

}
#------------------------------------------------------------------------------
# compute ATGC ratio
sub atgc_ratio {
        my ($ade, $thy, $gua, $cyt) = @_;

        return (($ade + $thy) / ($gua + $cyt));

}
#------------------------------------------------------------------------------
# variance
sub shift_data_variance {
        my (@data) = @_;

        if ($#data + 1 < 2) { return 0.0; }

        my $K = $data[0];
        my ($n, $Ex, $Ex2) = 0.0;

        for my $x (@data) {
                $n = $n + 1;
                $Ex += $x - $K;
                $Ex2 += ($x - $K) * ($x - $K);
        }

        my $variance = ($Ex2 - ($Ex * $Ex) / $n) / ($n); ## ($n - 1)

        return $variance;

}
#------------------------------------------------------------------------------
# nucle score
#sub nucle_score {
#       my ($variance, $gcPercent, $atgcRatio, $length) = @_;
#
#       return (($variance * $gcPercent * $atgcRatio) / $length);
#}
sub nucle_score {
        my ($variance, $gcPercent, $atgcRatio, $length) = @_;
        return log2(($variance * $gcPercent * $atgcRatio ** (3)) / sqrt($length));
}

#------------------------------------------------------------------------------
sub log2 {
  my $n = shift;
  return (log($n) / log(2));
}
