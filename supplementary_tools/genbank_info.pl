#!/usr/bin/perl
use strict;
use warnings;
#use 5.010;
#use autodie;
#use IO::Zlib;  #requires Compress::Zlib,

#variables
my $country = "";
my $pubmedId = "";
my $host = "";
my $isoSource = "";
my $geneNumber = "";
my $shape = "";

my $summary = "summary_GenBank_info.xls";


my $version = "1.0.1";

my $start = time();
my $command = "";

print "##################################################################\n";
print "# --> Welcome to $0 (version $version)!\n";
print "##################################################################\n";

#.gbff.gz
#Entete summary
open(SUM, ">", $summary) or die "Error writing file $!: ";
print SUM "File_Name\tCountry\tPUBMED_ID\tHost\tIsolation_Source\tGenes_Nb\tShape\n";

#FASTA/Q files
if(@ARGV){
  for my $arg (@ARGV){
    #if ($arg =~ m/.gbff/ or $arg =~ m/.gbff.gz/ ){  #
        if ($arg =~ m/.gz/ ){
          $command = "zcat $arg";
        }
        else{
          $command = "cat $arg";
        }
        print "GenBank Full format file: $arg\n";
        #my $genbankFile = $arg;
                my @array = split(/\./, $arg); #$array[0];

                open GBFF, "$command |";
                while (<GBFF>) {

                        chomp;
                        if ($_ =~ /\/country="(.*)"/) { $country = trim($1); }
                        if ($_ =~ /PUBMED(.*)/) {  $pubmedId = trim($1); }
                        if ($_ =~ /\/host="(.*)"/) {  $host = trim($1); }
                        if ($_ =~ /\/isolation_source="(.*)"/) {  $isoSource = trim($1); }
                        if ($_ =~ /\(Genes \(total\)\s+::(.*)/) { $geneNumber = trim($1); }
                        if ($_ =~ /LOCUS.*\s+([a-z]{1,})\s+[a-z]{1,}\s+[0-9]{2,}-[a-z]{1,}-[0-9]{4,}$/i) { $shape = trim($1); }

                }
                close(GBFF) or die "error close file $!:";

                print SUM "$array[0]\t$country\t$pubmedId\t$host\t$isoSource\t$geneNumber\t$shape\n";
                print "$array[0]\t$country\t$pubmedId\t$host\t$isoSource\t$geneNumber\t$shape\n";

    #}
  }
}

close(SUM) or die "error close file $!:";


my $end = time();

my $total = $end - $start;
my $min = $total / 60;
my $hrs = $min / 60;

print "\n\n";
print "***** Thank you for using $0! \n";
print "***** Total time: $total seconds OR $min minutes OR $hrs hours *****\n";


#-----------
# remove back and front spaces
sub trim {
my ($string) = @_;
$string =~ s/^\s+//;
$string =~ s/\s+$//;
return $string;
}
