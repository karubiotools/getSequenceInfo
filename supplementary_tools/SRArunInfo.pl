#!/usr/bin/perl

use strict;
use warnings;

my $start = time();

################################################################
# Script allowing to get SRA run information using accession ID
# wget command is required to use this script
################################################################

# example of use: perl SRArunInfo.pl SRR7693877,SRR9850824,SRR9850830 OR perl SRArunInfo.pl list_accessions.txt

# options


my $runs = $ARGV[0];

my @tabRuns = ();

if($ARGV[0] =~  m/,/ ){
  @tabRuns = split (/,/, $runs);
}
elsif(-e $runs){ 
  open my $handle, '<', $runs;
  chomp(@tabRuns = <$handle>);
  close $handle;
}
else{
  push(@tabRuns, $runs);
}

#@tabRuns = split (/,/, $runs) ;

my $summary = "summary_Runs.tsv";
my $country2 = "";
my @tabCSV = ();
my %hashCenter = ();

# Center names (e.g. abbreviation SC means: "The Wellcome Trust Sanger Institute")
$hashCenter{"SC"} = "The Wellcome Trust Sanger Institute";
$hashCenter{"BI"} = "Broad Institute";

open (SUM, ">$summary") or die "open : $!";

print SUM "Run\tRelease_Date\tBases (bp)\tAssembly_Name\tTaxonomyID\tScientific_Name (or species)\tCenter_Name\tConsent\tCountry\tLibrary_Strategy\tLibrary_Selection\tLibrary_Source\tLibrary_Layout\tPlatform\tModel\n";

for my $run (@tabRuns){
  #$str='String 1GIANT FISHString 2'
  #($country)= $str =~ /String 1(.*)String 2/

  my $first_Cmd = "wget -q -O ./$run.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$run' ";

  my $second_Cmd = "wget -q -O ./$run.xml \"https://www.ncbi.nlm.nih.gov//sra?term=($run)%20NOT%20cluster_dbgap%5BPROP%5D&report=FullXml\" ";

  system($first_Cmd);
  system($second_Cmd);

  open (XML, "<$run.xml") or die "open : $!";
  while (<XML>) {
    chomp();
    my $string1 = "Country&gt;</b>";
    my $string2 = "<b>";
    $_ =~ /$string1(.*?)$string2/; #/<div class="xml-tag"><b>\&lt;Country\&gt;<\/b>(.*)<b>\&lt;/;
    if($1){
      $country2 = $1; # my $country $country[0];
    }
    else{
      $country2 = "ND";
    }
    #open COUNTRY, "echo $_ | grep -oP '(?<=Country&gt;</b>).*?(?=<b>)' "; # Non-greedy match (Notice the '?' after '*' in .*)
    #while (<COUNTRY>) {
    #  chomp();
    #  print $_;
    #}
    #my ($substr) = ($string =~ /period_1_(.*)\.ssa/);
  }
  close (XML) or die "close file error : $!";

  open (RUNXML, "<$run.xml") or die "open : $!";
  if ($country2 eq "ND"){
   while (<RUNXML>) {
    chomp();
    my $ostring1 = "</b>geographic location (country and/or sea  region)<b>&lt;/TAG&gt;</b></div><div class=\"xml-tag\"><b>&lt;VALUE&gt;</b>";
    my $ostring2 = "<b>&lt;";
    $_ =~ /$ostring1(.*?)$ostring2/; #/<div class="xml-tag"><b>\&lt;Country\&gt;<\/b>(.*)<b>\&lt;/;
    if($1){
      print "Country2 = ".$1."\n";
      $country2 = $1; # my $country $country[0];
    }
    else{
      $country2 = "ND";
    }
    
  }
  }
  close (RUNXML) or die "close file error : $!";

  #</b>geographic location (country and/or sea  region)<b>&lt;/TAG&gt;</b></div><div class="xml-tag"><b>&lt;VALUE&gt;</b>
  #<b>&lt;
  #close (COUNTRY) or die "close file error : $!";

  open (CSV, "<$run.csv") or die "open : $!";
  while (<CSV>) {
    chomp();
    if ($_ =~  m/$run/) {
      @tabCSV = split (/,/, $_) ;
    }
  }
  close (CSV) or die "close file error : $!";

  my $tmpCenter = "";
  if($hashCenter{$tabCSV[41]}) { $tmpCenter = $hashCenter{$tabCSV[41]}; }

  print SUM "$run\t$tabCSV[1]\t$tabCSV[4]\t$tabCSV[8]\t$tabCSV[27]\t$tabCSV[28]\t$tabCSV[41] ($tmpCenter)\t$tabCSV[44]\t$country2\t$tabCSV[12]\t$tabCSV[13]\t$tabCSV[14]\t$tabCSV[15]\t$tabCSV[18]\t$tabCSV[19]\n";

}


close (SUM) or die "close file error : $!";

my $end = time();

my $total = $end - $start;

print "***** Total time (in seconds) is: $total *****\n";



