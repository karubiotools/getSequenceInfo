#!/usr/bin/bash

echo 'Welcome to getSequenceInfo Unix installer !';
echo 'cpan must first be installed on your computer';
echo '----------------------------------------------------------';
cpan -i Date::Calc;
cpan -i Bio::SeqIO;
cpan -i LWP::Simple;
cpan -i Data::Dumper;
cpan -i IO::Uncompress::Gunzip;
cpan -i IO::File;
cpan -i File::Log;
cpan -i Getopt::Long;
cpan -i Net::FTP;
cpan -i Tk;
echo 'end of install';

