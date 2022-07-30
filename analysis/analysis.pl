#!/usr/bin/perl
#get data from EVENTBIN generated by FLUKA
#EVENTBIN saves each event including hit energy and pixel.
print "Please input [mass]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 1){
    print "You miss arguments\n";
}elsif( $numArgs >1){
    print "You have too many arguments\n";
}
use warnings;
use strict;

for (my $j=10;$j<=15;$j++){
    my $mass = "1e".$j;
    for (my $i=-5; $i<=5;$i++){
	my $erg = "1e".$i;
	my $rootfile = "../data/data_scatter/monopole_".$mass."GeV_".$erg."GeV.root";
	my $outputcsv = "./data/analysis_".$mass."GeV_".$erg."GeV.csv";
	print($rootfile,"\n");
	system("rm monopole.root");
	system("ln -s $rootfile monopole.root");
	system("root -b -q analysis.cc");
	system("mv analysis.csv $outputcsv");
	#system("rm monopole.root");
    }
}
