#!/usr/bin/perl
#get data from EVENTBIN generated by FLUKA
#EVENTBIN saves each event including hit energy and pixel.
print "Please input [energy]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 1){
    print "You miss arguments\n";
}elsif( $numArgs >1){
    print "You have too many arguments\n";
}
my $ERG = $ARGV[0];
print($ERG,"\n");
use warnings;
use strict;
cddir "../monopole-build/";

my $filename = "./monopole.in";
my $outputfile = "./monopole_".$ERG."GeV.in";
print $filename,"\n";
open(my $fin,'<',$filename) or die $!;
open(my $fout,">",$outputfile);
my $str1 = "/gun/energy 1 GeV";
my $str2 = "/gun/energy ".$ERG." GeV";

while(my $line = <$fin>){
    
    if ($line =~ m/$str1/) {
	#print "match\n";
	my $old = $str1;
	my $new = $str2;
	$line =~ s/$old/$new/g;
	#print $line,"n";
	
    } else {
	#print '';
    }
    #print $line,"\n"
    print $fout $line;
}
close $fin;


