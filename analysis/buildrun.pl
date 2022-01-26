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
#my $mass = $ARGV[0];
#print($ERG,"\n");
use warnings;
use strict;

for (my $i=10; $i<=30;$i++){
    print($i);
    my $mass = "1e".$i."GeV";
    system("./build.pl $mass");
}
