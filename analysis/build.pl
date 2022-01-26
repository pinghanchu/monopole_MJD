#!/usr/bin/perl
print "Please input [mass] [energy]\n";
$numArgs = $#ARGV + 1;
if( $numArgs < 2){
    print "You miss arguments\n";
}elsif( $numArgs >2){
    print "You have too many arguments\n";
}

my $mass = $ARGV[0];
for (my $i=-5;$i<=10;$i++){
    my $energy = "1e".$i."GeV";
    my $outputcsv = "../data/data_".$mass."_".$energy.".csv";
    my $outputroot = "../data/monopole_".$mass."_".$energy.".root";
    print($outputroot,"\n");
    print($outputcsv,"\n");
    system("cp $outputroot monopole.root");
    system("root -b -q fill.cc");
    system("mv data.csv $outputcsv");
}
