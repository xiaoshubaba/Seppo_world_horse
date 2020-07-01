#!/usr/bin/perl -w
use strict;
my %mappingRatio;

open IN,"HorseV3.mapping.log" or die; <IN>;
while (<IN>){
chomp;
my @in = split /\t/,$_;
$mappingRatio{$in[0]} = $in[-1];
}
close (IN);

my %generichness;
my %CAZYenezy;
#sample  cohort  status  OBSenzyme       OBSgene
open I,"CAZY.status" or die;
while (<I>){
chomp;
my @a = split /\t/,$_;
$generichness{$a[0]} = $a[4];
$CAZYenezy{$a[0]} = $a[3];
}
close (I);

my %obs;
#sample  bacShannon      bacSimpson      bacObs  proShannon      proSimpson      proObs  wholeShannon    wholeSimpson    wholeObs geneRichness5   geneRichness    cohort  status
open A,"Mer.alpha"; <A>;
while (<A>){
chomp;
my @b = split /\t/,$_;
$obs{$b[0]} = $b[9];
}
close (A);

my %res;
open RES,"HorseV3.resgene.number" or die; <RES>;
while (<RES>){
chomp;
my @res = split /\t/,$_;
$res{$res[0]} = $res[-1];
}

&OT;
sub OT{
open I,"sample.metadata" or die;
my $h = <I>; $h =~ s/\n//;
print "$h\tObservedSpecies\tGeneRichness\tMappingRatio\tAntibiotisResistantGeneNumber\tCAZYgeneNumber\n";
while (<I>){
chomp;
my @I = split /\t/,$_;
print "$_\t$obs{$I[0]}\t$generichness{$I[0]}\t$mappingRatio{$I[0]}\t$res{$I[0]}\t$CAZYenezy{$I[0]}\n";
}
close (I);
}
