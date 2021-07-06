#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
#
#This is a script to generate quantity control script
#ChangeLog
#20160302 Add Universal&Indexed adaptors 
#20160908 Change structure


my $RawDataPath;
my $Outputdirectory;
my $Shell_name;
my $Split_num = 1;
my $cpu = 4;
# null paramters
my $hc = " "; # host contamination in QC
my $al = " "; # adaptor contamination in QC
my $index = " "; #if remove index not constant
#
my $Quality_cutoff; 
my $low_quality_length_cutoff;
my $ambugious_number_cutoff;
my $trim_bases_cutoff;
my $rawLength = 0;
my $clip = " ";
my $Reference = "null";
my $ADAPTOR = "null";
my $indexedAdaptor = "";
my $adaptorLength = 8;
my $removeNotConstantIndex = "N";
my $cfg;

#script cfg
my $part = 0;
my @binDir;
my $binDir = dirname($0);
my $adaptorContaminationPL = "$binDir/AdaptorContamination.pl";
my $sam2readable = "$binDir/sam2reads.pl";
my $QC = "$binDir/QC_V_2_1.pl";
my $Status = "$binDir/CleanData_status.pl";
push (@binDir,$adaptorContaminationPL);
push (@binDir,$sam2readable);
push (@binDir,$sam2readable);
push (@binDir,$QC);
push (@binDir,$Status);
&Check_bin;
#
my %opt = qw();
GetOptions(\%opt,"i:s","r:s","ot_dir:s","n:s","t:i","cfg:s");

my $use = "
This is a script to generate quantity control script
-i	String(M)  Input file, 7 column file, contains:
	Sample FlowCell Lane InsertSize FQ1 FQ2 IndexSequences
-ot_dir	String(M)  output directory, will gernate a shell/clean sub-directories.
-n	String(M)  Output Main shell name
-t	Int	Split shell into multiple ones
-cfg	file	paramters Configure files.
Example: perl $0 -i RawdataPath -od 00_Clean -n mainClean -t 10 -cfg clean.parameter.cfg &

Warning: Do add a & in the end and make sure it is nohuped.
         Do not add nohup in the begining!
";

if (scalar keys %opt == 0 ){
print DARK,GREEN,"$use\n",RESET;
exit;
}
my $error;
my @error;
my %fq1;
my %fq2;
my %index;
my %full_index;
&checkOptions;
&CFG;
&ERROR;
&GENER;
&SPLIT;

##################
sub GENER{
	open OT,">$Outputdirectory/$Shell_name.main.sh" or die;
	foreach my $m (keys %fq1){
	`mkdir $Outputdirectory/Clean/$m`;
	my $sampleDir = "$Outputdirectory/Clean/$m";
	my $fq1 = $fq1{$m};
	my $fq2 = $fq2{$m};
	my $index = $index{$m};
	open SH,">$Outputdirectory/Shell/$m.clean.sh";
		if ($Reference ne "null"){
		print SH "#host contamination removal\n";
		print SH "/gpfsdata/lia/BIN/bwa-0.6.1/bwa aln -t 4 $Reference $fq1 > $sampleDir/$m.aln.sai1 2> $sampleDir/$m.aln.sai1.log\n";
		print SH "/gpfsdata/lia/BIN/bwa-0.6.1/bwa aln -t 4 $Reference $fq2 > $sampleDir/$m.aln.sai2 2> $sampleDir/$m.aln.sai2.log\n";
		print SH "/gpfsdata/lia/BIN/bwa-0.6.1/bwa sampe $Reference $sampleDir/$m.aln.sai1 $sampleDir/$m.aln.sai2 $fq1 $fq2 > $sampleDir/$m.sam 2> $sampleDir/$m.sam.log\n";
		print SH "perl $sam2readable $sampleDir/$m.sam $sampleDir/$m.host.contamination\n";
		print SH "rm -f $sampleDir/$m.aln.sai1 $sampleDir/$m.aln.sai2 $sampleDir/$m.sam $sampleDir/$m.aln.sai1.log $sampleDir/$m.aln.sai2.log\n";
		$hc = "-hc $sampleDir/$m.host.contamination";
		}
		
		if ($ADAPTOR eq "T"){
		my $AdaptorForFq1; my $AdaptorForFq2;
		print SH "#index adaptor contamination removal\n";
			#configure index origin
			if( (length $index == 6) or (length $index == 8) ){ # single index
			die "$index do not have a sequence" unless exists $full_index{$index};
			$AdaptorForFq1 = $full_index{$index};
			die "$index do not have a sequence" unless exists $full_index{"-"};
			$AdaptorForFq2 = $full_index{"-"};
			print SH "#single index\t$m\t$index,$AdaptorForFq1\tUniversal,$AdaptorForFq2\n";
			}
			elsif( (length $index == 16) or (length $index == 12) ){ # duaal index
			my $indexForFq1 = substr($index,0,8);
			my $indexForFq2 = substr($index,8,8);
			die "$indexForFq1 do not have a sequence" unless exists $full_index{$indexForFq1};
			$AdaptorForFq1 = $full_index{$indexForFq1};
			die "$indexForFq2 do not have a sequence" unless exists $full_index{$indexForFq2};
			$AdaptorForFq2 = $full_index{$indexForFq2};
			print SH "#dual index\t$m\t$indexForFq1,$AdaptorForFq1\t$indexForFq2,$AdaptorForFq2\n";
			}
			else{
			die "$index do not have a 6/8, 12/16 length";
			}
		$AdaptorForFq2 = reverse $AdaptorForFq2;
		$AdaptorForFq2 =~ tr/ATGC/TACG/;
		print SH "cutadapt -b $AdaptorForFq1 --discard -O $adaptorLength -o $sampleDir/$m.temporary.1.fq  $fq1 > $sampleDir/$m.temporary.1.cutadaptor.log\n";
		print SH "cutadapt -b $AdaptorForFq2 --discard -O $adaptorLength -o $sampleDir/$m.temporary.2.fq  $fq2 > $sampleDir/$m.temporary.2.cutadaptor.log\n";
		print SH "perl $adaptorContaminationPL $fq1 $sampleDir/$m.temporary.1.fq $fq2 $sampleDir/$m.temporary.2.fq > $sampleDir/$m.adaptor.list\n";
		
		print SH "rm $sampleDir/$m.temporary.1.fq $sampleDir/$m.temporary.2.fq\n";
		$al = "-al $sampleDir/$m.adaptor.list";
		}

		if ($removeNotConstantIndex eq "Y"){
		$index = "-index $index";
		}
		
		print SH "perl $QC -fq1 $fq1 -fq2 $fq2 -Q $Quality_cutoff -lqc $low_quality_length_cutoff -TBc $trim_bases_cutoff -amc $ambugious_number_cutoff $clip $hc $al $index -fmtq -g -out $Outputdirectory/$m >$Outputdirectory/$m.clean.log\n";
		print SH "echo \"$m whole clean process is finished\"\n";
	close (SH);
	print OT "nohup sh $Outputdirectory/Shell/$m.clean.sh > $Outputdirectory/Shell/$m.clean.log\n";
	}
}
#################
sub SPLIT{
my $sample_number = scalar (keys %fq1);
my $number = int ( $sample_number / $Split_num );
   $number = 1 if $number == 0;
`nohup split --verbose -d -l $number $Outputdirectory/$Shell_name.main.sh $Outputdirectory/$Shell_name.FLAG >  $Outputdirectory/split.log`;
open IN,"$Outputdirectory/split.log" or die;
while (<IN>){
chomp;
my @information = split /\'|\`/;
`mv $information[-1] $information[-1].sh`;
$part ++;
}
close (IN);
`rm $Outputdirectory/split.log $Outputdirectory/$Shell_name.main.sh`;
my $flagShells = `ls $Outputdirectory/$Shell_name.FLAG*.sh`;
my @flagShells = split /\n/,$flagShells;
	foreach my $flag (@flagShells){
	`nohup sh $flag > $flag.log &`;
	}
while (1){
sleep(60);
my $finishedNumber = `grep 'whole clean process is finished' $Outputdirectory/Shell/*.clean.log | wc -l`;
last if $finishedNumber == $sample_number;
}
`perl $Status $RawDataPath $Outputdirectory $Outputdirectory/clean.status`;
}
#################
sub Check_bin{
	foreach my $m (@binDir){
	die "$m not E" unless -e $m;
	}
	
}
#################
sub checkOptions{
	if (exists $opt{"i"}){
	$RawDataPath = $opt{"i"};
		open IN,"$RawDataPath" or die $!; <IN>; # skip head
		while (<IN>){
		chomp;
		my @array = split /\t/,$_;
			die "rawdatapath must be a 7 column file!" unless scalar @array == 7;
			unless (  (-e $array[4]) && (-e $array[5]) ){
			$error = "Input file $array[4] or $array[5] does not exist";
			push (@error,$error);
			}
		my $TAG = "$array[0]-$array[1]-$array[2]-$array[6]"; #sample-flowcell-lane-index
		$fq1{$TAG} = $array[4];
		$fq2{$TAG} = $array[5];
		$index{$TAG} = $array[6];
		}
		close (IN);
	}
	else{
	$error = "-i is a mandotary option";
	push (@error,$error);
	}
	#
	if (exists $opt{"ot_dir"}){
	$Outputdirectory = $opt{"ot_dir"};
	`mkdir $Outputdirectory/Shell $Outputdirectory/Clean`;
		unless (-w $Outputdirectory){
		$error = "You do not have write permission in ourput directory $Outputdirectory";
		push (@error,$error);
		}
	}
	else{
	$error = "-ot_dir is a mandotary option";
	push (@error,$error);
	}
	#
	if (exists $opt{"n"}){
	$Shell_name = $opt{"n"};
	}
	else{
	$error = "-n is a mandotary option";
	push (@error,$error);
	}
	#
	if (exists $opt{"t"}){
	$Split_num = $opt{"t"};
	}
	#
	unless (exists $opt{"cfg"}){
	$error = "cfg must be initialized";
	push (@error,$error);
	}
}
########################
sub CFG{
	open IN,"$opt{cfg}" or die;
	while (<IN>){
	chomp;
	my @cfg = split /\t/,$_;
		if ($cfg[0] eq "Reference"){
		$Reference = $cfg[1];
			unless (-e $Reference){
			$error = "reference $Reference do not exists.";
			push (@error,$error);
			}
		}
		elsif ($cfg[0] eq "Quality_cutoff"){
		$Quality_cutoff = $cfg[1];
		}
		elsif ($cfg[0] eq "low_quality_length_cutoff"){
		$low_quality_length_cutoff = $cfg[1];
		}
		elsif ($cfg[0] eq "ambugious_number_cutoff"){
		$ambugious_number_cutoff = $cfg[1];
		}
		elsif ($cfg[0] eq "trim_bases_cutoff"){
		$trim_bases_cutoff = $cfg[1];
		}
		elsif ($cfg[0] eq "clip"){
		$clip = "-clip $cfg[1]";
		}
		elsif ($cfg[0] eq "rawLength"){
		$rawLength = $cfg[1];
		}
		elsif ($cfg[0] eq "adaptorListFile"){#adaptorListFile #adaptorListFile
		$ADAPTOR = "T";
		&READ_INDEX($cfg[1]);
		}
		elsif ($cfg[0] eq "adaptorLengthCutoff"){
		$adaptorLength = $cfg[1];
		}
		elsif ($cfg[0] eq "RemovenotConstantindex"){
		$removeNotConstantIndex = $cfg[1];
		}
	}
	close (IN);
}
########################
sub ERROR{
if (scalar @error > 0){
my $error = join ("\n",@error);
print DARK, RED,"Fatal error\n$error\n",RESET;
print DARK,GREEN,"$use\n",RESET;
exit;
}
}
########################
sub READ_INDEX{
open INDEX,"$_[0]" or die;
while (<INDEX>){
chomp;
unless ($_ =~ "^#"){
	my @array = split /\t/,$_;
	$index{$array[1]} = $array[0];
        	if ($array[1] ne "-"){
	        my @Full = split /\[|\]/,$array[2];
	        $full_index{$array[1]} = "$Full[0]$array[1]$Full[2]";
	        }
	        else{
	        $full_index{$array[1]} = $array[2];
	        }
}
}
close (INDEX);
}
