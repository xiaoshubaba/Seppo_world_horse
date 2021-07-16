######################################################
#program: Illumina high-throughpu paired end reads Quanlity Control Version 2.3
#Author: Li Ang (lia@zju.edu.cn)
#Change Log
#20120210	Version 1.0
#Initial script was finished
#20120216	Version 1.1
#Fix a avergae position calculation Error
#20120410	Version 2.0
#Make vital change of the script structure
#Add a parameter to contol the CASAVA quality value cutoff , 64/32
#Add some the coments
#20130106	Version 2.1
#Add functions to trim N in 5 end of each read
#20150103	Version 2.2
#Add clip, just use the first XX bases of all reads and quality
#20160906	Version 2.3
#Add if remove unmatched index
#20170704	#skip vesrion
#Check if input fq file, is a link, then return error
#20181012	Verion 2.4
#Add BGI sequencer way to calculate quality
#20190220	#skip vesrion
#Change suffix bug: fastaq to fastq
######################################################
#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Term::ANSIColor qw(:constants);

# files and options
my $fq1; # input fastaq file 1
my $fq2; # input fastaq file 2
my $fq;  # input fastaq file 
my $sequencing_index; #sequence index
my $clip; #clip the first XX reads
my $Q;   # the asiic code of the low quantity value cutoff , the real quantity is  ORD{"#"} - 33;
my $lqc; # Cutoff of  num of LOW Quantity in the whole read
my $TBc; #The num of LOW Quantity or ambugilous in the tail of read,if it is set as f, means all LQ in the tail will be trimmed, if it is set as a, means the average position will be used ,  when the tail of reads is "N" or the quantity value is lower than the setted "Q", it will be triggered
my $cc;  #Cutoff of continuous A T G C
my $amc; #Cutoff of N
my $al;  # input file of reads belong to adapter contamination
my $hc;  # input file of reads belong to host contamination
my $fmtq;# Triger the fastaq format output
my $g;   # Triger the compressed format output
my $out; # Output file prefix

my $trim_head_a; # Trim N in head of sequences if it is N before check the number of N in sequences
my $trim_head_b; # Trim N in head of sequences if it is N afer check the number of N in sequences

my %Pre_input_al; # Reads belong to adapter contamination will be stored into this hash
my %Pre_input_hc; # Reads belong to HOST contamination will be stored into this hash
# Statistics;

my $S_index1 = 0; my $S_index2 = 0; my $S_index = 0; #Number of sequencing indices
my $S_ac1 = 0;my $S_ac2 = 0;my $S_ac = 0; # Number of the adapter contamination reads
my $S_lqc1 = 0;my $S_lqc2 = 0;my $S_lqc = 0; # Number of the low quantity reads
my $S_TBc1 = 0;my $S_TBc2 = 0;my $S_TBc = 0; # Value to trim the tail,f you decide to trim the tail, means TBc is set, these three variables will work 
my $S_cc1 = 0;my $S_cc2 = 0;my $S_cc = 0; # Number of the continuous nucleitides
my $S_amc1 = 0;my $S_amc2 = 0;my $S_amc = 0; # Number of the ambugilous
my $S_hcc = 0;# Number of the host contamination , in pe sequencing , reads from another direction will also be seen as contamination
my $S_PE_number_FQ1 = 0; my $S_PE_number_FQ2 = 0; my $S_SE_number_FQ1 = 0; my $S_SE_number_FQ2 = 0; # Store the pe and se information
my $S_PE_number_FQ1_size = 0; my $S_PE_number_FQ2_size = 0; my $S_SE_number_FQ1_size = 0; my $S_SE_number_FQ2_size = 0; #Store the data size information

my $Reads_num = 0; # the whole reads 
my $length = 0; # a variable to store the sequencing length
my $chr; # the chr of the quantity value
my ($outputname,$directory,$suffix); # Parse the outputfile
my ($sum_position,$trim_initial,$last_position_N,$last_position_low,$index,$last_nuclietide,$last_quantity,$fatter,$initial); # variables to get the strcter trim position
my ($SizeTrimA,$SizeTrimB,$SizeTrim); # nuclietide number of trim
my $outputfile;
my $head_information; # information to store the sequencing information, file path , etc.
# Generate the options
my %opt = qw();
GetOptions(\%opt,"fq1:s","fq2:s","fq:s","index:s","clip:i","Q:i","lqc:i","TBc:s","cc:i","amc:i","al:s","hc:s","fmtq!","g!","out:s","tNa!","tNb!","help|h!");

my $usage = "
This is a perl script to fileter the host contamination in reads quantity control

   --fq1	(String)	Input Pair end File 1
   --fq2	(String)	Input Pair end File 2
   --fq		(String)	Input Single file
   --index	(String)	Input sequencint index (oligo-nuclietide, 6 or 8 bps)
   --Q		(Int)		Low quantity value,default:2
   --clip	(Int)		Just use the first XX bases of all reads and quanlity
   --lqc	(Int)		Cutoff of  num of LOW Quantity in the whole read
   --TBc	(String)	
   Note of TBc:	                The num of LOW Quantity in the tail of read,if it is set as \"f\", means all LQ in the tail will be trimmed, if it is setas \"a\", means the average position will be used to trim ,when the tail of reads is \"N\" or the quantity value is lower than the seted \"Q\", it will be triggered
   --cc		(Int)		Cutoff of continuous A T G C
   --amc	(Int)		Cutoff of N
   --al		(String)	Reads belong to adapter contamination
   --hc		(String)	Reads belong to host contamination
   --tNa			Trim N in 5 end of each reads before check number of ambigulous bases
   --tNb			Trim N in 5 end of each reads after check number of ambigulous bases
   --fmtq			Output file as a fastaq format
   --g				Gzip the output file
   --out	(String)	Output file, prefix
   --help or h			Show the help information
";

if ( (scalar (keys %opt) == 0) or (exists $opt{"h"}) or (exists $opt{"help"}) ){
print DARK,GREEN,"$usage\n",RESET;
exit;
}
############
my $error;
my @error;
############
&Check_files;
&Check_options;
&Output_ERROR;
&Read_al(\%Pre_input_al,$al) if defined $al;
&Read_hc(\%Pre_input_hc,$hc) if defined $hc;
&Standard_Output;

if ( defined $fq1){ # fq1 and fq2 both pass the check
	if ( (defined $TBc) && ($TBc eq "a") ){
	&Calculate($fq1,$S_TBc1);
	&Calculate($fq2,$S_TBc2);
	}
&Print_TBc;
&Write_PE($fq1,$fq2);
}
elsif (defined $fq){
	if ( (defined $TBc) && ($TBc eq "a") ){
	&Calculate($fq,$S_TBc);
	}
&Print_TBc;
&Write_SE($fq);
}
###################
print "\nOutput file list\n$outputfile\n\n########################################\nQC information\n\n"; 
my $clean_1;
my $clean_2;
my $clean;
if (defined $fq1){
print "\tReads_1\tReads_2\nRawReadsNum\t$Reads_num\t$Reads_num\n";
print "index\t$S_index1\t$S_index2\n" if defined $sequencing_index;
print "Adapter\t$S_ac1\t$S_ac2\n" if defined $al;
print "LowQuantity\t$S_lqc1\t$S_lqc2\n" if defined $lqc;
print "Continuous\t$S_cc1\t$S_cc2\n" if defined $cc;
print "Ambiguous\t$S_amc1\t$S_amc2\n" if defined $amc;
$S_hcc = $S_hcc / 2;
print "Host\t$S_hcc\t$S_hcc\n" if defined $hc;
$clean_1 = ($Reads_num - $S_index1 - $S_ac1 - $S_lqc1 - $S_cc1 - $S_amc1 - $S_hcc);
$clean_2 = ($Reads_num - $S_index2 - $S_ac2 - $S_lqc2 - $S_cc2 - $S_amc2 - $S_hcc);
print "CleanData(num)\t$clean_1\t$clean_2\n";
print "TrimmedBP(bp)\t$SizeTrimA\t$SizeTrimB\n" if defined $TBc;
my $left_1 = ($Reads_num - $S_ac1 - $S_lqc1 - $S_cc1 - $S_amc1 -$S_hcc) * $length - $SizeTrimA;
my $left_2 = ($Reads_num - $S_ac2 - $S_lqc2 - $S_cc2 - $S_amc2 -$S_hcc) * $length - $SizeTrimB;
print "CleanData(bp)\t$left_1\t$left_2\n";
print "########################################\n";
print "CleanData(numPE)\t$S_PE_number_FQ1\t$S_PE_number_FQ2\n";
print "CleanData(bpPE)\t$S_PE_number_FQ1_size\t$S_PE_number_FQ2_size\n";
print "CleanData(numSE)\t$S_SE_number_FQ1\t$S_SE_number_FQ2\n";
print "CleanData(bpSE)\t$S_SE_number_FQ1_size\t$S_SE_number_FQ2_size\n";
}
elsif (defined $fq){
print "\tReads\nNum\t$Reads_num\n";
print "index\t$S_index\n" if defined $S_index;
print "Adapter\t$S_ac\n" if defined $al;
print "LowQuantity\t$S_lqc\n" if defined $lqc;
print "Continuous\t$S_cc\n" if defined $cc;
print "Ambiguous\t$S_amc\n" if defined $amc;
print "Host\t$S_hcc\n" if defined $hc;
$clean = ($Reads_num - $S_ac - $S_lqc - $S_cc - $S_amc -$S_hcc);
print "CleanData(num)\t$clean\n";
print "TrimmedBP(bp)\t$SizeTrim\n" if defined $TBc;
my $left = ($Reads_num - $S_ac - $S_lqc - $S_cc - $S_amc -$S_hcc) * $length - $SizeTrim;
print "CleanData(bp)\t$left\n";
}
print "########################################\n";
###################	End


sub Print_TBc{
print "Tail low quantity trim length:\t";
	if ($TBc eq "a"){
	print "$S_TBc1&&$S_TBc2\n" if defined $fq1;
	print "$S_TBc\n" if defined $fq;
	}
	else{
	print "$TBc\n";
	}	
}

sub Read_al{
my ($a,$b) = @_;
	open AL,"$b" or die;
	while (<AL>){
	chomp;
	$$a{$_} = 1;
	}
	close (AL);
}

sub Read_hc{
my ($a,$b) = @_;
	open HC,"$b" or die;
	while (<HC>){
	chomp;
	$$a{$_} = 1;
	}
	close (HC);
}

sub Check_files{
	if(  ( (exists $opt{"fq1"}) and (exists $opt{"fq2"}) ) or (exists $opt{"fq"}) ){
		if ( (exists $opt{"fq1"}) and (exists $opt{"fq2"}) ){
			if ( (-e $opt{"fq1"}) and (-e $opt{"fq2"}) ){
			$fq1 = $opt{"fq1"};
			$fq2 = $opt{"fq2"};
			&Get_sequencing_length($fq1);
			$head_information = "Sequencing stratage:\tPair end Sequencing\nSequencing length:\t$length bp\n";
			$head_information .= "Input file Fastaq 1:\t$fq1\nInput file Fastaq 2:\t$fq2\n";
			}
			else{
			$error = "Input Pair-end file\n" .$opt{"fq1"}. "\nor\n" .$opt{"fq2"}. "\ndo not exist";
			push (@error,$error);
			}
		}
		elsif (exists $opt{"fq"}){
			if (-e $opt{"fq"}){
			$fq = $opt{"fq"};
			&Get_sequencing_length($fq);
			$head_information = "Sequencing stratage:\tSingle end Sequencing\t$length bp\n";
			$head_information .= "Input file:\tFastaq:\t$fq";
			}
			else{
			$error = "Input Single-end file \n".$opt{"fq"}."\n does not exist";
			push (@error,$error);
			}
		}
	}
	else{
	$error = "Input file error: you must input two pair-end files or one single-end file";
	push (@error,$error);
	}

	if(exists $opt{"al"}){
		if (-e $opt{"al"}){
		$al = $opt{"al"};
		$head_information .= "Input file adapter contamination reads\t$al\n";
		}
		else{
		$error = "Adapter contamination file ".$opt{"al"}." does not exist";
		push (@error,$error);
		}
	}

	if (exists $opt{"hc"}){
		if (-e $opt{"hc"}){
		$hc = $opt{"hc"};
		$head_information .= "Input file host contamination reads\t$hc\n";
		}
		else{
		$error = "Host contamination file ".$opt{"hc"}." does not exist";
		push (@error,$error);
		}
	}
}

sub Get_sequencing_length{
$length = `gzip -dc $_[0] |head -2| wc -L`; # Get the sequencing length
$length =~ s/\n//;
my @length = split /\s+/,$length;
$length = $length[0];
}


sub Check_options{
		#
		if (exists $opt{"clip"}){
			if ($opt{"clip"} >= 50){
			$clip = $opt{"clip"};
			$head_information .= "Just use the first 5' bases of sequences and quality\t$Q\n";
			}
			else{
			$error = "clip is no less than 50";
			push (@error,$error);	
			}
		}
		#
		if (exists $opt{"index"}){
		$sequencing_index = $opt{"index"}; 
		#sequencing_index
		$sequencing_index = uc($sequencing_index);
		$head_information .= "sequencing index:\t$sequencing_index\n";
		my @index = split /\t/,$index;
		foreach my $l_index (@index){
			unless ($index[$l_index] =~ m/[ATGC]/){
			$error = "index $index[$l_index] position $l_index is not [ATGC]";
			} 
		}
		}
		#
		if (exists $opt{"Q"}){
			if ($opt{"Q"} > 0){
			$Q = $opt{"Q"};
			my $new_Q = 33 + $Q;
			$chr = chr($new_Q); # Get the quantitu character
			$head_information .= "\nQuantity character:\t$chr\n";
			$head_information .= "Quantity value cutoff :\t$Q\n";
			}
			else{
			$error = "Q must be a positive integer";
			push (@error,$error);
			}
		}
		else{
		$Q = 2;
		$chr = "#";
		}
		if (exists $opt{"lqc"}){
			if ($opt{"lqc"} > 0){
			$lqc = $opt{"lqc"};
			$head_information .= "Low quantity cutoff:\t$lqc\n";
			}
			else{
			$error = "lqc must be a positive interger";
			push (@error,$error);
			}
		}
		if (exists $opt{"TBc"}){
			if ( ($opt{"TBc"} eq "f") || ($opt{"TBc"} eq "a") || ($opt{"TBc"} > 0) ){
			$TBc = $opt{"TBc"};
			}
			else{
			$error = "TBc must be a positive integer || f || a";
			push (@error,$error);
			} 
		# We should add the information when avergae position is calculated 
		}
		if (exists $opt{"cc"}){
			if ($opt{"cc"} > 0){
			$cc = $opt{"cc"};
			$head_information .= "Continuous nuclietides cutoff:\t$cc\n";
			}
			else{
			$error = "cc must be a positive integer";
			push (@error,$error);
			}
		}
		if (exists $opt{"amc"}){
			if ($opt{"amc"} > 0){
			$amc = $opt{"amc"};
			$head_information .= "Ambugilous nuclietides (N) cutoff:\t$amc\n";
			}
			else{
			$error = "amc must be a positive integer";
			push (@error,$error);
			}
		}
		if (exists $opt{"out"}){
		Parse_file($opt{"out"});
			if (-w $directory){
			$out = $opt{"out"};
			}
			else{
			$error = "you do not have a write permission in directory $directory";
			push (@error,$error);
			}
		}
		else{
		$error = "-out is missing, where should I output my clean data?";
		push (@error,$error);
		}
		if ($opt{"fmtq"}){
		$fmtq = "fastq";
		}
		else{
		$fmtq = "fasta";
		}
		$head_information .= "Output file format:\t$fmtq\n";
		if ($opt{"g"}){
		$g = 1;
		$head_information .= "Output file Compressed:\tYes\n";
		}
		else{
		$g = 0;
		}
		if ($opt{"tNa"}){
		$trim_head_a = 1;
		}
		if ($opt{"tNb"}){
		$trim_head_b = 1;
		}

		if ($opt{"tNa"} && $opt{"tNb"}) {
		$error = "you can not both choose tNa and tNb "; 
		push (@error,$error);
		}
}

sub Parse_file{
($outputname,$directory,$suffix) = fileparse($_[0]);
$head_information .= "Output directory:\t$directory\n";
$head_information .= "Output file prefix:\t$outputname\n";
}

sub Output_ERROR{
	if (scalar @error > 0){
	$error = join ("\n",@error);
	print DARK, RED,"Fatal error\n$error\n",RESET;
	print DARK,GREEN,"$usage\n",RESET;
	exit;
	}
}

sub Standard_Output{
print "\n########################################\nBasic information\n\n$head_information";
}


sub Calculate{
	my $file = `file $_[0]`;
	if ($file =~ "gzip compressed data"){
	open IN,"gzip -dc $_[0]|" or die $!;
	}
	else{
	open IN,"$_[0]" or die $!;
	}
	$sum_position = 0; # sum of the trim position
	$initial = 0; # a calculator to store the number of reads have a trim position
	while (<IN>){ # begin to read the fq file
	my $query = $_;
	chomp;
	my $sequence = <IN>; # read sequence 
	   $sequence =~ s/\n//; # read quantity
	<IN>;
	my $quantity = <IN>;
	   $quantity =~ s/\n//;
		if ( (exists $Pre_input_al{$query}) || (exists $Pre_input_hc{$query}) ){
		}# Do not need to calculte thie read cause it belongs to contamination or adapter 
		else{
		my $last_Q = substr($quantity,-1,1);
		my $last_QV = ord($last_Q) - 33;
			if  ( ($last_QV <= $Q) or ($sequence =~  "N\$") ){ 
			&Who_is_strict($quantity,$sequence,$initial);
			$sum_position += $fatter;
			}
		}
	}
	close (IN); # read all reads
	$_[1] = $sum_position / $initial; # calculate the average position
}

sub Who_is_strict{
my @quantity = split (//,$_[0]);
my @sequence = split (//,$_[1]);
$trim_initial = 0; # begin
$last_position_N = 0;
$last_position_low = 0;
	while ($trim_initial < $length){
	$index = -1 - $trim_initial;
	$last_nuclietide =  $sequence[$index];
	$last_quantity   =  ord($quantity[$index]) - 33;
		if ( ($last_nuclietide ne "N") and ($last_quantity > $Q) ){
		last;
		}
		else{
		$_[2] ++ if (defined $_[2]);
			if ($last_nuclietide eq "N"){
			$last_position_N = $trim_initial;
			}
			if ($last_quantity <= $Q){
			$last_position_low = $trim_initial;
			}
		}
	$trim_initial ++;
	}
	$last_position_N ++;
	$last_position_low ++;
	if ($last_position_N >= $last_position_low){
	$fatter = $last_position_N;
	}
	else{
	$fatter = $last_position_low;
	}
}

sub Write_PE{
my $file;

	$file = `file $_[0]`;
	if ($file =~ "gzip compressed data"){
	open IN_1,"gzip -dc $_[0]|" or die $!;
	}
	else{
	open IN_1,"$_[0]" or die $!;
	}
	$file = `file $_[1]`;
	if ($file =~ "gzip compressed data"){
	open IN_2,"gzip -dc $_[1]|" or die $!;
	}
	else{
	open IN_2,"$_[1]" or die $!;
	}
	if ($g == 1){
	open OT_A,"| gzip >$out.1.$fmtq.gz" or die $!;
	open OT_B,"| gzip >$out.2.$fmtq.gz" or die $!;
	open OT_S,"| gzip >$out.single.$fmtq.gz" or die $!;
	$outputfile = "$out.1.$fmtq.gz\n$out.2.$fmtq.gz\n$out.single.$fmtq.gz";
	}
	else{
	open OT_A,">$out.1.$fmtq" or die $!;
	open OT_B,">$out.2.$fmtq" or die $!;
	open OT_S,">$out.single.$fmtq" or die $!;
	$outputfile = "$out.1.$fmtq\n$out.2.$fmtq\n$out.single.$fmtq";
	} 
	while (<IN_1>){
	chomp;
	$Reads_num ++;
	my $skip_A = "F";
	my $skip_B = "F";
	my $query_A = $_;
	my @query_A = split (/\@/,$query_A);
	my $sequence_A = <IN_1>;
	   $sequence_A =~ s/\n//;
	my $trimmed_sequence_A = $sequence_A;
	my $direction_A = <IN_1>;
	   $direction_A =~ s/\n//;
	my $quantity_A = <IN_1>;
	   $quantity_A =~ s/\n//;
	my $trimmed_quantity_A = $quantity_A;

	my $query_B = <IN_2>;
	   $query_B =~ s/\n//;
	my @query_B = split (/\@/,$query_B);
	my $sequence_B = <IN_2>;
	   $sequence_B =~ s/\n//;
	my $trimmed_sequence_B = $sequence_B;
	my $direction_B = <IN_2>;
	   $direction_B =~ s/\n//;
	my $quantity_B = <IN_2>;
	   $quantity_B =~ s/\n//;
	my $trimmed_quantity_B = $quantity_B;
	
	if ($trim_head_a == 1){
	&TRIM_HEAD($sequence_A);
	&TRIM_HEAD($sequence_B);
	}

	if ( (defined $al) or (defined $hc) ){
	&CHECK_Q($query_A[1],$skip_A,$S_ac1,$S_hcc);
	&CHECK_Q($query_B[1],$skip_B,$S_ac2,$S_hcc);
	}
	
	if (defined $sequencing_index){
	&CheckIndex($query_A[1],$query_B[1],$skip_A,$skip_B,$S_index1,$S_index2);
	}

	if (defined $lqc){
	&CHECK_LQC($quantity_A,$skip_A,$S_lqc1) if ($skip_A eq "F");
	&CHECK_LQC($quantity_B,$skip_B,$S_lqc2) if ($skip_B eq "F");
	}
	if (defined $cc){
	&CHECK_CC($sequence_A,$skip_A,$S_cc1) if ($skip_A eq "F");
	&CHECK_CC($sequence_B,$skip_B,$S_cc2) if ($skip_B eq "F");
	}
	if (defined $amc){
	&CHECK_AMC($sequence_A,$skip_A,$S_amc1) if ($skip_A eq "F");
	&CHECK_AMC($sequence_B,$skip_B,$S_amc2) if ($skip_B eq "F");
	}
	
        if ($trim_head_b == 1){
        &TRIM_HEAD($sequence_A);
        &TRIM_HEAD($sequence_B);
        }
	


	if (defined $TBc){ #  trim only if N is at the end or lower is found 
	my $last_nucli_A = substr($trimmed_quantity_A,-1,1);
	my $last_nucli_B = substr($trimmed_quantity_B,-1,1);
	my $last_nucli_Q_A = ord($last_nucli_A) - 33;
	my $last_nucli_Q_B = ord($last_nucli_B) - 33;
	&TRIM($trimmed_sequence_A,$trimmed_quantity_A,$S_TBc1,$SizeTrimA) if ( ($skip_A eq "F") && ( ($last_nucli_Q_A <= $Q) || ($trimmed_sequence_A =~ "N\$") ) );
	&TRIM($trimmed_sequence_B,$trimmed_quantity_B,$S_TBc2,$SizeTrimB) if ( ($skip_B eq "F") && ( ($last_nucli_Q_B <= $Q) || ($trimmed_sequence_B =~ "N\$") ) );
	}

	if ( ($skip_A eq "F") and ($skip_B eq "F") ){

	if (exists $opt{clip}){
	if (length $trimmed_sequence_A > $clip){
	$trimmed_sequence_A = substr($trimmed_sequence_A,0,$clip);
	$trimmed_quantity_A = substr($trimmed_quantity_A,0,$clip);
	}
	if (length $trimmed_sequence_B > $clip){
	$trimmed_sequence_B = substr($trimmed_sequence_B,0,$clip);
	$trimmed_quantity_B = substr($trimmed_quantity_B,0,$clip);
	}

	}

	print OT_A "\@$query_A[1]\n$trimmed_sequence_A\n$direction_A\n$trimmed_quantity_A\n" if ($fmtq eq "fastq");
	print OT_A ">$query_A[1]\n$trimmed_sequence_A\n" if ($fmtq eq "fasta");
	print OT_B "\@$query_B[1]\n$trimmed_sequence_B\n$direction_B\n$trimmed_quantity_B\n" if ($fmtq eq "fastq");
	print OT_B ">$query_B[1]\n$trimmed_sequence_B\n" if ($fmtq eq "fasta");
	$S_PE_number_FQ1 ++;
	$S_PE_number_FQ2 ++;
	$S_PE_number_FQ1_size += (length $trimmed_sequence_B);
	$S_PE_number_FQ2_size += (length $trimmed_sequence_B);
	}
	elsif ( ($skip_A eq "T") and ($skip_B eq "F") ){
	if (exists $opt{clip}){
	if (length $trimmed_sequence_B > $clip){
	$trimmed_sequence_B = substr($trimmed_sequence_B,0,$clip);
	$trimmed_quantity_B = substr($trimmed_quantity_B,0,$clip);
	}
	}
	print OT_S "\@$query_B[1]\n$trimmed_sequence_B\n$direction_B\n$trimmed_quantity_B\n" if ($fmtq eq "fastq");
	print OT_S ">$query_B[1]\n$trimmed_sequence_B\n" if ($fmtq eq "fasta");
	$S_SE_number_FQ2 ++;
	$S_SE_number_FQ2_size += (length $trimmed_sequence_B);
	}
	elsif ( ($skip_A eq "F") and ($skip_B eq "T") ){
	if (exists $opt{clip}){
	if (length $trimmed_sequence_A > $clip){
	$trimmed_sequence_A = substr($trimmed_sequence_A,0,$clip);
	$trimmed_quantity_A = substr($trimmed_quantity_A,0,$clip);
	}
	}
	print OT_S "\@$query_A[1]\n$trimmed_sequence_A\n$direction_A\n$trimmed_quantity_A\n" if ($fmtq eq "fastq");
	print OT_S ">$query_A[1]\n$trimmed_sequence_A\n" if ($fmtq eq "fasta");
	$S_SE_number_FQ1 ++;
	$S_SE_number_FQ1_size += (length $trimmed_sequence_A);
	}
	
	}
	close (OT_A);
	close (OT_B);
	close (OT_S);
}
sub Write_SE{
my $file;
	$file = `file $_[0]`;
	if ($file =~ "gzip compressed data"){
	open IN,"gzip -dc $_[0] |" or die $!;
	}
	else{
	open IN,"$_[0]" or die $!;
	}
	if ($g == 1){
	open OT,"| gzip > $out.$fmtq.gz" or die $!;
	$outputfile = "$out.$fmtq.gz\n";
	}
	else{
	open OT,">$out.$fmtq" or die $!;
	$outputfile = "$out.$fmtq\n";
	}
	while (<IN>){
	chomp;
	$Reads_num ++;
	my $skip = "F";
	chomp;
	my $query = $_;
	my @query = split (/\@/,$query);
	my $sequence = <IN>;
	   $sequence =~ s/\n//;
	my $trimmed_sequence = $sequence;
	my $direction = <IN>;
	   $direction =~ s/\n//;
	my $quantity = <IN>;
	   $quantity =~ s/\n//;
	my $trimmed_quantity = $quantity; 
		if ( (defined $al) or (defined $hc) ){
		&CHECK_Q($query[1],$skip,$S_ac,$S_hcc);
		}
		if (defined $lqc){
		&CHECK_LQC($quantity,$skip,$S_lqc) if ($skip eq "F");
		}
		if (defined $cc){
		&CHECK_CC($sequence,$skip,$S_cc) if ($skip eq "F");
		}
		if (defined $amc){
		&CHECK_AMC($sequence,$skip,$S_amc) if ($skip eq "F");
		}
		if (defined $TBc){
		my $last_nucli = substr($trimmed_quantity,-1,1);
		my $last_nucli_quantity = ord($last_nucli) - 33;
		&TRIM($trimmed_sequence,$trimmed_quantity,$S_TBc,$SizeTrim) if ( ($skip eq "F") && ( ($trimmed_sequence =~ "N\$") || ($last_nucli_quantity <= $Q) ) );
		}
		if ($skip eq "F"){
		print OT "\@$query[1]\n$trimmed_sequence\n$direction\n$trimmed_quantity\n" if ($fmtq eq "fastq");
		print OT ">$query[1]\n$trimmed_sequence\n" if ($fmtq eq "fasta");
		}
	}
	close (IN);
}

sub CHECK_Q{
my ($a,$b,$c,$d) = @_;
my @query = split /\s+/,$a;
	if (exists $Pre_input_al{$query[0]}){
	$b = "T";
	$_[1] = $b;
	$c ++;
	}
	elsif (exists $Pre_input_hc{$query[0]}){
	$b = "T";
	$_[1] = $b;
	$d ++;
	}
$_[2] = $c;
$_[3] = $d;
}

sub CHECK_LQC{
my ($a,$b,$c) = @_;
my @Quantity = split (//,$a);
my $judge = 0;

	foreach my $m (@Quantity){
	my $value = ord($m) - 33;
		if ($value <= $Q){
		$judge ++;
		}
	}
 
	if ($judge >= $lqc){
	$b = "T";
	$_[1] = $b;
	$c ++;
	}
$_[2] = $c;
}
#&CheckIndex($query_A[2],$query_B[2],$skip_A,$skip_B,$S_index1,$S_index2);
sub CheckIndex{
my @queryA = split /\:/,$_[0];
my @queryB = split /\:/,$_[1];
	if ( ($queryA[-1] ne $sequencing_index) or ($queryB[-1] ne $sequencing_index) ){
	$_[2] = "T";
	$_[3] = "T";
	$_[4] ++;
	$_[5] ++;
	}
}

sub CHECK_CC{
my ($a,$b,$c) = @_;
	if ( ($a =~ ("A" x $cc)) || ($a =~ ("T" x $cc)) || ($a =~ ("G" x $cc)) || ($a =~ ("C" x $cc)) ){
	$b = "T";
	$_[1] = $b;
	$c ++;
	}
$_[2] = $c;
}

sub CHECK_AMC{
my ($a,$b,$c) = @_;
$a =~ s/A|T|G|C//g;
	if ( (length $a) >= $amc ){
	$b = "T";
	$_[1] = $b;
	$c ++;
	}
$_[2] = $c;
}

sub TRIM_HEAD{
my $sequence = $_[0];
$sequence =~ s/^N+//;
$_[0] = $sequence;
}

sub TRIM{
my ($a,$b,$c,$d) = @_;
my $judge_size;
	&Who_is_strict($b,$a);

	if ($TBc eq "f"){
	$_[3] += $fatter;
	$fatter = $length - $fatter;
	$a = substr ($a,0,$fatter);
	$b = substr ($b,0,$fatter);
	$_[0] = $a;
	$_[1] = $b;
	}
	else{
		if ($TBc eq "a"){
		$judge_size = $c;
		}
		else{
		$judge_size = $TBc;
		}
		if ($judge_size > $fatter){
		$judge_size = $fatter;
		}
		$_[3] += $judge_size;
		$judge_size = $length - $judge_size;
		$a = substr ($a,0,$judge_size);
		$b = substr ($b,0,$judge_size);
		$_[0] = $a;
		$_[1] = $b;
	}
}
