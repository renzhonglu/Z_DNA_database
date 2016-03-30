#/usr/bin/perl
use strict;
#use warnings;
##Edition: 1.0
##Author: Ren Zhonglu
##Organization: Department of Bioinformatics, School of Basic Medical Sciences, Southern Medical University, Guangzhou, China PR
##This programme scans the whole Chromosome Fasta format sequence file and returns the potential ZDRs results##
#Script usage:
#1. perl Zcatcher_genome.pl -i chrY.fa -o chrY_007_out.txt -s -0.07
#2. perl Zcatcher_genome.pl

#################################################
################Arguments mode###################
#################################################
my $denscut = -0.07;
my $inputfile = "";
my $outputfile = "out.txt";
my @arg = @ARGV;

if($#arg != -1)
{
	my $first = $arg[0];
	my $second = @arg;
	if($first eq "-h" or $first eq "--help"){
			print "##################################################################################################\n";
			print "This programme is used to predict the potential Z-DNA regions in given chromosome.\n";
			print "The programme needs three arguments:\n";
			print "1. Sigma0 value: behind -s, it is negative supercoiling density always negative number.\n";
			print "   In here, it must be equal to or more negative than -0.07, and default value is -0.07.\n";
			print "2. Input file name: behind -i, the file must be fasta format of chromosome, one chromosome one file.\n";
			print "3. Output file name: behind -o, you must assign a name for output file with ZDRs results. The default name is out.txt.\n";
			print "You can also use this programme without any arguments, then you will key in arguments following programme.\n";
			print "Then get it!\n";
			print "##################################################################################################\n";
		}
		else{
			my $arg_len = scalar(@arg);
			for(my $j=0; $j < $arg_len/2; $j++)
			{
				$first = shift @arg;
				$second = shift @arg;
				if($second ne ""){
					if($first eq "-i"){
						$inputfile = $second;
					}
					elsif($first eq "-o"){
						$outputfile = $second;
					}
					elsif($first eq "-s"){
						$denscut = $second;
					}
					else{
						print "Please key in correct arguments\n";
						last;
					}
				}
				else{
					print "Please key in correct arguments\n";
					last;
				}
			}
			&prediction($denscut,$inputfile,$outputfile);
		}
}
else{
	print "Using interactive mode.\n";
	print "Please key in the negative supercoiling density cutoff value for the prediction\n";
	print "It must be equal to or more negative than -0.07, otherwise please use Z-catcher_lowdensity.pl\n";
	print "Sigma0 is: ";
	$denscut = <STDIN>;
	chomp $denscut;
	print "Please key in the name of the file which contains a chromosome sequence in FASTA format\n";
	print "File name is: ";
	$inputfile = <STDIN>;
	chomp $inputfile;
	print "Plase key in the name of the output file\n";
	print "Output file name is: ";
	$outputfile = <STDIN>;
	&prediction($denscut,$inputfile,$outputfile);
}

#################################################
################Interaction End##################
#################################################



####################################################
###########First subroutine#########################
sub prediction {
	my $now1 = time;
	my ($denscut,$inputfile,$outputfile) = @_;
	my ($sumg,$dens,@dens,$segment,$segmentlen,$sign,$extent,$minidens,$drawback,$ntposition,$endposition,$segmentend,$stringlen,$position);
	##Deleting the N or n base which is in the fasta file and 
	##noting the start and end position of non-N sequences
	open (OUTONE, ">nonN_seq.txt")||die;
	open (CHRSEQ, "$inputfile")||die("The query sequence file does not exist\n");
	print "...Please wait...\n";
	#reading fasta file per line
	my $line = <CHRSEQ>;
	my $name = $line;
	chomp $name;
	my $seq = "";
	while ($line=<CHRSEQ>){
		chomp $line;
		$line = uc($line); #upper case
		$seq = $seq.$line;
	}
	close(CHRSEQ);

	my $seq_len = length($seq);
	my (@sta_pos,@end_pos,$nonn_end,$nonn_sta,$nonn_len,$sub_seq,$len,$end,$start);
	if($seq =~ /(N+)/i) {
		while ($seq =~ m/(N+)/ig) { 	#match one or more (n|N)
			$len = length($1); 	#return match length
			$end = pos($seq); 	# pos() to get the end position of this match
			$start = $end - $len + 1; 	#obtain start position of this match
			#print "$id\t$start\t$end\t$len\n"; 	#print results
			push @sta_pos, $start;
			push @end_pos, $end;
		}
		if($sta_pos[0] eq 1){
			for(my $i = 0;$i < $#sta_pos;$i++){
				$nonn_sta = $end_pos[$i] + 1;
				$nonn_end = $sta_pos[$i+1] - 1;
				$nonn_len = $nonn_end - $nonn_sta + 1;
				$sub_seq = substr($seq,$nonn_sta-1,$nonn_len);
				print OUTONE "$name\t$nonn_sta\n";
				print OUTONE "$sub_seq\n";
			}
			#last non-n seq position
			$nonn_sta = $end_pos[$#end_pos]+1;
			if($nonn_sta < $seq_len){
				#$nonn_len = $seq_len - $nonn_sta + 1;
				$sub_seq = substr($seq,$nonn_sta-1);
				print OUTONE "$name\t$nonn_sta\n";
				print OUTONE "$sub_seq\n";
			}
		}
		else{
			print OUTONE "$name\t1\n";
			$nonn_end = $sta_pos[0] - 1;
			$sub_seq = substr($seq,0,$nonn_end-1);
			print OUTONE "$sub_seq\n";
			for(my $i = 0;$i < $#sta_pos;$i++){
				$nonn_sta = $end_pos[$i] + 1;
				$nonn_end = $sta_pos[$i+1] - 1;
				$nonn_len = $nonn_end - $nonn_sta + 1;
				$sub_seq = substr($seq,$nonn_sta-1,$nonn_len);
				print OUTONE "$name\t$nonn_sta\n";
				print OUTONE "$sub_seq\n";
			}
			#last non-n seq position
			$nonn_sta = $end_pos[$#end_pos]+1;
			if($nonn_sta < $seq_len){
				$sub_seq = substr($seq,$nonn_sta-1);
				print OUTONE "$name\t$nonn_sta\n";
				print OUTONE "$sub_seq\n";
			}
		}
	}
	else{
		print OUTONE "$name\t1\n";
		print OUTONE "$seq\n";
	}
	close(OUTONE);
	print "The temporary file nonN_seq.txt is done...OK!\n";
	print "Please wait for the further analysis...\n";
	####nonN_seq.txt is temporary file for the next analysis####

	#############################################################################################
	##This is the result files that contains the predicted ZDRs under a certain sigma threshold##
	open (OUT, ">$outputfile");
	print OUT "cutoff supercoiling density is\t$denscut\n";
	print OUT "Chromosome\tID\tZDRstart\tZDRend\tZDRlength\tZDRsequence\n";
	
	open (SEQ, "nonN_seq.txt")||die;
	my $number = 1;
	while (my $file = <SEQ>)
	{
		chomp $file;
		if($file =~ />/){
			my @namearray = split(/\t/,$file);
			$name = uc($namearray[0]);
			$position = $namearray[1];
		}
		else{
			my $sequence = $file;
			##$len is the total length of this continuous sequence segment in the genome#
			$len = length($sequence);
			##the $stringlen is the position on the sequence itself start as $stringlen = 0#
			$stringlen = 0;
			while ($stringlen < $len - 12)
			{
				my $prescreen = 'impossible';
				while(($prescreen eq 'impossible') && ($stringlen < $len - 12))
				{
					$segment = substr($sequence, $stringlen, 12);
					$prescreen = &prescreen($segment);
					$stringlen++;
				}
			
				$stringlen = $stringlen - 1;
				##$segmentlen is the length of the segment under consideration, it is just part of the continuous sequence segment#
				$segmentlen = length($segment);
				$sign = &signas($segment, $segmentlen);
				$sumg = &sumg($segment, $sign);
				@dens = &dens($sumg, $segmentlen);
				$dens = $dens[0];
				$minidens = $dens[1];
				if ($dens < $denscut)
				{
					$stringlen++;
				}
				elsif ($dens >= $denscut)
				{
					$extent = 2;
					$segmentlen = 12 *$extent;
					$segmentend = $stringlen + $segmentlen;
					while (($dens >= $denscut) &&  ($segmentend < $len+12))
					{
						$segment = substr($sequence, $stringlen, $segmentlen);
						$segmentlen = length($segment);
						$sign = &signas($segment, $segmentlen);
						$sumg = &sumg($segment, $sign);
						@dens = &dens($sumg, $segmentlen);
						$dens = $dens[0];
						$minidens = $dens[1];
						$extent++;
						$segmentlen = 12 *$extent;
						$segmentend = $stringlen + $segmentlen;
					}
					$segmentlen = length($segment);
					##$drawback is the nt cut from the extent segment, just to go back by 2 nt each time#
					$drawback = 0;
					while (($dens < $denscut) && ($drawback <= 12))
					{
						$segmentlen = $segmentlen - 2;
						$segment = substr($sequence, $stringlen, $segmentlen);
						$sign = &signas($segment, $segmentlen);
						$sumg = &sumg($segment, $sign);
						@dens = &dens($sumg, $segmentlen);
						$dens = $dens[0];
						$minidens = $dens[1];
						$drawback += 2;
					}
					#ZDRstart position
					$ntposition = $stringlen + $position;
					#ZDRend position
					$endposition = $ntposition + $segmentlen - 1;
					print OUT "$name\t$number\t$ntposition\t$endposition\t$segmentlen\t$segment\n";
					$stringlen = $stringlen + $segmentlen;
					$number++;
				}
			}
		}
	}
	close SEQ;
	close OUT;

	#deleting temporary files
	my @temparray = ('nonN_seq.txt');
	unlink @temparray;

	#Getting time consuming
	my $now2 = time;
	my $diff = $now2-$now1;
	my $min = int($diff / 60);
	my $second = $diff % 60;
	print "The Programme is done.\n";
	print "$name: Time consuming is $min min, $second s.\n";
}

###############firtst subroutine end###################
#######################################################

sub prescreen
{
	#read the first argument in $seq
	my $seq = $_[0];
	my $translate = $seq;
	my $seqlength = length($seq);
	$translate =~ tr/ATGC/SASA/;
	my $prebase = chop($translate);
	my $score = 0;
	my  $threshold = $seqlength / 2;
	for(my $i=1; $i< $seqlength; $i++)
	{
			my $nextbase = chop($translate);
			if ($prebase ne $nextbase)
			{
				$score++;
			}
			$prebase = $nextbase;
	}
	
	if ($score >= $threshold)
	{
		return $seq;
	}
	if ($score < $threshold)
	{
		return 'impossible';
	}
}



sub signas
{
	my ($seq, $len)= @_;
	my $translate = $seq;
	$translate =~ tr/ATGC/SASA/;
	$translate =~ s/AAAAA/ASASA/g;
	$translate =~ s/SSSSS/SASAS/g;
	$translate =~ s/AAA/ASA/g;
	$translate =~ s/SSS/SAS/g;
	$translate =~ s/AASS/ASAS/g;
	$translate =~ s/SSAA/SASA/g;
	my $rvtrans = reverse $translate;
	my $prebase = chop($rvtrans);
	my $sign = $prebase;
	for(my $i=1; $i< $len; $i++)
	{
		my $nextbase = chop($rvtrans);
		if ($prebase ne $nextbase)
		{
			$sign = $sign."$nextbase";
		}
		if ($prebase eq $nextbase)
		{
			$sign = $sign.'Z';
		}
		$prebase = $nextbase;
	}
	return $sign;
}


sub sumg
{

	my %hashAS = ('CG' => 0.7,	'CA' => 1.3,	'CC' => 2.4,	'CT' => 3.4,
				'GC' => 4.0,	'GT' => 4.6,	'GG' => 2.4,	'GA' => 3.4,
				'AC' => 4.6, 'AT' => 5.9,	'AG' => 3.4,	'AA' => 3.9,
				'TG' => 1.3, 'TA' => 2.5,	'TC' => 3.4,	'TT' => 3.9);

	my %hashSA = ('CG' => 4.0,	'CA' => 4.6,	'CC' => 2.4,	'CT' => 3.4,
				'GC' => 0.7,	'GT' => 1.3,	'GG' => 2.4,	'GA' => 3.4,
				'AC' => 1.3, 'AT' => 2.5,	'AG' => 3.4,	'AA' => 3.9,
				'TG' => 4.6, 'TA' => 5.9,	'TC' => 3.4,	'TT' => 3.9);

	my %hashZZ = ('CG' => 4.0,	'CA' => 4.5,	'CC' => 4.0,	'CT' => 6.3,
				'GC' => 4.0,	'GT' => 4.5,	'GG' => 4.0,	'GA' => 6.3,
				'AC' => 4.5, 'AT' => 5.6,	'AG' => 6.3,	'AA' => 7.4,
				'TG' => 4.5, 'TA' => 5.6,	'TC' => 6.3,	'TT' => 7.4);

	my($seq, $assign) = @_;
	my $strlen = length($seq);
	my $dinucleotides = int($strlen/2);
	my $sumg = 0;
	for (my $i=0; $i< $dinucleotides; $i++)
	{
		my $offset = $i *2;
		if (substr($assign, $offset, 2) eq 'AS')
		{
			$sumg += $hashAS{substr($seq, $offset, 2)}; 
		}
		elsif (substr($assign, $offset, 2) eq 'SA')  ######2014.12.22########
		{
			$sumg += $hashSA{substr($seq, $offset, 2)}; 
		}
		else
		{
			$sumg += $hashZZ{substr($seq, $offset, 2)};
		}
	}
	return $sumg;
}

sub dens
{
	my ($sumg, $m) = @_;
	$sumg = $sumg + 10;
	$sumg = $sumg * 1000;
	my $a = 19957.428;
	my $c = $sumg * 4236;
	my $amm = $a * $m * $m;
	my $dens = 0 - $c - $amm;
	$dens = $dens/ 93932961.12;
	$dens = $dens / $m;
	my $minidens= sqrt($sumg/26092489.2);
	$minidens = 0 - $minidens;
	return my @array = ($dens,$minidens);
}


