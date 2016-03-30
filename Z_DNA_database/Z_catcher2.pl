#!/usr/bin/perl
use strict;
#use warnings;
##This programme scan the sequences in fasta format and then return the results##
#2012.12.26
#Usage: 1. perl Z_catcher.pl -i inputfile -o outputfile -s -0.07
#       2. perl Z_catcher.pl
#################################################
################Arguments mode###################
#################################################
my $denscut = -0.07;
my $inputfile ="";
my $outputfile ="out.txt";
my @arg = @ARGV;

if($#arg != -1)
{
	my $first = $arg[0];
	my $second = @arg;
	if($first eq "-h" or $first eq "--help"){
		print "\n##################################################################################################\n";
		print "This programme is used to predict the potential Z-DNA regions in given Fasta sequences.\n";
		print "The programme needs three arguments:\n";
		print "1. Sigma0 value: behind -s, it is negative supercoiling density always negative number.\n";
		print "   In here, it must be equal to or more negative than -0.07, and default value is -0.07.\n";
		print "2. Input file name: behind -i, the file must be FASTA format of sequences without n or N base and gap.\n";
		print "3. Output file name: behind -o, you must assign a name for output file with ZDRs results. The default output filename is Sigma-0.07_out.txt\n";
		print "You can also use this programme without any arguments, then you will key in arguments following programme.\n";
		print "Then get it!\n";
		print "##################################################################################################\n\n";
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

###############################################
###############Main subroutine################
###############################################
sub prediction {
	print "Please wait......\n";
	my ($denscut,$inputfile,$outputfile) = @_;
	#my $outputfilename = "Sigma".$denscut."_".$outputfile;
	my ($sumg,$dens,@dens,$segment,$segmentlen,$sign,$extent,$minidens,$drawback,$ntposition,$endposition,$segmentend);
	open(OUT,">$outputfile") || die "Please key in the output file name!";
	##the inputfile should be the name of the sequence that need to be predicted##
	open (SEQ, "$inputfile")||die("The query sequences file does not exist\n");
	#print OUT "cutoff supercoiling density is\t$denscut\n";
	print OUT "Sequence\tName\tZDRstart\tZDRend\tZDRlength\tZDRsequence\n";
	print "cutoff supercoiling density is\t $denscut\n";
	#print "Sequence\tName\tZDRstart\tZDRend\tZDRlength\tZDRsequence\n";
	

	#read query sequences file per line
	my $file = <SEQ>;
	chomp $file;
	while ($file ne '')
	{
		my $now1 = time;
		#upper case bases and fasta head
		$file = uc($file);
		my @namearray = split(/\t/, $file);
		my $name = $namearray[0];
		##$position is the starting point of this sequence in the genome, the position for a particular nt is $position + $stringlen##
		my $position = 1;
		my $sequence = '';
		$file = <SEQ>;
		chomp $file;
		#upper case bases
		$file = uc($file);
		
		my $sequence;
		#pasting multiple lines sequences into one line
		while(($file !~ />/)&&($file ne ''))
		{
			$sequence = $sequence.$file;
			$file = <SEQ>;
			chomp $file;
			$file = uc($file);
		}
		#####end pasting####

		##$len is the total length of this continuous sequence segment in the genome##
		##getting total length of the sequence
		my $len = length($sequence);
		##the $stringlen is the position on the sequence itself start as $stringlen = 0##
		##setting the start=0
		my $stringlen = 0;

		my $number = 1;
		while ($stringlen < $len - 12)
		{
			my $prescreen = 'impossible';
			while(($prescreen eq 'impossible') && ($stringlen < $len - 12))
			{
				#extracting the first 12bp segment for pre-screening
				$segment = substr($sequence, $stringlen, 12);
				$prescreen = &prescreen($segment);
				$stringlen++;
			}
			$stringlen = $stringlen - 1;
			##$segmentlen is the length of the segment under consideration, 
			##it is just part of the continuous sequence segment
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
				my $extent = 2;
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
				##$drawback is the nt cut from the extent segment, just to go back by 2 nt each time##
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
				#print "$name\t$number\t$ntposition\t$endposition\t$segmentlen\t$segment\n";
				$stringlen = $stringlen + $segmentlen;
				$number++;
			}
		}
		my $now2 = time;
		my $diff = $now2-$now1;
		my $min = int($diff / 60);
		my $second = $diff % 60;
		print "$name: Time consuming is $min min, $second s.\n";

	}
	close SEQ;
	close OUT;
	print "The Program is done.\n";
}


sub prescreen
{
	#read the first argument in $seq
	my $seq = $_[0];
	my $translate = $seq;
	my $seqlength = length($seq);
	$translate =~ tr/ATGC/SASA/;
	my $prebase = chop($translate);
	my $score = 0;
	my $threshold = $seqlength / 2;
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

	my($seq, $assign)= @_;
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
		if (substr($assign, $offset, 2) eq 'SA')
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
	my ($sumg, $m) =@_;
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


