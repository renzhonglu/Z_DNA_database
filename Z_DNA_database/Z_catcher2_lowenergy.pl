#!/usr/bin/perl
use strict;
#use warnings;
##This programme scan the DNA sequence fasta format sequence file and return the results###
##Please Choose the negative supercoiling density cutoff value for the prediction###
##It could be -0.035, -0.04, -0.045, -0.05, -0.055, -0.06, -0.065###
##Usage of this script is same with Z_catcher.pl
#################################################
################Arguments mode###################
#################################################

my ($denscut,$inputfile,$outputfile);
my @arg = @ARGV;
if($#arg != -1)
{
	my $first = $arg[0];
	my $second = @arg;
	if( $first eq "-h" or $first eq "--help" ){
		print "\n##################################################################################################\n";
		print "This programme is used to predict the potential Z-DNA regions in given Fasta sequences.\n";
		print "The programme needs three arguments:\n";
		print "1. Sigma value: behind -s, it is negative supercoiling density always negative number.\n";
		print "   In here, it must be equal to or less negative than -0.07.\n";
		print "   It could be -0.035, -0.04, -0.045, -0.05, -0.055, -0.06, -0.065.\n";
		print "2. Input file name: behind -i, the file must be fasta format of sequences without n or N base and gap.\n";
		print "3. Output file name: behind -o, you must assign a name for output file with ZDRs results.\n";
		print "\n";
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
	###get Sigma###
	print "Please Choose the negative supercoiling density cutoff value for the prediction\n";
	print "It could be -0.035, -0.04, -0.045, -0.05, -0.055, -0.06, -0.065\n";
	print "Sigma is: ";
	$denscut = <STDIN>;
	chomp $denscut;
	###get fasta file name###
	print "Please key in the name of the file which contains sequences in FASTA format\n";
	print "File name is: ";
	$inputfile = <STDIN>;
	chomp $inputfile;
	###get output file name###
	print "Plase key in the name of the output file\n";
	print "Output file name is: ";
	$outputfile = <STDIN>;
	################################################
	&prediction($denscut,$inputfile,$outputfile);
}

#################################################
################Interaction End##################
#################################################


sub prediction{
	my ($denscut, $inputfile, $outputfile) = @_;
	my $now1 = time;
	print "...Please wait...\n";
	&getTempa($inputfile);
	#########OPEN goldensum.txt for the next analysis###
	#########Read in correspond value for analysis###
	########Values in hash structure named %lowdens###
	open (SUM, "goldensum.txt")||die ("goldensum.txt file does not exist\n");
	my ($header, $entry, @keyvalues, $minilength, $maxsegmentlen, %lowdens, $key, $value);
	#changing the separator from "\n" to ">"
	local $/ = ">";
	#Jump the first one
	my $line = <SUM>;
	while ( $line =<SUM> ) {
		$line =~ s/\r?\n>//;
		($header, $entry) = split "\n", $line, 2;
		if($header =~ /^$denscut$/ ) {
			@keyvalues = split "\n",$entry;
			for(my $a=0;$a<@keyvalues;$a++){
				my ($key,$value) = split "\t",$keyvalues[$a];
				$lowdens{$key} = $value;
				$minilength = $key if($a == 0);
				$maxsegmentlen = $key if($a == @keyvalues - 1)
			}
		}
		
		#foreach my $key  (sort keys(%lowdens)){
			#print "$key: $lowdens{$key}\n";
		#}
	}
	#Changing back defaults separator "\n"
	$/ = "\n"; 
	if($minilength eq "" or $maxsegmentlen eq ""){
		print "\n\nplease set the cutoff density as -0.035, -0.04, -0.045, -0.05, -0.055, -0.06, -0.065\n";
		exit;
	}
	close SUM;
	print "The file: goldensum.txt is read in Hash.\n";

	####Using TEMPA temporary file to generate tempb.txt
	print "...Please wait... tempb.txt\n";

	open (TEMPB, ">tempb.txt")||die ("tempb.txt file does not exist\n");
	open (SEQ, "tempa.txt")||die ("tempa.txt file does not exist\n");
	#my $sequence = <SEQ>;
	#chomp $sequence;
	while (my $sequence =<SEQ>){
		chomp $sequence;
		my @array = split(/\t/,$sequence);
		my $seqname = $array[0];
		my $position = $array[1];
		my $len = $array[2];
		#filtering ZDR length using $minilength in golden file
		if ($len < $minilength)
		{
			next;
			$sequence = <SEQ>;
			#chomp $sequence;
		}
		else
		{ 
			my $zseq = $array[3];
			my $segmentlen = $maxsegmentlen; 
			if ($segmentlen > $len)
			{
				$segmentlen = $len;
			}
			while ($segmentlen >= $minilength)
			{
				my $sumgmax = $lowdens{$segmentlen};
				my $stringlen = 0;
				while ($stringlen <= $len - $segmentlen)
				{
					my $segment = substr($zseq, $stringlen, $segmentlen);
					my $sign = &signas($segment, $segmentlen);
					my $sumg = &sumg($segment, $sign);
					if ($sumg > $sumgmax)
					{
						$stringlen++;
					}
					if ($sumg <= $sumgmax)
					{
						my $ntposition = $position + $stringlen;
						my $ntpositionend = $ntposition + $segmentlen - 1;
						print TEMPB "$seqname\t$ntposition\t$ntpositionend\t$segmentlen\t$segment\n";
						$stringlen += $segmentlen;
					}
				}
				$segmentlen = $segmentlen - 1;
			}
		#}#################################################
		#$sequence = <SEQ>;
		#chomp $sequence;
		}
	}
	close SEQ;
	close TEMPB;
	print "The temporary file: tempb.txt is generated. OK!\n";
	&getResult($outputfile, $denscut);

	my @tempfiles = ('tempa.txt','tempb.txt');
	unlink @tempfiles;

	my $now2 = time;
	my $diff = $now2-$now1;
	my $min = int($diff / 60);
	my $second = $diff % 60;
	print "Time consuming is $min min, $second s.\n";
}

###############################################################
sub getResult{
	my ($outputfile,$denscut) = @_;
	open (OUT, ">$outputfile")||die;
	open (SEQ, "tempb.txt")||die;
	#print OUT "The cutoff density is\t$denscut\n";
	print  "The cutoff density is\t$denscut\n";
	print OUT "SequenceName\tZDRstart\tZDRlength\tZDRsequence\n";
	my $file = <SEQ>;
	chomp $file;
	while($file ne '')
	{
		my $sequence1 = $file;
		my @array1 = split(/\t/, $sequence1);
		print OUT "$array1[0]\t$array1[1]\t$array1[3]\t$array1[4]\n";
		my $name1 = $array1[0];
		my $start1 = $array1[1];
		my $end1 = $array1[2];
		$file = <SEQ>;
		
		if ($file eq '')
		{
			last;
		}
		chomp $file;
		my $sequence2 = $file;
		my @array2 = split(/\t/, $sequence2);
		my $name2 = $array2[0];
		my $start2 = $array2[1];
		my $end2 = $array2[2];
		while(($start2 <= $end1)&&($name1 eq $name2))
		{
			$file = <SEQ>;
			
			if ($file eq '')
			{
				last;
			}
			chomp $file;
			$sequence2 = $file;
			@array2 = split(/\t/, $sequence2);
			$name2 = $array2[0];
			$start2 = $array2[1];
			$end2 = $array2[2];
		}
	}
	close SEQ;
	close OUT;
}


sub getTempa {
	my ($segment, $segmentend, $name, @namearray);
	#Here firstly created a file using -0.09 as the cutoff of ZDRs, and named tempa.txt#
	my $inputfile = $_[0];
	print "...Please wait... tempa.txt\n";
	open (SEQ, "$inputfile")||die("The query sequences file does not exist\n");
	open (TEMPA, ">tempa.txt");
	my $tempdenscut = 0 - 0.09;
	my $file = <SEQ>;
	chomp $file;
	while ($file ne '')
	{
		my @namearray = split(/\t/, $file);
		$name = $namearray[0];
		#$position is the starting point of this sequence in the genome, the position for a particular nt is $position + $stringlen#
		my $position = 1;
		my $sequence = '';
		$file = <SEQ>;
		chomp $file;
		$file = uc($file);
		#my $sequence;
		while(($file !~ />/)&&($file ne ''))
		{
			$sequence = $sequence.$file;
			$file = <SEQ>;
			chomp $file;
			$file = uc($file);
		}
		#####end pasting######
	
		#$len is the total length of this continuous sequence segment in the sequence#
		my $len = length($sequence);
		#the $stringlen is the position on the sequence itself start as $stringlen = 0#
		my $stringlen = 0;
	
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
			#$segmentlen is the length of the segment under consideration, it is just part of the continuous sequence segment#
			my $segmentlen = length($segment);
			my $sign = &signas($segment, $segmentlen);
			my $sumg = &sumg($segment, $sign);
			my @dens = &dens($sumg, $segmentlen);
			my $dens = $dens[0];
			my $minidens = $dens[1];
			if ($dens < $tempdenscut)
			{
				$stringlen++;
			}
			elsif ($dens >= $tempdenscut)
			{
				my $extent = 2;
				$segmentlen = 12 *$extent;
				$segmentend = $stringlen + $segmentlen;
				while (($dens >= $tempdenscut) &&  ($segmentend < $len+12))
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
				#$drawback is the nt cut from the extent segment, just to go back by 2 nt each time#
				my $drawback = 0;
				while ($dens < $tempdenscut && $drawback <= 12)
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
				my $ntposition = $stringlen + $position;
				print TEMPA "$name\t$ntposition\t$segmentlen\t$segment\n";
				$stringlen = $stringlen + $segmentlen;
			}
		}
	}
	close SEQ;
	close TEMPA;
	print "The temporary file: tempa.txt is generated. OK!\n";
}


###########################################################
###########################################################
sub prescreen
{
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
		elsif (substr($assign, $offset, 2) eq 'SA')
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


