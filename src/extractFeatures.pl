#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use Bio::DB::Fasta; # Library for parsing Fasta files
my %fastaHandles = ("path" => "/users/nilsahin/Desktop/morris/hg19"); # human reference genome hg19

###############################################################
my $caller = $ARGV[0]; # caller definition: muse/caveman/dkfz..
my $vcfFile = $ARGV[1];
###############################################################

# Assigning 96 trinucleotide context into a hash
my $line;
my %counts;
my $countsFile = "trinucleotide.txt";
open (COUNTSFILE, $countsFile) or die "Couldn't open $countsFile\n";
while ($line = <COUNTSFILE>) # initializing 96 trinucleotide counts to 0
{
	chomp $line;
	$counts{"$line"} = 0;
}
close COUNTSFILE;

# Reading the vcf file
my $chr;
my $start;
my $end;
my $fastaHandles;
my $seq;
open (VCFFILE, $vcfFile) or die "Couldn't open $vcfFile\n";
while ($line = <VCFFILE>)
{
	if ($line !~ /^#/)
	{
		chomp $line;
		my @fields = split(/\t/, $line);
		$chr = $fields[0];
		my $pos = $fields[1];
		$start = $pos - 1;
		$end = $pos + 1;
		my $ref = $fields[3];
		my $alt = $fields[4];
		# used if 2 alternative nucleotides are present such as G,A in ALT column in the vcf File
		my $alt1;
        my $alt2;
        
        # Check on the accepted mutations in the vcf file as 
        # Mutations found by the muse caller are all passed
        my $checkPass=0;
		if (($caller eq "muse") or ($fields[6] eq "PASS"))
		{
			$checkPass++;
		}
        
        # Check on the chromosome containing the particular mutation is in the human reference genome hg19
        my $checkChr=0;
    	if (($chr eq 'M') or ($chr eq 'X') or ($chr eq 'Y'))
    	{
    		$checkChr++;
    	}
    	elsif (isdigit($chr))
    	{
    		if (($chr > 0) and ($chr < 23))
    		{
    			$checkChr++;
    		}
    	}
    	
    	# The checks should be passed
		# The reference should be a single nucleotide for the particular mutation to count as an SNV
		if ((length($ref) == 1) and ($checkChr == 1) and ($checkPass == 1))
		{
			my $checkSNP=0; # Check on the SNP conditions are passed
			my $checkGA=0; # Check on the nucleotide as a Guanine or Adenine for further editing
			
			my $comma = index($alt,',');
			if ($comma != -1) # 2 alternative nucleotides are present such as G,A in ALT column in the vcf File
			{
				$checkSNP++;
               	$alt1 = substr $alt, $comma-1, 1;
        		$alt2 = substr $alt, $comma+1, 1;
				if ($ref =~ s/G/C/) # The complementary of the reference nucleotide is entered: G -> C
       			{
       				# The complementary of the alternative nucleotides are entered
       				if ($alt1 =~ s/T/A/) {}
               		elsif ($alt1 =~ s/C/G/) {}
                	elsif ($alt1 =~ s/A/T/) {}
               		if ($alt2 =~ s/T/A/) {}
                	elsif ($alt2 =~ s/C/G/) {}
                	elsif ($alt2 =~ s/A/T/) {}
					$checkGA++;
                }
                elsif ($ref =~ s/A/T/) # The complementary of the reference nucleotide is entered: A -> T
                {
                	# The complementary of the alternative nucleotides are entered
                	if ($alt1 =~ s/T/A/) {}
                	elsif ($alt1 =~ s/C/G/) {}
                	elsif ($alt1 =~ s/G/C/) {}
                	if ($alt2 =~ s/T/A/) {}
                   	elsif ($alt2 =~ s/C/G/) {}
       	 	        elsif ($alt2 =~ s/G/C/) {}
        		   	$checkGA++;
				}
			}
			elsif (length($alt) == 1)
			{
				$checkSNP++;
				if ($ref =~ s/G/C/) # The complementary of the reference nucleotide is entered: G -> C
				{
					# The complementary of the alternative nucleotide is entered
					if ($alt =~ s/T/A/) {}
					elsif ($alt =~ s/C/G/) {}
					elsif ($alt =~ s/A/T/) {}
					$checkGA++;
				}
				elsif ($ref =~ s/A/T/) # The complementary of the reference nucleotide is entered: A -> T
                {
                	# The complementary of the alternative nucleotide is entered
                    if ($alt =~ s/T/A/) {}
                    elsif ($alt =~ s/C/G/) {}
                    elsif ($alt =~ s/G/C/) {}
					$checkGA++;
                }
			}
			if ($checkSNP == 1)
			{
				$seq = getSequence($chr, $start, $end, \%fastaHandles); # Trinucleotide context is found
				$seq = uc $seq; # Uppercase
       			my $first = substr $seq, 0, 1; # 5' nucleotide of the trinucleotide context
        		my $third = substr $seq, 2, 1; # 3' nucleotide of the trinucleotide context
        		if ($checkGA == 1) # The complementary of the trinucleotide context sequence is entered
				{
					if ($first =~ s/G/C/) {}
                	elsif ($first =~ s/C/G/) {}
        			elsif ($first =~ s/A/T/) {}
         			elsif ($first =~ s/T/A/) {}
          			if ($third =~ s/G/C/) {}
       		 		elsif ($third =~ s/C/G/) {}
           			elsif ($third =~ s/A/T/) {}
           			elsif ($third =~ s/T/A/) {}
           			$seq = $third.''.$ref.''.$first; # The sequence is updated if the complementary is entered
      			}
      			my $NinSeq = index($seq,'N'); # The sequences containing N for a nucleotide are ignored
				if ($NinSeq == -1)
				{
					# trinucleotide counts are incremented
					if ($comma != -1)
					{
						$counts{"$ref\t$alt1\t$seq"}++;
						$counts{"$ref\t$alt2\t$seq"}++;
					}
					else
					{
						$counts{"$ref\t$alt\t$seq"}++;
					}
				}
			}
		}
	}
}
close VCFFILE;

# The vcf file and the trinucleotide counts are outputted
print $vcfFile;
for my $variant (sort keys %counts)
{
	print "\t$counts{$variant}";
}
print "\n";

# Function that retrieves the sequence from the human reference genome hg19
sub getSequence
{
	$chr = "chr".$_[0];
	$start = $_[1];
	$end = $_[2];
	$fastaHandles = $_[3];

	unless (exists ($fastaHandles->{$chr})) # create a handle for the chromosome fasta, if it doesn't exist
	{
		$fastaHandles->{$chr} = Bio::DB::Fasta->new("$fastaHandles{path}/$chr.fa");
	}

	$seq = $fastaHandles{$chr}->seq($chr, $start, $end);
	return $seq;
}
