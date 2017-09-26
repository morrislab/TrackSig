#!/usr/bin/perl
# Make mutation counts by 100 from vcf and phi file. Also script dumps the ordered mutations into MutOrderFile (useful to reorder mutation assignments in the same order).

# AUTHOR: Amit Deshwar

use strict;
use warnings;
use POSIX;
use Bio::DB::Fasta; # Library for parsing Fasta files
use List::Util qw(sum);
sub mean {
    return sum(@_)/@_;
}
sub min ($$) { $_[$_[0] > $_[1]] }

my %fastaHandles = ("path" => "annotation/hg19/"); # human reference genome hg19

###############################################################
my $caller = $ARGV[0]; # caller definition: muse/caveman/dkfz..
my $vcfFile = $ARGV[1];
my $phiFile = $ARGV[2];
my $MutOrderFile = $ARGV[3];
###############################################################

# Assigning 96 trinucleotide context into a hash
my $line;
my %counts;
my $countsFile = "annotation/trinucleotide.txt";
open (COUNTSFILE, $countsFile) or die "Couldn't open $countsFile\n";
while ($line = <COUNTSFILE>) # initializing 96 trinucleotide counts to 0
{
	chomp $line;
	$counts{"$line"} = 0;
}
close COUNTSFILE;

# Leave only mutations that are present in vcf file
my @filtered_mutations;
open (VCFFILE, $vcfFile) or die "Couldn't open $vcfFile\n";
while ($line = <VCFFILE>)
{
        if ($line !~ /^#/)
        {
                chomp $line;
                my @fields = split(/\t/, $line);
                my $chr = $fields[0];
                if ($chr  =~ /^chr[\d]*/)
                {
                        $chr = substr $chr, 3;
                }
                my $pos = $fields[1];

                push @filtered_mutations, "$chr\t$pos";

        }
}
close VCFFILE;
my %filtered_mutations_hash = map { $_ => 1 } @filtered_mutations;


# Ordering the phi file
my %counts_phi; # hash for the specific mutation with the phi value as the key
my $mutation_count=0; # stores the mutation counts in the phi file
open (PHIFILE, $phiFile) or die "Couldn't open $phiFile\n";
while ($line = <PHIFILE>)
{
	chomp $line;
	my $char1 = index($line, '_');
	my $char2 = index($line, '=');
	my $chr_phi = substr $line, 0, $char1;
	my $pos_phi = substr $line, $char1+1, $char2-$char1-1;
	my $phi = substr $line, $char2+1;
	if (exists($filtered_mutations_hash{"$chr_phi\t$pos_phi"}))
	{	
		$counts_phi{"$chr_phi\t$pos_phi"} = $phi;
		$mutation_count++
	}
}
close PHIFILE;
# Sorting the phis of mutations in decreasing order
my @phi_array = sort {$counts_phi{$b} <=> $counts_phi{$a}} keys %counts_phi;


open(my $fh, '>', $MutOrderFile) or die "Could not open file '$MutOrderFile' $!";
print $fh "$_\n" for @phi_array;
close $fh;

my @phi_values = $counts_phi{@phi_array};

# Reading the vcf file
my $chr;
my $start;
my $end;
my $fastaHandles;
my $seq;
my %mutation_types;
open (VCFFILE, $vcfFile) or die "Couldn't open $vcfFile\n";
while ($line = <VCFFILE>)
{
	if ($line !~ /^#/)
	{
		chomp $line;
		my @fields = split(/\t/, $line);
		$chr = $fields[0];
		if ($chr  =~ /^chr[\d]*/)
		{
			$chr = substr $chr, 3;
		}

		my $pos = $fields[1];
		$start = $pos - 1;
		$end = $pos + 1;
		my $ref = uc($fields[3]);
		my $alt = uc($fields[4]);
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
    	elsif ($chr =~ /^\d+?$/)
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
			elsif (length($alt) == 1) # 1 alternative nucleotide is present
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
						$mutation_types{"$chr\t$pos"} = "$ref\t$alt1\t$seq";
						# $counts{"$ref\t$alt1\t$seq"}++;
						# $counts{"$ref\t$alt2\t$seq"}++;
					} 
					else
					{
						$mutation_types{"$chr\t$pos"} = "$ref\t$alt\t$seq";
					}
				}
			}
		}
	}
}
close VCFFILE;

for my $variant (@phi_array)
{
	if (exists $mutation_types{$variant}) {
        	print "$variant\t$counts_phi{$variant}\t$mutation_types{$variant}\n";
	}
}

# Function that retrieves the sequence from the human reference genome hg19
sub getSequence
{
	my $chr_ = "chr".$chr;
	$start = $_[1];
	$end = $_[2];
	$fastaHandles = $_[3];

	unless (exists ($fastaHandles->{$chr_})) # create a handle for the chromosome fasta, if it doesn't exist
	{
		$fastaHandles->{$chr_} = Bio::DB::Fasta->new("$fastaHandles{path}/$chr_.fa");
	}

	$seq = $fastaHandles{$chr_}->seq($chr_, $start, $end);
	return $seq;
}
