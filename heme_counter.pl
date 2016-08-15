#!/usr/bin/perl
use strict;
use Getopt::Std;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Count CXXCH and CXXXCH motifs in each entry of an amino acid fasta file.
#Reminder: Output will be three columns: header to space, CXXCH count, CXXXCH count. The header on your fasta file should be gene_id followed by a space, such as what is exported by IMG: >2264878274 A35ODRAFT_02484 4-hydroxy-...

# - - - - - U S E R    V A R I A B L E S - - - - - - - -


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %options=();
getopts("f:h", \%options);
my %Sequences;

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
        print "-f = amino acid fasta file\n";
        print "-h = This help message\n\n";
	die;
    }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
&FASTAread($options{f}, 1);

foreach my $i (sort keys %Sequences)
	{	my $aminoseq = $Sequences{$i}{'gappy-ntseq'};
		my $aminoseq2 = $aminoseq;
		$Sequences{$i}{'CXXCH'} = $aminoseq =~ s/C..CH/C..CH/g;
		$Sequences{$i}{'CXXXCH'} = $aminoseq2 =~ s/C...CH/C...CH/g;
	}
foreach my $k (sort keys %Sequences)
	{	unless ($Sequences{$k}{'CXXCH'} >= 1)
		{$Sequences{$k}{'CXXCH'} = 0;}
		unless ($Sequences{$k}{'CXXXCH'} >= 1)
		{$Sequences{$k}{'CXXXCH'} = 0;}
	}

open(OUT1, ">$options{f}"."_hememotifs.txt");
print OUT1 "gene_id\tCXXCH_count\tCXXXCH_count\n";
foreach my $j (sort keys %Sequences)
	{	print OUT1 "$Sequences{$j}{'SHORT_HEAD'}\t$Sequences{$j}{'CXXCH'}\t$Sequences{$j}{'CXXXCH'}\n";
	}


close(OUT1);


print "\n\n* * * * * * * * D O N E * * * * * * * *\n\n";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub FASTAread
{	print "   Reading file . . . \n";
	# 1. Load FIlE . . . . . . . . . .
	$/=">";                                     # set input break string
	my $infile = $_[0];
        my $filenumber = $_[1];
	open(IN, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
	my @DATA = <IN>; close(IN); shift(@DATA);	
	# 2. Parse sequence data . . . . . . . . . . . . .
	my $unid = $filenumber.10001;                           # string to generate unique ids
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		my $seq = '';
		foreach my $i (1..$#data)
		{	$seq .= $data[$i];  }
		$seq =~ s/>//;
		$Sequences{$unid}{'HEAD'}    = $data[0];       # store header
		my @shorthead = split(' ', $data[0]);
		$Sequences{$unid}{'SHORT_HEAD'} = $shorthead[0];
		$Sequences{$unid}{'gappy-ntseq'}   = uc($seq);       # store aligned sequence
		$Sequences{$unid}{'SIZE'}    = length($seq);   # store length
		$seq =~ s/\.//;
                $seq =~ s/\-//;
                $Sequences{$unid}{'degapped-ntseq'} = uc($seq);     # store degapped sequence
                $Sequences{$unid}{'filenumber'} = $filenumber;
                $unid += 1;
	}
	$/="\n";
}
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
