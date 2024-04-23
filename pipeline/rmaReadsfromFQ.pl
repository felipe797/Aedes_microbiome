### Description: Extract FASTQ reads taking sequence headers as inputs

### Usage: <this script> <header file as input> <read1.fq> <read2.fq>
### Example: perl rmaReadsFromFQ.pl elizabethkingia.fasta reads_r1.fastq reads_r2.fastq

#!usr/bin/perl
use strict;
use IO::Uncompress::Bunzip2 '$Bunzip2Error';
use Cwd;

# accept input files required from the command line
my $infile = $ARGV[0]; chomp $infile;	# headers list file (e.g. elizabethkingia.fasta)
my $r1fq = $ARGV[1]; chomp $r1fq;	# R1 file
my $r2fq = $ARGV[2]; chomp $r2fq;	# R2 file
my %heads = ();
print "$infile - $r1fq - $r2fq\n";
# read headers file
open IN,"$infile" or die "Cannot open $infile!\n";
while(<IN>)
{
	chomp;
	next if(/^\s*$/);
	if(/^>/)
	{
		# changing the '>' symbol to '@'. this is to match with the header from original FASTQ file.
		$_=~s/>/@/; #>D00480L:369:H555TBCX2:2:1101:9836:96824/1
		
		# store the headers (only) in a hash array
		$heads{$_} = 0;
	}
}
close IN;
print scalar keys (%heads)."\n"; # if you want to count total number of headers, uncomment this line
#foreach my $x(keys %heads)
#{ print $x."\n"; }

# call subroutine to search for the headers in R1 and R2 files
# if a match is found then an output FASTQ file with the corresponding reads (sequence, quality scores) will be created
&geth($r1fq);
&geth($r2fq);

sub geth
{
	my $fq = shift;
	my @temp1 = split("\/",$fq,999);
	my $name = $1 if($temp1[8]=~/(.*)\.+bz2$/);
	my $outdir = Cwd::cwd();
	my $out = "$outdir\/$name.reads";
	open OUT,">$out" or die "Cannot open file for writing!\n";
print "$name\n$outdir\n$out\n"; 	
	my $ZF = IO::Uncompress::Bunzip2->new( $fq, {
	    AutoClose   => 1,
	    Transparent => 1,
	} ) or die "IO::Uncompress::Bunzip2 failed: $Bunzip2Error\n";

	while(my $line1 = <$ZF>)
	{
		chomp($line1);
		my $line2 = <$ZF>; chomp($line2);
		my $line3 = <$ZF>; chomp($line3);
		my $line4 = <$ZF>; chomp($line4);
		
		# read file format: @D00480L:375:H7TNCBCX2:1:1101:6690:2712 1:N:0:TCCGCGAA+AGGCGAAG
		# header file format: @D00480L:369:H555TBCX2:2:1101:9836:96824/1
		my @tmp1 = split(" ",$line1);
		my @tmp2 = split(":",$tmp1[1]);
		my $newh = "$tmp1[0]/$tmp2[0]";
		if(defined($heads{$newh}))
		{
			print OUT "$line1\n$line2\n$line3\n$line4\n";
		}
		else
		{
			#		
		}
	}
	close OUT;
}
