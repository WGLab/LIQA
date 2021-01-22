#!/usr/bin/perlq  
 
###### 

use Getopt::Long;
use Pod::Usage;

my $gp1;
my $gp2;
my $output;

GetOptions('gp1=s'=>\$gp1,'gp2=s'=>\$gp2,'o=s'=>\$output);

if((!($gp1))||(!($output))||(!($gp2))){
    pod2usage();
}

my @file1;
my @file2; 

open FP, "$gp1";
while(<FP>) {
    chomp();
    @file1 = (@file1, $_);
}
close FP;

open FP, "$gp2";
while(<FP>) {
    chomp();
    @file2 = (@file2, $_);
}
close FP;



#my $allfile = (@file1, @file2);

my %alldata;
my %spid;

for my $i (0..$#file1) {
    my %isgenedone;
    open FP, "$file1[$i]";
    while(<FP>) {
	chomp();
	next if /GeneName/;
	my @a = split("\t");
	my $tmpgene = $a[0];
	my $tmpiso = $a[1];
	my $tmpest = $a[3];
	$isgenedone{$tmpgene} = 1;
    }
    close FP;
    foreach my $gene (keys %isgenedone) {
	$spid{$gene}++;
    }
    open FP, "$file1[$i]";
    while(<FP>) {
        chomp();
	next if /GeneName/;
        my @a = split("\t");
        my $tmpgene = $a[0];
        my $tmpiso = $a[1];
        my $tmpest = $a[3];
	$alldata{$tmpgene}{$tmpiso}{$spid{$tmpgene}} = "1;".$tmpest;
    }
    close FP;
}

for my $i (0..$#file2) {
    my %isgenedone;
    open FP, "$file2[$i]";
    while(<FP>) {
        chomp();
	next if /GeneName/;
        my @a = split("\t");
        my $tmpgene = $a[0];
        my $tmpiso = $a[1];
        my $tmpest = $a[3];
        $isgenedone{$tmpgene} = 1;
    }
    close FP;
    foreach my $gene (keys %isgenedone) {
        $spid{$gene}++;
    }
    open FP, "$file2[$i]";
    while(<FP>) {
        chomp();
	next if /GeneName/;
        my @a = split("\t");
        my $tmpgene = $a[0];
        my $tmpiso = $a[1];
        my $tmpest = $a[3];
        $alldata{$tmpgene}{$tmpiso}{$spid{$tmpgene}} = "2;".$tmpest;
    }
    close FP;
}

open OUT, ">$output";
foreach my $gene (keys %alldata) {
    foreach my $iso (keys %{$alldata{$gene}}) {
	for my $i (1..$spid{$gene}) {
	    my $out = $alldata{$gene}{$iso}{$i};
	    my @a = split(";", $out);
	    print OUT "$gene\t$iso\t$a[0]\t$a[1]\tinformation\t$i\n";
	}
    }
}
close OUT;

=head1 SYNOPSIS

    -gp1 ---estimation file for group 1
    
    -gp2 ---estimation file for group 2

    -o ---The file name that you want to save the results
