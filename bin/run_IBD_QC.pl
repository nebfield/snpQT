#!/usr/bin/perl

use strict;

# Store first argument (plink.imiss to variable $imiss and pihat_min0.125.genome to $pihat_file, from the command line
my ($imiss_file, $pihat_file)=  @ARGV;

# Checking if the arguments are passed
if (not defined $imiss_file) {
  die ".imiss file is not defined\n";
}

if (not defined $pihat_file) {
  die ".genome file is not defined\n";
}

my %imiss;
my %removed;


open IMISS, '<', $imiss_file
        or die "Cannot open genotypes file (pihat_0.125.genome): $!\n";
print "Reading PLINK .imiss file plink.imiss\n";
while(<IMISS>){
	s/^\s+//;
    my @fields = split /\s+/, $_;
    $imiss{$fields[0]}{$fields[1]} = $fields[5];
}

open GENOME, '<',  $pihat_file
        or die "Cannot open genotypes file (pihat_0.125.genome): $!\n";
open OUT, '>', "pihat_failed_samples.txt";
print "Reading PLINK .genome file pihat_min0.125.genome\n";

print "For pair of samples which are unrelated check which sample has the highest missingness and print it to pihat_failed_samples.txt...\n";
while(<GENOME>){
    s/^\s+//;
    my @fields = split /\s+/, $_;
	if($fields[9] > 0.125){
 		if($imiss{$fields[0]}{$fields[1]}>$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[0]}{$fields[1]}){
 				print OUT "$fields[0] $fields[1]\n";
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 		elsif($imiss{$fields[0]}{$fields[1]}<$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[2]}{$fields[3]}){
 				print OUT "$fields[2] $fields[3]\n";
 				$removed{$fields[2]}{$fields[3]} = 1;
 			}
 		}
 		else{
 			unless($removed{$fields[0]}{$fields[1]}){
 				print OUT "$fields[0] $fields[1]\n";
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 	}
}
print "Done!\n";
