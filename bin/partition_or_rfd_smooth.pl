#!/usr/bin/env perl

# Original script by Maria Dalby <mdalbydk@gmail.com> and Robin Andersson <robin@bio.ku.dk> (https://github.com/anderssonlab/Replication_SCARseq/blob/78ffbabbfc8149cdc5bc5b5717878b26ebd4d8f1/README.md)
# Firstly modified by Nicolas Alcaraz <nicolas.alcaraz@cpr.ku.dk> (https://github.com/grothlab/SCARseq_Pipeline/blob/a4b327f1901ae6a980767d05ec7af79896a604c9/libs/partition_smooth.pl)
# Most recently modified by Samuel Ruiz-Pérez <samper@cancer.dk> (2025)

## Usage:
##   ./partition_or_RFD_smooth.pl <RFD|partition> <F.tab> <R.tab> <radius> <dradius> <zradius>
##
## Arguments:
##   <RFD|partition>   : Specify 'RFD' for OK-seq or 'partition' for SCAR-seq signal calculation.
##   <F.tab>           : Input file with forward strand counts (tab-delimited).
##   <R.tab>           : Input file with reverse strand counts (tab-delimited).
##   <radius>          : Smoothing window radius (integer, >= 1).
##   <dradius>         : Derivative window radius (integer, <= radius).
##   <zradius>         : Zero-crossing window radius (integer, <= radius).
##
## Input:
##   F.tab and R.tab list the same windows in the same order, one per line, with the window ID
##   (<chrom>_<window number>) in column 1 and the count in column 4. Each chromosome's lines
##   must be together, but can be in any order: windows are smoothed in order of their number.
##   One chromosome is held in memory at a time.
##
## Output:
##   Prints tab-delimited columns: ID, RFD/partition value, smoothed value, derivative, boundary score, zero-crossing derivative.

use strict;
use warnings;
use List::Util qw(sum);

sub timestamp {
    my @t = localtime();
    return sprintf("[%04d-%02d-%02d %02d:%02d:%02d]", $t[5]+1900, $t[4]+1, $t[3], $t[2], $t[1], $t[0]);
}

my $RFD_or_partition = $ARGV[0];
if (!defined $RFD_or_partition || ($RFD_or_partition ne "RFD" && $RFD_or_partition ne "partition")) {
    die timestamp() . " Invalid argument: $RFD_or_partition. Use 'RFD' or 'partition'.\n";
}

open(my $fin ,"$ARGV[1]") or die timestamp() . " Could not open $ARGV[1]: $!\n";
open(my $rin, "$ARGV[2]") or die timestamp() . " Could not open $ARGV[2]: $!\n";

my $radius = $ARGV[3];
my $dradius = $ARGV[4];
my $zradius = $ARGV[5];

if (!defined $radius || !defined $dradius || !defined $zradius) {
    die timestamp() . " Missing radius/dradius/zradius arguments.\n";
}

if ($dradius > $radius) {
    die timestamp() . " dradius ($dradius) cannot be greater than radius ($radius).\n";
}

if ($zradius > $radius) {
    die timestamp() . " zradius ($zradius) cannot be greater than radius ($radius).\n";
}

print STDERR timestamp() . " Starting processing with parameters: mode = $RFD_or_partition, radius = $radius, dradius = $dradius, zradius = $zradius\n";
print STDERR timestamp() . " Input files: F = $ARGV[1], R = $ARGV[2]\n";

## initiate derivative weights as in Nonparametric Derivative Estimation (Brabanter et al.)
my @w = ();
for (my $i=0; $i<=$dradius; $i++) {
    push(@w, $i**2);
}
my $wsum = sum(@w);
for (my $i=1; $i<=$dradius; $i++) {
    $w[$i] /= $wsum;
}

## Windows of the chromosome being read: ID, RFD and window number, in input order
my @ID = ();
my @RFD = ();
my @number = ();
my $current_chr = undef;
my %finished = ();
my $chr_count = 0;
my $total_lines = 0;

## Smooth one chromosome's windows in order of their number and print them in input order
sub process_chromosome {
    my ($chr) = @_;

    # Order the windows along the genome
    my @order = sort { $number[$a] <=> $number[$b] } (0 .. $#number);
    for (my $k=1; $k<=$#order; $k++) {
        if ($number[$order[$k]] == $number[$order[$k-1]]) {
            die timestamp() . " Window $chr\_$number[$order[$k]] appears more than once.\n";
        }
    }
    my @r = map { $RFD[$_] } @order;
    my $n = scalar(@r);

    # Results per input line; undefined values are printed as NA
    my @RFD_smooth = ();
    my @RFD_deriv = ();
    my @score = ();
    my @zero_crossing_derivative = ();

    ## calculate smoothed RFD (uniform blur), for windows with $radius windows on both sides
    my @s = ();
    for (my $j=$radius; $j<$n-$radius; $j++) {
        my $RFD_smooth = 0;
        for (my $i=$j-$radius; $i<=$j+$radius; $i++) {
            $RFD_smooth += $r[$i];
        }
        $RFD_smooth /= $radius*2+1;
        $s[$j] = $RFD_smooth;
        $RFD_smooth[$order[$j]] = $RFD_smooth;
    }

    ## derivative, boundary score and zero crossings, for windows with smoothed values $radius windows on both sides
    for (my $j=2*$radius; $j<$n-2*$radius; $j++) {
        ## calculate RFD derivative
        my $RFD_deriv = 0;
        for (my $i=1; $i<=$dradius; $i++) {
            $RFD_deriv += $w[$i] * ($s[$j+$i]-$s[$j-$i]) / (2*$i+1);
        }

        ## calculate RFD boundary score
        my $boundary_score = $RFD_deriv*(1.0-abs($s[$j]));

        ## restrict derivative reporting at zero crossings (zradius)
        my $min_RFD=1;
        my $max_RFD=-1;
        for (my $i=-1*$zradius; $i<=$zradius; $i++) {
            if ($s[$j+$i]>$max_RFD) {
                $max_RFD=$s[$j+$i];
            }
            if ($s[$j+$i]<$min_RFD) {
                $min_RFD=$s[$j+$i];
            }
        }

        $RFD_deriv[$order[$j]] = $RFD_deriv;
        $score[$order[$j]] = $boundary_score;
        $zero_crossing_derivative[$order[$j]] = $RFD_deriv if ($min_RFD<0 && $max_RFD>0);
    }

    ## print, in the input order
    for (my $l=0; $l<=$#ID; $l++) {
        print join("\t", $ID[$l], $RFD[$l], map { defined $_ ? $_ : "NA" } ($RFD_smooth[$l], $RFD_deriv[$l], $score[$l], $zero_crossing_derivative[$l])) . "\n";
    }

    $chr_count++;
    $total_lines += $n;
    print STDERR timestamp() . " Finished chromosome $chr ($n windows)\n";
}

## Read both tables line by line, one chromosome at a time
while (my $fline = <$fin>) {
    my $rline = <$rin>;
    die timestamp() . " $ARGV[2] has fewer lines than $ARGV[1].\n" if !defined $rline;

    my @F = split(/\t/,$fline);
    my @R = split(/\t/,$rline);

    if ($F[0] ne $R[0]) {
        die timestamp() . " Line $.: window $F[0] in $ARGV[1] but $R[0] in $ARGV[2]. Both must list the same windows in the same order.\n";
    }

    # Window IDs are <chrom>_<window number>; the chromosome name can itself contain underscores
    my ($chr, $number) = $F[0] =~ /^(.*)_(\d+)$/;
    if (!defined $chr) {
        die timestamp() . " Could not extract chromosome and window number from $F[0]\n";
    }

    ## calculate RFD
    my $RFD;
    if ($RFD_or_partition eq "partition") {
        # For partition signal (SCAR-seq), we use the formula:
        # Partition = (F - R)/(F + R)
        $RFD = ($F[3] - $R[3])/($F[3] + $R[3]);
    } else {
        # For RFD signal (OK-seq), we use the formula:
        # RFD = (R - F)/(R + F)
        $RFD = ($R[3] - $F[3])/($F[3] + $R[3]);
    }

    # A new chromosome: finish the previous one
    if (!defined($current_chr) || $chr ne $current_chr) {
        if (defined $current_chr) {
            process_chromosome($current_chr);
            $finished{$current_chr} = 1;
        }
        if ($finished{$chr}) {
            die timestamp() . " Line $.: the windows of chromosome $chr are not all together in the input.\n";
        }
        @ID = ();
        @RFD = ();
        @number = ();
        $current_chr = $chr;
    }

    push(@ID, $F[0]);
    push(@RFD, $RFD);
    push(@number, $number);
}
die timestamp() . " $ARGV[2] has more lines than $ARGV[1].\n" if defined <$rin>;
process_chromosome($current_chr) if defined $current_chr;
close($fin);
close($rin);

print STDERR timestamp() . " Processing complete. Chromosomes processed: $chr_count\n";
print STDERR timestamp() . " Total lines/windows processed: $total_lines\n";
