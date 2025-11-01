#!/usr/bin/env perl

# Original script by Maria Dalby <mdalbydk@gmail.com> and Robin Andersson <robin@bio.ku.dk> (https://github.com/anderssonlab/Replication_SCARseq/blob/78ffbabbfc8149cdc5bc5b5717878b26ebd4d8f1/README.md)
# Firstly modified by Nicolás Alcaraz <nicolas.alcaraz@cpr.ku.dk> (https://github.com/grothlab/SCARseq_Pipeline/blob/a4b327f1901ae6a980767d05ec7af79896a604c9/libs/partition_smooth.pl)
# Finally modified by Samuel Ruiz-Pérez <samper@cancer.dk> (2025)

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

## Keep a stack of RFD values, smoothed RFD values and window IDs for derivative calculation
my @RFD_stack = ();
my @RFD_smooth_stack = ();
my @ID_stack = ();
my $size = 0;
my $size_smooth = 0;

## initiate derivative weights as in Nonparametric Derivative Estimation (Brabanter et al.)
my @w = ();
for (my $i=0; $i<=$dradius; $i++) {
    push(@w, $i**2);
}
my $wsum = sum(@w);
for (my $i=1; $i<=$dradius; $i++) {
    $w[$i] /= $wsum;
}

# Subroutine to reset stacks and counters
sub reset_stacks {
    @RFD_stack = ();
    @RFD_smooth_stack = ();
    @ID_stack = ();
    $size = 0;
    $size_smooth = 0;
}

# Subroutine to print remaining elements for the current chromosome
sub print_remaining {
    for (my $i=$radius; $i<(2*$radius); $i++) {
        last if $i-$radius > $#ID_stack || $i > $#RFD_smooth_stack;
        print "$ID_stack[$i-$radius]\t$RFD_stack[$i-$radius]\t$RFD_smooth_stack[$i]\tNA\tNA\tNA\n";
    }
    for (my $i=$radius; $i<(2*$radius); $i++) {
        last if $i > $#ID_stack;
        print "$ID_stack[$i]\t$RFD_stack[$i]\tNA\tNA\tNA\tNA\n";
    }
}

my $prev_chr = undef;
my $line_count = 0;
my $chr_count = 0;
my $total_lines = 0;  # Counter for total lines/windows processed

while(my $fline = <$fin>) {
    my $rline = <$rin>;
    if (!defined $rline) {
        warn timestamp() . " Warning: $ARGV[2] ended before $ARGV[1].\n";
        last;
    }

    my @F = split(/\t/,$fline);
    my @R = split(/\t/,$rline);

    # Extract chromosome name (everything before last underscore)
    my ($curr_chr) = $F[0] =~ /^(.*)_/;
    if (!defined $curr_chr) {
        die timestamp() . " Could not extract chromosome from $F[0]\n";
    }

    # If chromosome changed, print remaining and reset stacks/counters
    if (!defined($prev_chr) || $curr_chr ne $prev_chr) {
        if (defined $prev_chr) {
            print STDERR timestamp() . " Finished chromosome $prev_chr (lines processed: $line_count)\n";
            print_remaining();
        }
        reset_stacks();
        $prev_chr = $curr_chr;
        $chr_count++;
        $line_count = 0;
        print STDERR timestamp() . " Starting chromosome $curr_chr (chromosome #$chr_count)\n";
    }

    $line_count++;
    $total_lines++;  # Increment total lines/windows processed
    print STDERR timestamp() . " Processed $line_count lines for $curr_chr\r" if ($line_count % 100000 == 0);

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

    ## Push new item values to the stacks
    push(@RFD_stack, $RFD);
    push(@ID_stack, $F[0]);
    $size++;

    if ($size < ($radius*2+1)) {
        if ($size < ($radius+1)) {
            print "$ID_stack[$size-1]\t$RFD_stack[$size-1]\tNA\tNA\tNA\tNA\n";
        }
        next;
    }

    ## calculate smoothed RFD (uniform blur)
    my $RFD_smooth = 0;
    for (my $i=0; $i<($radius*2+1); $i++) {
        $RFD_smooth += $RFD_stack[$i];
    }
    $RFD_smooth /= $radius*2+1;

    # Push new value to the stack
    push(@RFD_smooth_stack, $RFD_smooth);
    $size_smooth++;

    if ($size_smooth < ($radius*2+1)) {
        if ($size_smooth < ($radius+1)) {
            print "$ID_stack[$radius+$size_smooth-1]\t$RFD_stack[$radius+$size_smooth-1]\t$RFD_smooth_stack[$size_smooth-1]\tNA\tNA\tNA\n";
        }
    }
    else {
        ## calculate RFD derivative
        my $RFD_deriv = 0;
        for (my $i=1; $i<=$dradius; $i++) {
            $RFD_deriv += $w[$i] * ($RFD_smooth_stack[$radius+$i]-$RFD_smooth_stack[$radius-$i]) / (2*$i+1);
        }

        ## caclculate RFD boundary score
        my $score = $RFD_deriv*(1.0-abs($RFD_smooth_stack[$radius]));

        ## restrict derivative reporting at zero crossings (zradius)
        my $zero_crossing_derivative = "NA";
        my $min_RFD=1;
        my $max_RFD=-1;
        for (my $i=-1*$zradius; $i<=$zradius; $i++) {
            if ($RFD_smooth_stack[$radius+$i]>$max_RFD) {
                $max_RFD=$RFD_smooth_stack[$radius+$i];
            }
            if ($RFD_smooth_stack[$radius+$i]<$min_RFD) {
                $min_RFD=$RFD_smooth_stack[$radius+$i];
            }
        }
        if ($min_RFD<0 && $max_RFD>0) {
            $zero_crossing_derivative = $RFD_deriv;
        }

        ## print
        print "$ID_stack[0]\t$RFD_stack[0]\t$RFD_smooth_stack[$radius]\t$RFD_deriv\t$score\t$zero_crossing_derivative\n";

        shift @RFD_smooth_stack;
        $size_smooth--;
    }

    ## Remove first items from the stack
    shift @RFD_stack;
    shift @ID_stack;
    $size--;
}

# Print remaining elements for the last chromosome
print STDERR timestamp() . " Finished chromosome $prev_chr (lines processed: $line_count)\n" if defined $prev_chr;
print_remaining();

print STDERR timestamp() . " Processing complete. Chromosomes processed: $chr_count\n";
print STDERR timestamp() . " Total lines/windows processed: $total_lines\n";

close($fin);
close($rin);