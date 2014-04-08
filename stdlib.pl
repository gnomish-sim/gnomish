# Gnomish std functions module
#
# File: Stdlib.pm
# Author: Hai Pham, hpham42@gmail.com
# Last Modified: Jan 21, 2003
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */
#

package Gnomish::Stdlib;

use strict;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(initialize computeAlleleStateSet gaussianRand 
		explode createRecombinationDistanceString);

# put stuff in here that need to be initialized at the start of the simulation
sub initialize
{
#	print "initializing random number generator\n";
	
	# initialize the random number generator;
	srand( time() ^ ($$ + ($$ << 15)) );
	
}

# the 'AlleleSet' is an array that contains the valid allele values in the
# genome.  For example, if an allele can have a value of 0, 1, 2, or 3
# then AlleleSet is (0, 1, 2, 3).  This function computes the
# various pairwise permutations that are possible and their ordinate value.
# for example, the set 0, 1, 2, 3 can have pairwise relationships like:
# 00, 01, 02, 03, 10, 11, 12, 13, 20, 21, 22, 23, .....
# however, 10 is the same as 01, 20 is the same as 02, etc.
# This function returns a look up table that maps the ordinate value
# of a pair. 
sub computeAlleleStateSet
{
	my %params = @_;

	croak "'AlleleSet' missing" unless exists $params{AlleleSet};

	my $s = $params{AlleleSet};
	if (not ref $s)
	{
		croak "'AlleleSet' must be passed by reference";	
	}
	
	my @set = @$s;
	my $alleles = scalar @set;
	
	my ($i, $j, %allele_states, $a, $b);
	for ($i = 0; $i < $alleles; $i++)
	{
		for ($j = 0; $j < $alleles; $j++)
		{
			$a = $set[$i];
			$b = $set[$j];
			
			$allele_states{"$a$b"} = "$b$a";
			$allele_states{"$b$a"} = "$b$a";
		}
	}
	
	return %allele_states;
}


# compute a random number according to a guasian distribution
sub gaussianRand
{
	my %params = @_;
	
	croak "'Variance' missing" unless exists $params{Variance};
	my $variance = $params{Variance};
	
	if ($variance == 0)
	{
		return 0;	
	}
	
	my $sd = sqrt ($variance);
	
	my ($u1, $u2, $w);
	do
	{
		$u1 = 2.0 * rand() - 1.0;
		$u2 = 2.0 * rand() - 1.0;
		$w = ($u1 * $u1) + ($u2 * $u2);
	} until ( $w < 1.0 ); # $w should never get below 0.0, there is no need to check $w > 0
	
	my $v = sqrt ( (-2.0 * log ($w)) / $w );
	
	return ($u1 * $v * $sd);
}


# explode a string into an array, with one character of the string
# per array element.
sub explode ($) 
{
	my $string = shift;

	my $len = length $string;
	
	my ($i, @str);
	for ($i = 0; $i < $len; $i++)
	{
		push @str, substr $string, $i, 1;
	}
 	
	return @str;	
}


# create a new recombination distnace string with the provided parameters
# Chromosomes: number of chromosomes in the genome
# ChromosomeLengths: comma separated list of the number loci in each chromosome
# RecombinationDistances: comma separated list of the recombination distance at each chromosome
sub createRecombinationDistanceString
{
	my %params = @_;
	
	croak "'Chromosomes' missing"
		unless exists $params{Chromosomes};
	croak "'ChromosomeLengths' missing"
		unless exists $params{ChromosomeLengths};
	croak "'RecombinationDistances' missing"
		unless exists $params{RecombinationDistances};
		
	my $chromosomes = $params{Chromosomes};
	
	my $cl = $params{ChromosomeLengths};
	$cl =~ s/[\s\n\r\t]//g;
	my @c_lengths = split ',', $cl;
	
	my $rd = $params{RecombinationDistances};
	$rd =~ s/[\s\n\r\t]//g;
	my @rdists = split ',', $rd;
	
	if ( (scalar @c_lengths) != $chromosomes )
	{
		croak "Number of chromosome lenghts not equal to number of chromosomes";
	}
	
	if ( (scalar @rdists) != $chromosomes )
	{
		croak "Number of recombination distances not equal to number of chromosomes";
	}
	
	my @r_dist;
	for (my $c = 0; $c < $chromosomes; $c++)
	{
		for (my $l = 0; $l < ($c_lengths[$c] - 1); $l++)
		{
			push @r_dist, "$rdists[$c]";
		}
		push @r_dist, "0.5";
	}
	
	# remove the last '0.5'
	pop @r_dist;
	
	my $r = join ',', @r_dist;
	return $r;
}

1;
