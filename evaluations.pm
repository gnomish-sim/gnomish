# Evaluations module for Gnomish
#
# File: Evaluations.pm
# Author: Hai Pham, hpham42@gmail.com
# Last Changed: Mar 6, 2003
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


package Gnomish::Evaluations;

use strict;
use Carp;

require Exporter;
our @ISA = qw (Exporter);
our @EXPORT = qw (computeTraitPhenotypeMean computeTraitPhenotypeVariance
		 selectGenotypeArray 
		 countAlleleStatesInGenotype getAlleleStateCountInGenotype
		 countAlleleStatesInPopulation getAlleleStateCountInPopulation);

		 
use Gnomish::Stdlib;

# Compute the phenotypic mean value for a particular trait.
# Make sure that the phenotypic values for the particular TraitName has been
# computed first!
sub computeTraitPhenotypeMean
{
	my %params = @_;

	croak "'Population' missing" unless exists $params{Population};
	croak "'TraitName' missing" unless exists $params{TraitName};

	my $pop = $params{Population};
	if (ref $pop eq "Gnomish::Population")
	{
		croak "Population object must be passed by reference";
	}
	my $trait_name = $params{TraitName};
	
	my $nmembers = $$pop->nMembers();
	my ($n, $m);
	my $sum = 0;
	for ($n = 0; $n < $nmembers; $n++)
	{
		$m = $$pop->getMemberRef (Member => $n);
		$sum += $$m->getPhenotypicValueForTrait (TraitName => $trait_name); 
	}
	
	return ($sum / $nmembers);
}

# compute the phenotypic variance for the named trait.
sub computeTraitPhenotypeVariance
{
	my %params = @_;

	croak "'Population' missing" unless exists $params{Population};
	croak "'TraitName' missing" unless exists $params{TraitName};

	my $pop = $params{Population};
	if (ref $pop eq "Gnomish::Population")
	{
		croak "Population object must be passed by reference";
	}
	my $trait_name = $params{TraitName};

	my $nmembers = $$pop->nMembers();
	my ($n, $m, $pval);
	my $sum_x = 0;
	my $sum_x2 = 0;
	for ($n = 0; $n < $nmembers; $n++)
	{
		$m = $$pop->getMemberRef (Member => $n);
		$pval = $$m->getPhenotypicValueForTrait (TraitName => $trait_name);
		
		$sum_x += $pval;
		$sum_x2 += ($pval * $pval);
	}
	
	my $v = ($sum_x2 - (($sum_x * $sum_x) / $nmembers)) / ($nmembers - 1.0);
	
	return $v;
}


# Select genotypes from a population with a particular locus state at specified locus
# Return an array with the usual suspects.
sub selectGenotypeArray
{
	my %params = @_;
	
	croak "'Population' missing" unless exists $params{Population};
	croak "'Locus' missing" unless exists $params{Locus};
	croak "'LocusState' missing" unless exists $params{LocusState};
	
	my $pop = $params{Population};
	if (ref $pop eq "Gnomish::Population")
	{
		croak "Population must be passed by reference";	
	}
	
	my $locus = $params{Locus};
	my $locus_state = $params{LocusState};
	
	my $members = $$pop->nMembers();
	my ($i, $m, @mlist, $state);
	for ($i = 0; $i < $members; $i++)
	{
		$m = $$pop->getMemberRef (Member => $i);
		$state = $$m->getLocusState (Locus => $locus);
		if ($state eq $locus_state)
		{
			push @mlist, $m;	
		}
	}
	
	return @mlist;
}

# compute the number of different allele states in the allele strings of the genotype.
# return the result as a hash object.
# $genotype should be a reference to the genotype! 
sub countAlleleStatesInGenotype
{
	my %params = @_;
	
	croak "'Genotype' missing" unless exists $params{Genotype};
	croak "'OrdinateAlleleStates' missing" unless exists $params{OrdinateAlleleStates};
	my $genotype = $params{Genotype};
	
	if ( ref $genotype eq "Gnomish::Genotype" )
	{
		croak "Gnotype object must be passed by reference";
	}
	
	my $ordinate_states = $params{OrdinateAlleleStates};
	if ( not ref $ordinate_states )
	{
		croak "The OrdinateAlleleStates hash must be passed by reference";
	}
	
	# get the two allele strings from the genotype
	my @a_string0 = $$genotype->getGameteAsArray ( Gamete => 0 );
	my @a_string1 = $$genotype->getGameteAsArray ( Gamete => 1 );
	
	my $loci = scalar @a_string0;
	my (%allele_states, $allele_state, $ordinate_state);
	for (my $i = 0; $i < $loci; $i++)
	{
		$allele_state = $a_string0[$i] . $a_string1[$i];
		$ordinate_state = $$ordinate_states{$allele_state};
		
		$allele_states{$ordinate_state}++;
	}
	
	return %allele_states;
}

# return the number of times a particular allele state appears in a genotype
sub getAlleleStateCountInGenotype
{
	my %params = @_;
	
	croak "'Genotype' missing" unless exists $params{Genotype};
	croak "'AlleleState' missing" unless exists $params{AlleleState};
	
	my $genotype = $params{Genotype};
	my $allele_state = $params{AlleleState};
	
	if ( ref $genotype eq "Gnomish::Genotype" )
	{
		croak "Genotype object must be passed by reference";
	}
	
	# compute the ordinate allele state set
	# first, get the allele set
	my $genome = $$genotype->getGenomeRef();
	my @allele_set = split ',', $$genome->getAlleleSet();
	my %ordinate_states = computeAlleleStateSet (AlleleSet => \@allele_set);
	
	my %allele_states = countAlleleStatesInGenotype ( Genotype => $genotype,
				OrdinateAllelestates => \%ordinate_states );
	
	return $allele_states{$allele_state};
}

# count the occurances of all allele states in the population
sub countAlleleStatesInPopulation
{
	my %params = @_;

	croak "'Population' missing" unless exists $params{Population};
	
	my $population = $params{Population};
	if ( ref $population eq "Gnomish::Population" )
	{
		croak "Population object must be passed by reference";
		
	}
	
	# compute the ordinate allele state set
	# first, get the allele set
	my $g = $$population->getMemberRef (Member => 0);
	my $genome = $$g->getGenomeRef();
	my @allele_set = split ',', $$genome->getAlleleSet();
	my %ordinate_states = computeAlleleStateSet (AlleleSet => \@allele_set);
	
	my %allele_states;
	for (my $i = 0; $i < $$population->nMembers(); $i++)
	{
		my $genotype = $$population->getMemberRef ( Member => $i );
		
		my %g_allele_states = countAlleleStatesInGenotype ( Genotype => $genotype,
					OrdinateAlleleStates => \%ordinate_states );
		
		foreach my $state (keys %g_allele_states)
		{
			$allele_states{$state} += $g_allele_states{$state};
		}
	}
	
	return %allele_states;
}

# get the count of the number of times a particular allele state occurs in the genotypes in the population
sub getAlleleStateCountInPopulation
{
	my %params = @_;
	
	croak "'Population' missing" unless exists $params{Population};
	croak "'AlleleState' missing" unless exists $params{AlleleState};
	
	my $population = $params{Population};
	if ( ref $population eq "Gnomish::Population" )
	{
		croak "Population object must be passed by reference";
	}
	
	my $allele_state = $params{AlleleState};
	
	my %allele_states = countAlleleStatesInPopulation ( Population => $population );
	
	return $allele_states{$allele_state};
}

1;
