# Gnomish::Trait class definition
#
# File: Trait.pm
# Author: Hai Pham, hpham42@gmail.com
# Last Changed: Jan 23, 2003
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

package Gnomish::Trait;

use strict;
use Carp;

use Gnomish::Stdlib;

sub new
{
	my $class = shift;
	my %params = @_;

	my $self = {};
	bless $self, $class;

	$self->{Name} = "";
	if (exists $params{Name})
	{
		$self->{Name} = $params{Name};	
	}
	
	return $self;
}

# initialize the trait with the provided properties
# inputs:
#  - TraitMean:Real (for evaluating, this value is added to the trait phenotype)
#  - TraitVariance:Real (for evaluating, we add a random number from a Gausian distribution
#                        having this variance to the phenotype.)
#  - TraitLoci:Int[..] list of loci in the genome that affect the trait in an additive/dominant manner
#  - TraitAlleleEffects:Int[..][..] this array contains the genetic value for each allele state
sub makeTrait
{
	my $self = shift;
	my %params = @_;

	croak "'Name' missing" unless exists $params{Name};
	croak "'TraitMean' missing" unless exists $params{TraitMean};
	croak "'TraitVariance' missing" unless exists $params{TraitVariance};
	croak "'TraitAlleleEffects' missing" unless exists $params{TraitAlleleEffects};
	croak "'AlleleSet' missing" unless exists $params{AlleleSet}; 
	
	$self->{Name} = $params{Name};
	
	$self->{_TraitMean} = $params{TraitMean};
	$self->{_TraitVariance} = $params{TraitVariance};

	my $tae = $params{TraitAlleleEffects};
	if (not ref $tae)
	{
		croak "'TraitAlleleEffects' must be passed by reference";	
	}
	$self->setTraitAlleleEffects (TraitAlleleEffects => $tae);
	
	$self->_setAlleleSet (AlleleSet => $params{AlleleSet});
	
}



# the expected reference should be an array of hashes
# the row index corresponds to the loci
# ("2" => \( "00" => 2, "01" => 1...),  
# "4" => \( "00" => 4, "01" => 5...) )
sub setTraitAlleleEffects
{
	my $self = shift;
	my %params = @_;

	croak "'TraitAlleleEffects' missing" unless exists $params{TraitAlleleEffects};
	
	my $tae = $params{TraitAlleleEffects};
	if (not ref $tae)
	{
		croak "'TraitAlleleEffects' must be passed by reference";	
	}
	
	my (%effects_table, $locus);
	foreach $locus (keys %$tae)
	{
		my %state_effects = %{ $$tae{$locus} };
		$effects_table{$locus} = \%state_effects;
	}
	
	$self->{_TraitAlleleEffects} = \%effects_table;
}


# set the lookup table for various allele
# the Allele set must be specified as a comma seperated string
sub _setAlleleSet
{
	my $self = shift;
	my %params = @_;

	croak "'AlleleSet' missing" unless exists $params{AlleleSet};
	my $set = $params{AlleleSet};
	
	# make sure that there are no spaces in the allele set string
	$set =~ s/[\s\n\r\t]//g;
	
	my @allele_set = split ',', $set;
	
	$self->{_AlleleSet} = \@allele_set;
	
	my %allele_states = computeAlleleStateSet (AlleleSet => \@allele_set);
	
	$self->{_AlleleStateSet} = \%allele_states;
}

# return the ordinate state value for the particular state.
sub _getAlleleStateOrd
{
	my $self = shift;
	my $astate = shift;

	return $self->{_AlleleStateSet}->{$astate};	
}

sub printTraitAlleleEffects
{
	my $self = shift;

	my $effects_table = $self->{_TraitAlleleEffects};

	my ($locus, $allele);
	foreach $locus (sort keys %$effects_table)
	{
		print "locus: $locus ";
		my %state_effects = %{ $$effects_table{$locus} };
		foreach $allele (sort keys %state_effects)
		{
			print "$allele:", $state_effects{$allele}, " ";	
		}
		print "\n";
	}
}

# return the array of possible effect values at the locus
sub getLocusEffects
{
	my $self = shift;
	my %params = @_;
	
	croak "'Locus' missing" unless exists $params{Locus};
	my $locus = $params{Locus};
	
	my $effects_table = $self->{_TraitAlleleEffects};
	my %state_effects = %{ $$effects_table{$locus} };
	
	return %state_effects;
}

# get the effect value of a particular allele state
sub getLocusEffectValue
{
	my $self = shift;
	my %params = @_;
	
	croak "'Locus' missing" unless exists $params{Locus};
	croak "'AlleleState' missing" unless exists $params{AlleleState};

	my $locus = $params{Locus};
	my $allele_state = $params{AlleleState};
	
	my $effects_table = $self->{_TraitAlleleEffects};
	my $state_ord = $self->_getAlleleStateOrd ($allele_state);

	return $effects_table->{$locus}{$state_ord};
}

# return an array loci that affect the trait
sub getActiveLoci
{
	my $self = shift;
	
	my $effects_table = $self->{_TraitAlleleEffects};
	my @loci = keys %$effects_table;
	
	return @loci;
}



sub setName
{
	my $self = shift;
	my %params = @_;
	
	croak "'Name' missing" unless exists $params{Name};
	$self->{Name} = $params{Name};
}


sub getName
{
	my $self = shift;

	return $self->{Name};	
}

sub getTraitVariance
{
	my $self = shift;

	return $self->{_TraitVariance};	
}

sub getTraitMean
{
	my $self = shift;

	return $self->{_TraitMean};	
}


1;

