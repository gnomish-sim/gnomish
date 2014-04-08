# Gnomish class Population
# This module incoporates both the class definition for the Population class
# and regular functions for operation on population objects.  Note that these
# functions are *EXTERNAL* to the class definition!
#
# File: Population
# Author: Hai Pham, hpham42@gmail.com
# Last Modified: June 10, 2003
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
 
package Gnomish::Population;


# this for the export of population operation functions that are not part of the Population
# class definition
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw (randomMate chainCrosses fullDiallel NC_DesignII backCrosses
		doubledHaploids mergePopulations
		selfedProgeny clonalProgeny singleSeedDescent
		specificHybrids
		makeRandomAllelePopulation makeRandomHomozygousPopulation
		makeHomozygousPopulation );

use strict;
use Carp;

use Gnomish::Genotype;

############## Functions that operate on the Population Class ################

# Make a new population where one or more doubled haploid progeny is created 
# from each member of a population.
# Input:
# - ParentPop:Population*
# - NProgeny:Int  (number of progenies *per parent*)
sub doubledHaploids
{
	my %params = @_;
	
	croak "'NProgeny' missing" unless exists $params{NProgeny};
	croak "'ParentPop' missing" unless exists $params{ParentPop};
	
	my $n_progeny    = $params{NProgeny};
	my $parent_pop   = $params{ParentPop};
	
	if (ref $parent_pop eq "Gnomish::Population")
	{
		croak "'ParentPop' must be passed by reference";
	}
	
	my $parents = $$parent_pop->{NMembers};
	croak "Parent population needs to have at least one member"  
		unless $parents > 0;
	
	my $new_pop = new Gnomish::Population;

	my ($i, $j, $parent);
	for ($i = 0; $i < $parents; $i++)
	{
		$parent = $$parent_pop->getMemberRef (Member => $i);
		for ($j = 0; $j < $n_progeny; $j++)
		{
			my $progeny = new Gnomish::Genotype;
			$progeny->makeWithSelfCrossedParent (Parent => $parent);
			$new_pop-> addMember (Genotype => $progeny);
		}
	}
	
	return $new_pop;
}


# Create a new population in which each parent in the parent population is crossed
# with itself one or times.
# Inputs:
# - ParentPop:Population*
# - NProgeny:Int  (number of progenies per parent)
# - NGenerations: Int;
sub selfedProgeny
{
	my %params = @_;

	croak "'ParentPop' missing" unless exists $params{ParentPop};
	croak "'NProgeny' missing" unless exists $params{NProgeny};
	croak "'NGenerations' missing" unless exists $params{NGenerations};
	
	my $parent_pop = $params{ParentPop};
	my $n_progeny  = $params{NProgeny};
	my $n_generations = $params{NGenerations};
	
	
	if (ref $parent_pop eq "Gnomish::Population")
	{
		croak "ParentPop must be passed by reference";	
	}

	for (my $n = 0; $n < $n_generations; $n++)
	{
		my $new_pop = new Gnomish::Population;
	
		my ($i, $j, $parent);	
		for ($i = 0; $i < $$parent_pop->nMembers; $i++)
		{
			$parent = $$parent_pop->getMemberRef (Member => $i);
			for ($j = 0; $j < $n_progeny; $j++)
			{
				my $progeny = new Gnomish::Genotype;
				$progeny->makeWithParents (Parent1 => $parent, Parent2 => $parent);
			
				$new_pop->addMember (Genotype => $progeny);
			}
		}
		
		$parent_pop = \$new_pop;
		
	}
		
	
	return $$parent_pop;
}

# Create a new population in which each parent is cloned one or more times.
# Inputs:
#  - ParentPop:Population*
#  - NProgeny:Int
sub clonalProgeny
{
	my %params = @_;

	croak "'ParentPop' missing" unless exists $params{ParentPop};
	croak "'NProgeny' missing" unless exists $params{NProgeny};

	my $parent_pop = $params{ParentPop};
	my $n_progeny  = $params{NProgeny};
	
	if (ref $parent_pop eq "Gnomish::Population")
	{
		croak "ParentPop must be passed by reference";	
	}
	
	my $new_pop = new Gnomish::Population;
	
	my ($i, $j, $parent);
	for ($i = 0; $i < $$parent_pop->{NMembers}; $i++)
	{
		# get a copy of the parent.
		# we use a copy of the parent, not use the reference!
		$parent = $$parent_pop->getMember (Member => $i);
		for ($j = 0; $j < $n_progeny; $j++)
		{
			# make a copy of the parent
			my $progeny = $parent;
			$new_pop->addMember (Genotype => $progeny);
		}
	}
	
	return $new_pop;
}

# Perform N generations of single seed descent.  This function creates one
# selfed progengy for each member of the parent population, and then
# each member in the new population is replaced by one of its selfed progeny.
# The process is repeated until N generations have been reached.
sub singleSeedDescent
{
	my %params = @_;

	croak "'ParentPop' missing" unless exists $params{ParentPop};
	croak "'NGenerations' missing" unless exists $params{NGenerations};
	
	my $parent_pop = $params{ParentPop};
	if (ref $parent_pop eq "Gnomish::Population")
	{
		croak "'ParentPop' must be passed by reference";	
	}

	my $n_generations = $params{NGenerations};
	
	# create the F1 generation
	my $population = new Gnomish::Population;
	my ($i, $parent);
	for ($i = 0; $i < $$parent_pop->{NMembers}; $i++)
	{
		$parent = $$parent_pop->getMemberRef (Member => $i);		
		
		my $progeny = new Gnomish::Genotype;
		$progeny->makeWithParents (Parent1 => $parent, Parent2 => $parent);
		
		$population->addMember (Genotype => $progeny);
	}
	
	my ($n, $parents, $kids);
	$parents = $population;
	for ($n = 0; $n < $n_generations; $n++)
	{
		 $kids = new Gnomish::Population;
		 
		 for ($i = 0; $i < $parents->{NMembers}; $i++)
		 {
			$parent = $parents->getMemberRef (Member => $i);
			
			my $progeny = new Gnomish::Genotype;
			$progeny->makeWithParents (Parent1 => $parent, Parent2 => $parent);
			
			$kids->addMember (Genotype => $progeny);
		 }
		 
		 # now the kids become the parents of the next generation
		$parents = $kids;				
	}	
	
	return $parents;
}



# Randomly mate genotypes in the population, creating new population
# Take existing population, make a new population by crossing random parents
# from the existing population (ie. mate two random parents NProgeny times).  Each cross
# results in one new member of the the new population.
# Input: 
# - ParentPop:Population*
# - NProgeny:Int
sub randomMate
{
	my %params = @_;
	
	croak "'NProgeny' missing" unless exists $params{NProgeny};
	croak "'ParentPop' missing" unless exists $params{ParentPop};
	croak "'NGenerations' missing" unless exists $params{NGenerations};
		
	if ( (ref $params{ParentPop}) eq "Gnomish::Population" )
	{
		croak "ParentPop must be passed by reference";	
	}
	
	my $n_progeny    = $params{NProgeny};
	my $parent_pop   = $params{ParentPop};
	my $n_generations = $params{NGenerations};

	if ($$parent_pop->{NMembers} < 1)
	{
		croak "Parent population must have at least one member for mating to work";	
	}
	
	for (my $n = 0; $n < $n_generations; $n++)
	{
	
		my $new_pop = new Gnomish::Population;
	
		my $parents = $$parent_pop->{NMembers};
	
		my ($p1, $p2, $n);
		my ($mom, $dad);
		for ($n = 0; $n < $n_progeny; $n++)
		{
			# pick 2 random parents
			$p1 = int rand($parents);
			$mom = $$parent_pop->getMemberRef (Member => ($p1 - 1));
			$p2 = int rand($parents);
			$dad = $$parent_pop->getMemberRef (Member => ($p2 - 1));
		
			my $progeny = new Gnomish::Genotype;
				
				$progeny->makeWithParents( Parent1 => $mom,
					Parent2 => $dad);
					$new_pop->addMember (Genotype => $progeny);
		}
		
		$parent_pop = \$new_pop;
	}
	
	return $$parent_pop;
}


# Create a new population by sequentially crossing members of the population.
# This procedure crosses the first individual with the second, the second with the third, etc.
# The last individual is crossed with the first so that the number of crosses made is equal
# to the number of individuals in the parental population.  One progeny is made from each cross.
#
# Inputs:
#  - ParentPop:Population*
sub chainCrosses
{
	my %params = @_;
	
	croak "'ParentPop' missing" unless exists $params{ParentPop};
	
	if ( (ref $params{ParentPop}) eq "Gnomish::Population" )
	{
		croak "ParentPop object must be passed by reference";	
	}
 
	my $parent_pop = $params{ParentPop};
	
	my $parents = $$parent_pop->{NMembers};
	if ($parents < 2)
	{
		croak "Parent population must have at least two members";	
	}
	
	my $new_pop = new Gnomish::Population;
	
	my ($mom, $dad, $i);
	for ($i = 0; $i < ($parents - 1); $i++)
	{
		$mom = $$parent_pop->getMemberRef (Member => $i);
		$dad = $$parent_pop->getMemberRef (Member => $i + 1);
		
		my $progeny = new Gnomish::Genotype;
		$progeny->makeWithParents (Parent1 => $mom, Parent2 => $dad);
		$new_pop->addMember (Genotype => $progeny);
	}
	
	# now do the crossing of the last member of the pop with the first
	$mom = $$parent_pop->getMemberRef (Member => ($parents - 1));
	$dad = $$parent_pop->getMemberRef (Member => 0);
	my $progeny = new Gnomish::Genotype;
	$progeny->makeWithParents (Parent1 => $mom, Parent2 => $dad);
	$new_pop->addMember (Genotype => $progeny);

	return $new_pop;	
}

# Create a new population in which every individual is crossed with every others in the
# population, including itself.
# Input:
# - ParentPop:Population*
sub fullDiallel
{
	my %params = @_;
	
	croak "'ParentPop' parameter missing" unless exists $params{ParentPop};
	
	if (ref $params{ParentPop} eq "Gnomish::Population")
	{
		croak "'ParentPop' object must be passed by reference";
	}

	my $parent_pop = $params{ParentPop};
	my $parents = $$parent_pop->{NMembers};
	croak "Parent population must have at least 1 member" unless $parents > 0;
	
	return NC_DesignII (ParentPop1 => $params{ParentPop}, 
				ParentPop2 => $params{ParentPop});
}


# This mates all individuals in one population with all members of another population, making
# one progeny from each cross.  
# Inputs:
# - ParentPop1:Population*
# - ParentPop2:Population* 
sub NC_DesignII
{
	my %params = @_;
	
	croak "'ParentPop1' missing" unless exists $params{ParentPop1};
	croak "'ParentPop2' missing" unless exists $params{ParentPop2};
	
	if (ref $params{ParentPop1} eq "Gnomish::Population" or
		ref $params{ParentPop2} eq "Gnomish::Population")
	{
		croak "ParentPop objects must be passed by reference";	
	}
	
	my $parent_pop1 = $params{ParentPop1};
	my $parent_pop2 = $params{ParentPop2};
	
	my $parents1 = $$parent_pop1->{NMembers};
	croak "ParentPop1 must have at least one member" unless $parents1 > 0;
	my $parents2 = $$parent_pop2->{NMembers};
	croak "ParentPop2 must have at least one member" unless $parents2 > 0;
	
	my $new_pop = new Gnomish::Population;
	
	my ($m, $d, $mom, $dad);
	for ($m = 0; $m < $parents1; $m++)
	{
		$mom = $$parent_pop1->getMemberRef (Member => $m);
		
		for ($d = 0; $d < $parents2; $d++)
		{
			$dad = $$parent_pop2->getMemberRef (Member => $d);
			
			my $progeny = new Gnomish::Genotype;
			$progeny->makeWithParents (Parent1 => $mom, Parent2 => $dad);
			$new_pop->addMember (Genotype => $progeny);
		}
	}
	
	return $new_pop;
}

# This performs crosses between all individuals in one population (the donor population) and
# a single individual (the recurent parent) in a second population.  One or more progenies
# can be produced from each cross.
# Inputs:
#  - RecurentParent:Genotype*
#  - DonorPop:Population*
sub backCrosses
{
	my %params = @_;
	
	croak "'RecurrentParent' missing" unless exists $params{RecurrentParent};
	croak "'DonorPop' missing" unless exists $params{DonorPop};
	
	if (ref $params{RecurrentParent} eq "Gnomish::Genotype")
	{
		croak "'RecurrentParent' must be passed by reference";
	}
	
	if (ref $params{DonorPop} eq "Gnomish::Population")
	{
		croak "'DonorPop' must be passed by reference";	
	}

	my $dad = $params{RecurrentParent};
	my $donor_pop = $params{DonorPop};
	my $moms = $$donor_pop->{NMembers};
	
	my $new_pop = new Gnomish::Population;
	
	my ($m, $mom);
	for ($m = 0; $m < $moms; $m++)
	{
		$mom = $$donor_pop->getMemberRef (Member => $m);
		
		my $progeny = new Gnomish::Genotype;
		$progeny->makeWithParents (Parent1 => $dad, Parent2 => $mom);
		
		$new_pop->addMember ( Genotype => $progeny);
	}
	
	return $new_pop;
}

# This allows the user to make specific crosses.  The user must specify the parents of each
# cross manually. 
sub specificHybrids
{
	die "specificHybrids crossing has not been implemented yet!";
}


# Merge two specified populations.
# Inputs:
#  - Pop1:Population*
#  - Pop2:Population*
#  - Shuffle:BOOLEAN  (optional)
sub mergePopulations
{
	my %params = @_;

	croak "'Pop1' missing" unless exists $params{Pop1};
	croak "'Pop2' missing" unless exists $params{Pop2};
	
	my $pop1 = $params{Pop1};
	my $pop2 = $params{Pop2};
	if (ref $pop1 eq "Gnomish::Population" or
		ref $pop2 eq "Gnomish::Population")
	{
		croak "Population objects must be passed by reference";	
	}
	
	my $new_pop = new Gnomish::Population;
	
	# add the first population
	my ($i, $member);
	for ($i = 0; $i < $$pop1->{NMembers} ; $i++)
	{
		$member = $$pop1->getMemberRef (Member => $i);
		$new_pop->addMember (Genotype => $$member);
	}
	
	# add the second population
	for ($i = 0; $i < $$pop2->{NMembers}; $i++)
	{
		$member = $$pop2->getMemberRef (Member => $i);
		$new_pop->addMember (Genotype => $$member);
	}
	
	return $new_pop;
}

##### Population creation routines ########

# create a random population from scratch with random alleles 
sub makeRandomAllelePopulation
{
	my %params = @_;
	
	croak "'NMembers' missing" unless exists $params{NMembers};
	croak "'Genome' missing" unless exists $params{NMembers};
	
	my $genome = $params{Genome};
	if ( ref $genome ne "Gnomish::Genome" )
	{
		croak "Invalid genome object";
	}
	
	my $population = new Gnomish::Population;
	
	for ( my $i = 0; $i < $params{NMembers}; $i++)
	{
		my $genotype = new Gnomish::Genotype;
		$genotype->makeRandomAlleleGenotype ( Genome => $genome );
		$population->addMember ( Genotype => $genotype );
	}
	
	return $population;
}

# create a population with random homozygous alleles
sub makeRandomHomozygousPopulation
{
	my %params = @_;
	
	croak "'NMembers' missing" unless exists $params{NMembers};
	croak "'Genome' missing" unless exists $params{NMembers};
	
	my $genome = $params{Genome};
	if ( ref $genome ne "Gnomish::Genome" )
	{
		croak "Invalid genome object";
	}
	
	my $population = new Gnomish::Population;
	
	for ( my $i = 0; $i < $params{NMembers}; $i++)
	{
		my $genotype = new Gnomish::Genotype;
		$genotype->makeRandomHomozygousGenotype ( Genome => $genome );
		$population->addMember ( Genotype => $genotype );
	}
	
	return $population;
}

# create a population homozygous with specified allele
sub makeHomozygousPopulation
{
	my %params = @_;
	
	croak "'NMembers' missing" unless exists $params{NMembers};
	croak "'Genome' missing" unless exists $params{Genome};
	croak "'Allele' missing" unless exists $params{Allele};
	
	my $genome = $params{Genome};
	if ( ref $genome ne "Gnomish::Genome" )
	{
		croak "Invalid genome object";
	}
	
	my $allele = $params{Allele};
	if ( not $genome->validAllele ( Allele => $allele ) )
	{
		croak "Specified allele is not valid for the genome";
	}
	
	my $population = new Gnomish::Population;
	for (my $i = 0; $i < $params{NMembers}; $i++)
	{
		my $genotype = new Gnomish::Genotype;
		$genotype->makeHomozygousGenotype ( Genome => $genome,
					Allele => $allele );
		$population->addMember ( Genotype => $genotype );
	}
	
	return $population
}

##################### POPULATION CLASS DEF ########################

sub new
{
	my $class = shift;
	my %params = @_;
	
	my $self = {};
	bless $self, $class;
	
	if (exists $params{Name})
	{
		$self->{Name} = $params{Name};	
	}
	else
	{
		$self->{Name} = "";	
	}
	
	$self->{NMembers} = 0;
	my @population;
	$self->{_population} = \@population;

	return $self;	
}


# set the Population name
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
	
	if (not exists $self->{Name})
	{
		return 0;
		
	}
	return $self->{Name};
}

# Return the number of members in the population.
# This is here for consistency.
sub nMembers
{
	my $self = shift;
	
	return $self->{NMembers};
}


# Add a new genotype to the population.
# Input:
#  - Genotype:Genotype*
sub addMember
{
	my $self = shift;
	my %params = @_;

	croak "'Genotype' missing" unless exists $params{Genotype};
	
	my $genotype = $params{Genotype};
	if ( ref $genotype ne "Gnomish::Genotype" )
	{
		croak "Invalid Genotype object";	
	}
	
	my $population = $self->{_population};
	push @$population, \$genotype;
	
	# this is really only written this way for clarity...
	$self->{NMembers} = $self->{NMembers} + 1;
}


# Return a reference to a particular member of the population
# Input: Member
sub getMemberRef
{
	my $self = shift;
	my %params = @_;
	
	croak "'Member' missing" unless exists $params{Member};
	my $member = $params{Member};
	
	if ($member >= $self->{NMembers})
	{
		croak "Member index is out of bounds";	
	}
	
	my $pop = $self->{_population};
	return $$pop[$member];
}

sub getMember
{
	my $self = shift;
	my %params = @_;

	croak "'Member' missing" unless exists $params{Member};
	my $member = $params{Member};
	
	if ($member >= $self->{NMembers})
	{
		croak "Member index is out of bounds";	
	}

	my $pop = $self->{_population};
	my $memb = $$pop[$member];

	return $$memb; # de-reference the Genotype reference and return a copy 	
	
}

# Randomly shuffle the members in the population
sub shuffle
{
	my $self = shift;
	
	my $pop = $self->{_population};

	my ($i, $r, $cur, $dest);
	for ($i = 0; $i < $self->{NMembers}; $i++)
	{
		# pick a random member of the population to swap places with
		$r = int rand ($self->{NMembers});
		
		# save the original references into temp storage
		# to make sure they don't get lost		
		$cur  = $$pop[$i];
		$dest = $$pop[$r];
		
		# swap places
		$$pop[$i] = $dest;
		$$pop[$r] = $cur;
	}
	
	$self->{_population} = $pop;
}

# compute the phenotypic values for each of the genotypes in the population.
sub computePhenotypes
{
	my $self = shift;

	my $p = $self->{_population};
	my @pop = @$p;

	my $genotype;
	foreach $genotype (@pop)
	{
		$$genotype->computePhenotypicValues();	
	}
}



1;
