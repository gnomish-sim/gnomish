# BarleyHT tutorial
# In this tutorial, we want to simulate the breeding of two populations,
# one homozygous for the 'H' allele, and one homozygous for the 'T' allele.
#
# Author: Hai Pham, hpham42@gmail.com
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


use strict;

use Gnomish::Genome;
use Gnomish::Population;
use Gnomish::Genotype;
use Gnomish::Stdlib;
use Gnomish::Evaluations;


# Initialize random number generator
initialize();

### Setup the simulation ###

# First, create the genome
my $genome = new Gnomish::Genome;

# Create a string with the recombination distances. 
my $recombination_dists = createRecombinationDistanceString (
				Chromosomes => 7,
				ChromosomeLengths => "161,214,155,165,237,162,187",
				RecombinationDistances => "0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01" );

$genome->makeGenome ( Name => "BarleyHT",
			Nloci => 1281,
			RecombinationDists => $recombination_dists,
			AlleleSet => "H,T" );

# Next, we define a new Trait
my %trait_allele_effects;
$trait_allele_effects{717}{HH} = -1.1;
$trait_allele_effects{717}{TT} =  1.1;
$trait_allele_effects{717}{HT} =  0;
$trait_allele_effects{806}{HH} =  -1;
$trait_allele_effects{806}{TT} =  1;
$trait_allele_effects{806}{HT} =  0;
$trait_allele_effects{1190}{HH} = -1.3;
$trait_allele_effects{1190}{TT} = 1.3;
$trait_allele_effects{1190}{HT} = 0;

# Add it to the genome's definition.
$genome->addTrait ( Name => "Height",
                        TraitMean => 89,
                        TraitVariance => 5,
                        TraitAlleleEffects => \%trait_allele_effects);

# Create the genotypes
my $g0 = new Gnomish::Genotype;
$g0->makeHomozygousGenotype (Name => "Harrington",
                                Genome => $genome,
                                Allele => 'H');
                                
my $g1 = new Gnomish::Genotype;
$g1->makeHomozygousGenotype (Name => 'TR306',
                                Genome => $genome,
                                Allele => 'T');
 
# cross Harrington with TR306 to create the F1 population   
my $f1 = new Gnomish::Genotype;
$f1->makeWithParents (Parent1 => \$g0,
                        Parent2 => \$g1);
                                
# create the original f1 population
my $p0 = new Gnomish::Population;
$p0->addMember (Genotype => $f1);

# ... and generate 150 doubled haploids from the F1
my $population = doubledHaploids (NProgeny => 150,
              	                  ParentPop => \$p0);

# compute the phenotypic values.
$population->computePhenotypes();

# compute the mean of the trait "Height" and report
my $mean_height = computeTraitPhenotypeMean ( Population => \$population,
						TraitName => "Height" );
print "Mean Height: $mean_height\n";

# compute the variance of the trait "Height" and report
my $height_variance = computeTraitPhenotypeVariance (Population => \$population,
							 TraitName => "Height");
print "Height Variance: $height_variance\n";


# Find the number of doubled haploids with allele state 'HH'
# at locus 1175 and locus 1195.

# Get number of individuals in the population
my $members = $population->nMembers();  

# We will use this array to keep track of the number of members
# of the population with locus state HH at locus 1175 and 1195 
my @mlist;

# Iterate over the population
for (my $i = 0; $i < $members; $i++)
{
	# Get a copy of the individual from the population for us to examine
	my $m = $population->getMember ( Member => $i );
	
	# Get it's locus state at locus 1175 and 1195
	my $locus_state1175 = $m->getLocusState ( Locus => 1175 );
	my $locus_state1195 = $m->getLocusState ( Locus => 1195 );
	
	# Check to see if the locus state at both locus is 'HH'
	if ( ($locus_state1175 eq "HH") and ($locus_state1195 eq "HH") )
	{
		# if the locus state is HH at both locus, save the individual to a list
		push @mlist, $m;
	}
}

# Count the number of individuals found with the locus state
# 'HH' at locus 1175 and 1195, and report.
my $n = scalar @mlist;
print "Number of doubled happloids with two H alleles at locus 1175 and 1195: $n\n";

# Compute the mean height of the individuals found
my $sum = 0;
foreach my $m (@mlist)
{
	my $height = $m->getPhenotypicValueForTrait ( TraitName => "Height" );
	$sum = $sum + $height;
}
my $mean = $sum / $n;
print "Mean height: $mean\n";


# Next, find the number of doubled haploids with
# two T alleles at locus 1175 AND locus 1195

# Reset the member list, throwing out the old HH doubled haploids..
undef @mlist; 

for (my $i = 0; $i < $members; $i++)
{
	my $m = $population->getMember ( Member => $i );
	
	my $locus_state1175 = $m->getLocusState ( Locus => 1175 );
	my $locus_state1195 = $m->getLocusState ( Locus => 1195 );
	
	# Check and see if the state at both locus is "TT"
	if ( ($locus_state1175 eq "TT") and ($locus_state1195 eq "TT") )
	{
		push @mlist, $m;
	}
}

# Count the number of individuals found with the locus state 'TT' at locus 1175 and 1195
$n = scalar @mlist;
print "Number of doubled haploids with two T alleles at locus 1175 and locus 1195: $n\n";

# Compute the mean height of the individuals found.
$sum = 0;
foreach my $m (@mlist)
{
	my $height = $m->getPhenotypicValueForTrait ( TraitName => "Height" );
	$sum = $sum + $height;
}
$mean = $sum / $n;
print "Mean height: $mean\n";


