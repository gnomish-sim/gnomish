# RICE 
# This is a simulation of a simple cross
# SSD-derived poulation of 500 RILs segregating for four TRAITS,
# simple genetic model as suggested by Diane Mather
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



initialize();

### Setup the simulation ###

#  the genome
my $genome = new Gnomish::Genome;

# recombination distances. 
my $recombination_dists = createRecombinationDistanceString (
                                Chromosomes => 12,
                                ChromosomeLengths => "15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15",
                                RecombinationDistances => "0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1" );

$genome->makeGenome ( Name => "Rice",
                        Nloci => 180,
			RecombinationDists => $recombination_dists,
                        AlleleSet => "A,B" );

# TRT1 one QTL additive gene effect
my %trait_allele_effects;
$trait_allele_effects{25}{AA} = 1.0;
$trait_allele_effects{25}{BB} =  -1.0;
$trait_allele_effects{25}{AB} =  0;

# add to genome.
$genome->addTrait ( Name => "TRT1",
                        TraitMean => 100,
                        TraitVariance => 8,
                        TraitAlleleEffects => \%trait_allele_effects);

# TRT2 one QTL  greater effect and greater variance 
my %trait_allele_effects;
$trait_allele_effects{25}{AA} = 3.0;
$trait_allele_effects{25}{BB} =  -3.0;
$trait_allele_effects{25}{AB} =  0;

# add to genome.
$genome->addTrait ( Name => "TRT2",
                        TraitMean => 100,
                        TraitVariance => 25,
                        TraitAlleleEffects => \%trait_allele_effects);


#  TRT3 one QTL ADitive, greater variance than TRT2
my %trait_allele_effects;
$trait_allele_effects{25}{AA} = 3.0;
$trait_allele_effects{25}{BB} =  -3.0;
$trait_allele_effects{25}{AB} =  0;

# add to genome.
$genome->addTrait ( Name => "TRT3",
                        TraitMean => 100,
                        TraitVariance => 25,
                        TraitAlleleEffects => \%trait_allele_effects);

#  TRT4 one QTL, greater variance than TRT3 
my %trait_allele_effects;
$trait_allele_effects{25}{AA} = 3.0;
$trait_allele_effects{25}{BB} =  -3.0;
$trait_allele_effects{25}{AB} =  0;

# add to genome.
$genome->addTrait ( Name => "TRT4",
                        TraitMean => 100,
                        TraitVariance => 25,
                        TraitAlleleEffects => \%trait_allele_effects);


# Create the PARENTS
my $g0 = new Gnomish::Genotype;
$g0->makeHomozygousGenotype (Name => "PA",
                                Genome => $genome,
                                Allele => 'A');
                                
my $g1 = new Gnomish::Genotype;
$g1->makeHomozygousGenotype (Name => 'PB',
                                Genome => $genome,
                                Allele => 'B');
 

# HAI'S ADDITION
# To make the program loop n number of times, you can just use a 'for' loop that spans
# the part of the simulation that changes from run to run.
# So to make the simulation loop 50 times we can do:

for (my  $n = 0; $n < 50; $n++)
{
	print "Run: $n\n";		
						
	# cross PA with PB to generate the F1 population   
	my $f1 = new Gnomish::Genotype;
	$f1->makeWithParents (Parent1 => \$g0,
        	                Parent2 => \$g1);
                                
	# create the original f1 population
	my $p0 = new Gnomish::Population;
	$p0->addMember (Genotype => $f1);

	# clone F1 to get 500 individuals
	my $p1 = clonalProgeny (ParentPop =>\$p0, NProgeny => 500);

	#generate a population of 500 SSD lines
	my $population = singleSeedDescent (NGenerations => 4, ParentPop => \$p1);

	# compute the phenotypic values.
	$population->computePhenotypes();

	# visit each genotype in the population, one at a time
	
	# CODE ALTERED. Instead of printing the locus states of loci 0 to 179, we
	# print the values: 2 for state AA, 0 for state BB, 1 for state AB
	#for (my $i =0; $i<$population->nMembers(); $i++)
	#{
        #	my $genotype = $population->getMember (Member => $i);
	#
        #	# Get fragment of the genotypes genes starting at locus 1 to locus 1900
        #       my $fragment = $genotype->getFragment (StartLocus => 0,
        #       	                        StopLocus => 179);
	#
       	#	print "RIL$i $fragment\n";
	#}
	for (my $i =0; $i<$population->nMembers(); $i++)
	{
		print "RIL$i ";
		
        	my $genotype = $population->getMember (Member => $i);
		
		for (my $j = 0; $j < 180; $j++)
		{
			my $state = $genotype->getLocusState (Locus => $j);
			
			if ($state eq "AA")
			{
				print "2,";
			}
			elsif ($state eq "BB")
			{
				print "0,";
			}
			elsif ($state eq "AB" || $state eq "BA")
			{
				print "1,";
			}
			else
			{	# This is an invalid locus state, print a warning
				print "-1,";
			}
		}
		print "\n";
	}


	
	#get phenotypic value for TRAITs
	for (my $i =0; $i<$population->nMembers(); $i++)
	{
        	my $genotype = $population->getMember (Member => $i);

        	# get values
                	my $v1 = $genotype->getPhenotypicValueForTrait (TraitName => "TRT1");
                	my $v2 = $genotype->getPhenotypicValueForTrait (TraitName => "TRT2");
                	my $v3 = $genotype->getPhenotypicValueForTrait (TraitName => "TRT3");
                	my $v4 = $genotype->getPhenotypicValueForTrait (TraitName => "TRT4");
                	my $g1 = $genotype->getGenotypicValueForTrait (TraitName => "TRT1");
                	my $g2 = $genotype->getGenotypicValueForTrait (TraitName => "TRT2");
                	my $g3 = $genotype->getGenotypicValueForTrait (TraitName => "TRT3");
                	my $g4 = $genotype->getGenotypicValueForTrait (TraitName => "TRT4");

       		print "RIL$i $v1 $v2 $v3 $v4 $g1 $g2 $g3 $g4 \n";
       
	}


	# compute the variance of the traits and report
	my $variance_trt1 = computeTraitPhenotypeVariance (Population => \$population,
                                                    	TraitName => "TRT1");
	my $variance_trt2 = computeTraitPhenotypeVariance (Population => \$population,
        	                                            TraitName => "TRT2");
	my $variance_trt3 = computeTraitPhenotypeVariance (Population => \$population,
        	                                            TraitName => "TRT3");
	my $variance_trt4 = computeTraitPhenotypeVariance (Population => \$population,
        	                                            TraitName => "TRT4");


    	print "Variance: $variance_trt1 $variance_trt2 $variance_trt3 $variance_trt4\n";
}
