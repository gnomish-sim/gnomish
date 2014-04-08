# Hardy-Weinberg Equilibrium test
# Hai Pham
# Last Changed: Mar 11, 2003
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

# Load the Perl module that enforces strict syntax checking
# and scoping rules.
use strict;

# Load Gnomish Perl modules.
use Gnomish::Population;
use Gnomish::Genotype;
use Gnomish::Genome;
use Gnomish::Stdlib;
use Gnomish::Evaluations;


# Initialize the random number generator.
# Make sure that this is executed at the
# beginning of every script!
initialize();

# First, create a genome object.
my $genome = new Gnomish::Genome;

# Define a genome with a single locus.
# Since there is only one locus, the recombination 
# distance array is empty
$genome->makeGenome ( Name => "HW",
			 Nloci => 1,
			 RecombinationDists => "",
			 AlleleSet => "a,b" );

# Make a population of 100 genotypes that are completely
# homozygous for the 'a' allele
my $p0 = makeHomozygousPopulation (NMembers => 100,
				Genome => $genome,
				Allele => "a");

# Make a population of 900 genotypes that are completely
# homozygous for the 'b' allele
my $p1 = makeHomozygousPopulation (NMembers => 900,
			Genome => $genome,
			Allele => "b");
			
# Create a merged population of these two
my $f1 = mergePopulations (Pop1 => \$p0, Pop2 => \$p1);

# Subject the merged population to one generation of random mating,
# creating a new population of 1000 members that replaces the original
# parent (F1) population.
my $f2 = randomMate (ParentPop => \$f1, NProgeny => 1000, NGenerations => 1);

# Count the occurrence of allele states 'aa', 'ab', and 'bb'
my %allele_states_count = countAlleleStatesInPopulation (
				Population => \$f2 );

# Get the number of members found in the '$f2' population.
my $members = $f2->nMembers();

# Calculate the frequencies of each allele state and display
print "Population Members: $members\n";
foreach my $state (keys %allele_states_count)
{
	my $freq = $allele_states_count{$state} / $members;
	print "Allele state: $state, total: $allele_states_count{$state}, frequency: $freq\n";
}


