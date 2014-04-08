# Gnomish class object Genotype
# This object defines the genotypic properties of an individial in a population
#
# File: Genotype.pm
# Author: Hai Pham, hpham42@gmail.com
# Last changed: Feb 18, 2003
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


# NOTE: This class was implemented assuming that each Gamete would have 1 allele string (diploid organism)
#       This is true of most organisms, but there are lots of plants that are polyploid!
#       We may want to spend some effort extending this class to understand polyploid organisms
#       in the future.


package Gnomish::Genotype;

use strict;
use Carp;

use Gnomish::Stdlib;

my $DEBUG = 0;

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
	
	return $self;	
}

# create a new genotype with homozyous alleles
sub makeHomozygousGenotype
{
	my $self = shift;
	my %params = @_;
	
	croak "'Allele' missing" unless exists $params{Allele};
	croak "'Genome' missing" unless exists $params{Genome};
	
	if ( (ref $params{Genome}) ne "Gnomish::Genome")
	{
		croak "Genome object must be passed by value";                      
        }
	my $genome = $params{Genome};
        $self->{_Genome} = \$genome;
	if ( !$genome->validAllele(Allele => $params{Allele}) )
	{
		croak "The allele $params{Allele} is not valid for the specified genome";	
	}
	
	if (exists $params{Name})
	{
		$self->{Name} = $params{Name};	
	}		
	
	$self->_createAlleleStringsWithValue ($params{Allele});
}


# create a new genotype with a self-crossed parent.
sub makeWithSelfCrossedParent
{
	my $self = shift;
	my %params = @_;
	
	croak "'Parent' missing" unless exists $params{Parent};
	
	if (ref $params{Parent} eq "Gnomish::Genotype")
	{
		croak "Parent object must be passed by reference";
	}
	
	$self->_makeDoubledHaploidGenotype ($params{Parent});
}


# create new genotype with two parents.
sub makeWithParents
{
	my $self = shift;
	my %params = @_;
	
	croak "'Parent1' missing" unless exists $params{Parent1};
	croak "'Parent2' missing" unless exists $params{Parent2};
	
	if ( (ref $params{Parent1} eq "Gnomish::Genotype") or
		(ref $params{Parent2} eq "Gnomish::Genotype") )
	{
		croak "Parent objects must be passed by reference";
	}
	
	$self->_makeGenotypeWithParents ($params{Parent1}, $params{Parent2});
	
}

# create new genotype by specifying the allele strings and genome
sub makeToSpec
{
	my $self = shift;
	my %params = @_;
	
	croak "'Genome' missing" unless exists $params{Genome};
	croak "'AlleleString0' missing" unless exists $params{AlleleString0};
	croak "'AlleleString1' missing" unless exists $params{AlleleString1};
	
	if ( ref $params{Genome} ne "Gnomish::Genome" )
	{
		croak "Genome object must be passed by value";
	}
	
	my $genome = $params{Genome};
	$self->{_Genome} = \$genome;
	
	my $genotypic_values;
	if (defined $params{GenotypicValues})
	{
		$genotypic_values = $params{GenotypicValues};
	}
	
	my $phenotypic_values;
	if (defined $params{PhenotypicValues})
	{
		$phenotypic_values = $params{PhenotypicValues};
	}
	
	$self->_makeGenotypeFromSpecs ($params{AlleleString0},
					$params{AlleleString1},
					$genotypic_values,
					$phenotypic_values);
}

# Different sorts of MakeGenotypes:
# - make from parents
# - self-crossed parent (make from single parent)
# - create with a 'create method'
# - create from specs (given genome and allele strings)




# Create a new Genotype either from thin-air or with two parents.
# Inputs (creation with a self-crossed parent -- generating a doubled haploid progeny)
# - Parent:Genotype*
#
# Inputs (creation with two parents):
# - Parent1:Genotype*
# - Parent2:Genotype*
#
# Inputs (create from thin-air)
# - CreateMethod:String (one of the 6 approved methods) 
# - Genome:Genome*
#sub makeGenotype
#{
#	my $self = shift;
#	my %params = @_;

#	croak "'Genome' missing" unless exists $params{Genome};
	
#	if ( (ref $params{Genome}) ne "Gnomish::Genome")
#	{
#		croak "Genome object must be passed by value";			
#	}
#	
#	my $genome = $params{Genome};
#	$self->{_Genome} = \$genome;
#	$self->_makeGenotypeFromThinAir ($params{CreateMethod});
#	
#}

# create new genotype with randomly generated allele strings
sub makeRandomAlleleGenotype
{
	my $self = shift;
	my %params = @_;
	
	croak "'Genome' missing" unless exists $params{Genome};
	
	my $genome = $params{Genome};
	if ( ref $genome ne "Gnomish::Genome" )
	{
		croak "Invalid genome object";
	}
	
	$self->{_Genome} = \$genome;
	$self->_createRandomAlleleStrings();
}

# initialize the genotype with homozygous alleles that are randomly chosen at each locus
sub makeRandomHomozygousGenotype
{
	my $self = shift;
	my %params = @_;
	
	croak "'Genome' missing" unless exists $params{Genome};
	
	my $genome = $params{Genome};
	if ( ref $genome ne "Gnomish::Genome" )
	{
		croak "Invalid genome object";
	}
	
	$self->{_Genome} = \$genome;
	$self->_createRandomHomozygousAlleleStrings();
}

##### UTILITY methods ########

sub getName
{
	my $self = shift;

	if (exists $self->{Name})
	{
		return $self->{Name};	
	}
	
	return;
}

# return the genome name that the genotype belongs to
sub getGenomeName
{
	my $self = shift;
	
	my $genome = $self->{_Genome};
	my $name = $$genome->getName();

	return $name;
}

# return the Genotype's genome
sub getGenome
{
	my $self = shift;

	my $genome = $self->{_Genome};

	return $$genome;	
}

# return a reference to the genotype's genome
sub getGenomeRef
{
	my $self = shift;
	
	return $self->{_Genome};	
}


# Go through the genotype and compute a sum at each loci.
# Eg. the output from a genotype with Allele strings:
#    00111001
#    11110001
# Sum:
#    11221002
# The result is returned in an array.  
#sub sumLoci
#{
#	my $self = shift;
#
#	my $a_str0 = $self->{_AlleleStrings}[0];
#	my $a_str1 = $self->{_AlleleStrings}[1];
#	
#	my $genome = $self->{_Genome};
#	my $nloci = $$genome->getNloci ();
#
#	my ($i, @sums);
#	for ($i = 0; $i < $nloci; $i++)
#	{
#		push @sums, ($$a_str0[$i] + $$a_str1[$i]);
#	}
#
#	return @sums;
#}


# return the locus state at the particular locus...
sub getLocusState
{
	my $self = shift;	
	my %params = @_;
	
	croak "'Locus' missing" unless exists $params{Locus};
	my $locus = $params{Locus};
	
	my $a_str0 = $self->{_AlleleStrings}[0];
	my $a_str1 = $self->{_AlleleStrings}[1];
	
	
	my $genome = $self->{_Genome};
	my $nloci = $$genome->getNloci ();
	if ( $locus >= $nloci )
	{
		croak "'Locus' out of range";	
	}
	
	my $state = $$a_str0[$locus] . $$a_str1[$locus];
	
	return $state;
}

# get a fragment of the genotype's genes
sub getFragment
{
	my $self = shift;
	my %params = @_;
	
	croak "'StartLocus' missing" unless exists $params{StartLocus};
	croak "'StopLocus' missing" unless exists $params{StartLocus};
	
	my $start = $params{StartLocus};
	my $stop = $params{StopLocus};
	
	my ($i, @fragment);
	for ($i = $start; $i <= $stop; $i++)
	{
		push @fragment, $self->getLocusState (Locus => $i);
	}
	
	my $s = join ',', @fragment;
	return $s;
}


# Return the specified gamete
# Input:
#  - Gamete:Int (valid value is 0 or 1)
sub getGamete
{
	my $self = shift;
	my %params = @_;
	
	croak "'Gamete' missing" unless exists $params{Gamete};
	
	my $gamete = $params{Gamete};
	if ($gamete < 0 or $gamete > 1)
	{
		croak "Gamete value must be either '0' or '1'";	
	}
	
	my $g = $self->{_AlleleStrings}[$gamete];
	my $gs = join ',', @$g;
		
	return $gs;
}

# return the specified gamete in an array
sub getGameteAsArray
{
	my $self = shift;
	my %params = @_;
	
	croak "'Gamete' missing" unless exists $params{Gamete};
	
	my $gamete = $params{Gamete};
	if ($gamete < 0 or $gamete > 1)
	{
		croak "Gamete value must be either '0' or '1'";	
	}
	
	my $g = $self->{_AlleleStrings}[$gamete];
	
	return @$g;
}

# draw an ascii representation of the genotype
sub drawInAscii
{
	my $self = shift;

	my $gnome_name = $self->getGenomeName();
	my $allele_strings = $self->{_AlleleStrings};
		
	my ($allele_string, $allele, $i);
	$i = 0;
	foreach $allele_string (@$allele_strings)
	{
	#	print "#Allele $i\n";
		foreach $allele (@$allele_string)
		{
			if (defined $allele)
			{
				print $allele;
			}
		}
		$i++;
		print "\n";
	}
}


sub computeGenotypicValues
{
	my $self = shift;

	my $genome = $self->{_Genome};
	my @tnames = $$genome->getTraitNames();
	
	my ($t, $trait);
	my ($effect_sum, $effect_val);
	my (@loci, $locus, $locus_state);
	my ($g_value, %genotypic_values);
	foreach $t (@tnames)
	{
		$trait = $$genome->getTraitRef (Name => $t);
		
		# get array of loci that affect the trait
		@loci = $$trait->getActiveLoci();
		$effect_sum = 0;
		
		# now that we have the list of loci that affect the trait, go down the
		# list and see what their genotypic effect values are
		foreach $locus (@loci)
		{
			$locus_state = $self->getLocusState (Locus => $locus);
			$effect_val = $$trait->getLocusEffectValue ( Locus => $locus, 
								AlleleState => $locus_state );
			$effect_sum += $effect_val;
		}
		
		$g_value = $effect_sum + $$trait->getTraitMean();
		$genotypic_values{$t} = $g_value;
	}
	
	$self->{_GenotypicValues} = \%genotypic_values;
}

sub _getGenotypicValues
{
	my $self = shift;

	my $gv = $self->{_GenotypicValues};
	return %$gv;	
}

sub getGenotypicValueForTrait
{
	my $self = shift;
	my %params = @_;

	croak "'Name' missing" unless exists $params{TraitName};
	my $name = $params{TraitName};
		
	my $gv = $self->{_GenotypicValues};
	croak "No value for '$name'" unless exists $$gv{$name};
	
	return $$gv{$name};	
}


sub computePhenotypicValues
{
	my $self = shift;
	
	# make sure the genotypic values have been computed first
	if ( not exists $self->{_GenotypicValues} )
	{
		$self->computeGenotypicValues();	
	}
	
	my $genotypic_values = $self->{_GenotypicValues};
	my $genome = $self->{_Genome};
	my @tnames = $$genome->getTraitNames();
	
	my ($trait_name, $trait);
	my ($trait_variance, $rand_variance);
	my (%phenotypic_values, $p_value);
	my ($g_value);
	foreach $trait_name (@tnames)
	{
		$trait = $$genome->getTraitRef ( Name => $trait_name );
		$trait_variance = $$trait->getTraitVariance();
		$rand_variance = gaussianRand ( Variance => $trait_variance );
		$g_value = $$genotypic_values{$trait_name};
		$p_value = $g_value + $rand_variance;
		$phenotypic_values{$trait_name} = $p_value;
	}

	$self->{_PhenotypicValues} = \%phenotypic_values;	
}


sub _getPhenotypicValues
{
	my $self = shift;

	my $pv = $self->{_PhenotypicValues};
	return %$pv;	
}

sub getPhenotypicValueForTrait
{
	my $self = shift;
	my %params = @_;

	croak "'TraitName' missing" unless exists $params{TraitName};
	my $name = $params{TraitName};
		
	my $pv = $self->{_PhenotypicValues};
	croak "No value for '$name'" unless exists $$pv{$name};
	
	return $$pv{$name};	
}

# return the name of of the traits in the phenotypic & genotypic values
sub getTraitNames
{
	my $self = shift;
	
	my $pv = $self->{_PhenotypicValues};

	return (keys %$pv);	
}





#### Private class methods ##########

# create a new random genotype.
# Params:
#    - CreateMethod
# "CreateMethod" is an int or a string identifier that specifies one of the following
# CreateMethod:Int  (the method used to create the specific genotype)
# 1 = Random allele at each locus: we randomly assign 1's and 0's to each allele in the string
# 2 = Random homozygous pair at each locus: we randomly assign 1's and 0's at each locus (the allele pair are assigned the same value)
# 3 = 'solid' alleles at every locus: assign 1's to every allele in the genome
# 4 = 'dotted' alleles at every locus: assign 0's to every allele in the genome 
# 5 = heterozygous, complete coupling: we randomly assign all alleles in one string to be either
#                                      1 or 0, and we set the other allele string to it's opposite value.
#                                      ie. we randomly make the first allele string to be either 
#                                      all 1's or all 0's, and then we set the other allele 
#                                      string to it's opposite value.
# 6 = heterozygous, random coupling: at every locus, we randomly assign 1 or 0 to one allele and 
#                                    we set the other to the opposite value. ie. if we assign a 1 
#                                    to the first allele then we assign 0 to the other.
# CreateMethod: String (this function can also be activated with a string create method identifier
# "RANDOM_ALLELE" = 1
# "RANDOM_HOMOZYGOUS" = 2
# "SOLID_ALLELES" = 3
# "DOTTED_ALLELES" = 4
# "HETEROZYGOUS_COMPLETE_COUPLING" = 5
# "HETEROZYGOUS_RANDOM_COUPLING" = 6
sub _makeGenotypeFromThinAir
{
	my $self = shift;
	my $create_method = shift;

	if ($create_method eq "RANDOM_ALLELE" || $create_method eq "1")
	{
		# create genotype with random allele
		$self->_createRandomAlleleStrings();
	}
	elsif ($create_method eq "RANDOM_HOMOZYGOUS" || $create_method eq "2")
	{
		# create genotype with random hymozygous alleles
		$self->_createRandomHomozygousAlleleStrings();
	}
	elsif ($create_method eq "SOLID_ALLELES" || $create_method eq "3")
	{
		# create genotype with alleles set to '1'
		die "create method currently undefined";
		$self->_createAlleleStringsWithValue(1);
	}
	elsif ($create_method eq "DOTTED_ALLELES" || $create_method eq "4")
	{
		# create genotype with alleles set to '0'
		die "create method currently undefined";
		$self->_createAlleleStringsWithValue(0);
	}
	elsif ($create_method eq "HETEROZYGOUS_COMPLETE_COUPLING" || $create_method eq "5")
	{
		$self->_createHeterozygousAlleleStringsWithCompleteCoupling();	
	}
	elsif ($create_method eq "HETEROZYGOUS_RANDOM_COUPLING" || $create_method eq "6")
	{
		$self->_createHeterozygousAlleleStringsWithRandomCoupling();	
	}
	else
	{
		croak "Invalid CreateMethod in creation of new genotype!";	
	}
}

# create genotype with random allele at each locus: we randomly
# assign 1's and 0's to each allele in the string
sub _createRandomAlleleStrings
{
	my $self = shift;
	
	my $genome = $self->{_Genome};
	my $nloci = $$genome->getNloci;
	
	my @allele_set = split ',', $$genome->getAlleleSet();
	my $alleles = scalar @allele_set;
	
	my ($i, $a, $allele);
	for ($a = 0; $a < 2; $a++)
	{
		my @allele_string;
		for ($i = 0; $i < $nloci; $i++)
		{
			$allele = int rand($alleles);
			push @allele_string, $allele_set[$allele];
		}
		$self->{_AlleleStrings}[$a] = \@allele_string;
	}
}

# Create a new genotype with random homozygous allele strings.
# At each locus, randomly assign a 0 or a 1
sub _createRandomHomozygousAlleleStrings
{
	my $self = shift;
	
	my $genome = $self->{_Genome};
	my $nloci  = $$genome->getNloci();

	my @allele_set = split ',', $$genome->getAlleleSet();
	my $alleles = scalar @allele_set;
	
	# first we need to create some null array references
	my (@allele_str0, @allele_str1);
	my ($i, $a);
	for ($i = 0; $i < $nloci; $i++)
	{
		$a = int rand($alleles);
		push @allele_str0, $allele_set[$a];
		push @allele_str1, $allele_set[$a];
	}
	$self->{_AlleleStrings}[0] = \@allele_str0;
	$self->{_AlleleStrings}[1] = \@allele_str1;
}

# create Allele strings with all alleles set to the given value
sub _createAlleleStringsWithValue
{
	my $self = shift;
	my $value = shift;
	
	my $genome = $self->{_Genome};
	my $nloci  = $$genome->getNloci();
	
	if ( not $$genome->validAllele( Allele => "$value") )
	{
		croak "The allele '$value' is not valid for the genotype's Genome";
	}
	
	my ($i, $a);
	for ($a = 0; $a < 2; $a++)
	{
		my @allele_string;
		for ($i = 0; $i < $nloci; $i++)
		{
			push @allele_string, $value;
		}
		$self->{_AlleleStrings}[$a] = \@allele_string;
	}
}

# Create a new Genotype with heterozygous, complete coupling: we randomly assign all
# alleles in one string to be either 1 or 0, and we set the other allele string 
# to it's opposite value.  ie. we randomly make the first allele string to be either 
# all 1's or all 0's, and then we set the other allele string to it's opposite value.
# NOTE: this function assumes that there are ONLY TWO ALLELE STRINGS!
#       This is because it's behavior is undefined when there are more than 2 allele strings
sub _createHeterozygousAlleleStringsWithCompleteCoupling
{
	
	my $self = shift;
	
	my $genome = $self->{_Genome};
	my $nloci  = $$genome->getNloci();
	
	my @allele_set = split ',', $$genome->getAlleleSet();
	if ( (scalar @allele_set) != 2)
	{
		croak "'HETEROZYGOUS_COMPLETE_COUPLING' rule only works with genomes that have 2 alleles";
	}
	
	die;
	
	# define values for strings 1 & 2
	my $a1 = int rand(2); # randomly choose 1 or 0 for string1
	my $a2; 
	if ($a1 == 0) # assign the opposite value to string2
	{
		$a2 = 1;
	}
	else
	{
		$a2 = 0;	
	}
	
	my ($i, @string1, @string2);
	for ($i = 0; $i < $nloci; $i++)
	{
		push @string1, $a1;
		push @string2, $a2;
	}
	
	$self->{_AlleleStrings}[0] = \@string1;
	$self->{_AlleleStrings}[1] = \@string2;	
}

# create allele strings with heterozygous, random coupling: at every locus, we randomly 
# assign 1 or 0 to one allele and we set the other to the opposite value. 
# ie. if we assign a 1 to the first allele then we assign 0 to the other.
sub _createHeterozygousAlleleStringsWithRandomCoupling
{
	my $self = shift;
	
	die "the behaviour of this function is currently undefined!!";

	my $genome = $self->{_Genome};
	my $nloci  = $genome->getNloci();
	
	my ($i, $a1, $a2, @string1, @string2);
	for ($i = 0; $i < $nloci; $i++)
	{
		$a1 = int rand(2);
		if ($a1 == 0)
		{
			$a2 = 1;	
		}
		else
		{
			$a2 = 0;	
		}
		
		push @string1, $a1;
		push @string2, $a2;
	}
	
	$self->{_AlleleStrings}[0] = \@string1;
	$self->{_AlleleStrings}[1] = \@string2;	
}


# create a new doubled haploid progeny from one parent
# Input:
# - Parent:Genotype*
sub _makeDoubledHaploidGenotype
{
	my $self = shift;
	my $parent = shift;

	# define the genome of the progeny
	$self->{_Genome} = $$parent->{_Genome};

	# compute the gamete that is contibuted by the parent
	my $gamete = $self->_createNewGamete ($parent);
	
	# assign the same gamete to both haploids in the genotype
	$self->{_AlleleStrings}[0] = $gamete;
	$self->{_AlleleStrings}[1] = $gamete;
}


# create a new progeny from two parents.
# Input:
#  - Parent1
#  - Parent2
sub _makeGenotypeWithParents
{
	my $self = shift;
	my ($parent1, $parent2) = @_;
	
	# define the genome of the progeny
	$self->{_Genome} = $$parent1->{_Genome};

	# compute the allele string that is contributed by parent1
	$self->{_AlleleStrings}[0] = $self->_createNewGamete($parent1); 	

	# compute the allele string that is contributed by parent2
	$self->{_AlleleStrings}[1] = $self->_createNewGamete($parent2);
}

# Compute the allele sring that is contributed by the given parent to the progeny.
# NOTE: this function assumes that the genome of the parents only have 2 allele strings!
sub _createNewGamete
{
	my $self = shift;
	my $parent = shift;

	my $genome = $self->{_Genome};
	my $recomb_dists = $$genome->getRecombinationDistancesRef();
	my $nloci = $$genome->getNloci();
	
	my $allele_strs = $$parent->{_AlleleStrings};
	
	# first decide the allele string to start on
	my $str = int rand(2);	
	my @new_str;
	
	# the first loci is always copied
	push @new_str, $$allele_strs[$str][0];
	
	# now generate the new allele string
	my ($i, $rdist, $switch_prob, $r);
	my $loci = 1;
	foreach $rdist (@$recomb_dists)
	{
		$r = rand(1);
		if ($r < $rdist)
		{
			# switch the string to read from
			if ($str == 0)
			{
				$str = 1;
			}
			else
			{
				$str = 0;
			}
		}
		
		push @new_str, $$allele_strs[$str][$loci];
		$loci++;
	}
	
	# sanity check here for debugging purposes
	die unless ($loci == $nloci);
	
	return \@new_str;
}


sub _makeGenotypeFromSpecs
{
	my $self = shift;
	my ($a_str0, $a_str1, $genotypic_values, $phenotypic_values) = @_;

	# delete spaces from the strings
	$a_str0 =~ tr/ //d;
	$a_str1 =~ tr/ //d;		
	
	croak "Allele string lenghts do not match" 
		unless ( length ($a_str0) == length ($a_str1) );
	
	my @allele_string0 = split ',', $a_str0;
	my @allele_string1 = split ',', $a_str1;
	
	if (defined $genotypic_values)
	{
		if (not ref $genotypic_values)
		{
			croak "Genotypic values hash must be passed by reference";
		}
		
		$self->{_GenotypicValues} = $genotypic_values;	
	}
	
	if (defined $phenotypic_values)
	{
		if (not ref $phenotypic_values)
		{
			croak "Phenotypic values hash must be passed by reference";
		}
		
		$self->{_PhenotypicValues} = $phenotypic_values;
	}

	$self->{_AlleleStrings}[0] = \@allele_string0;
	$self->{_AlleleStrings}[1] = \@allele_string1;
}

1;
