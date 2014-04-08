# Gnomish class object Genome. 
# This class defines the genomic properties of the organism. eg. oats
#
# File: Genome.pm
# Author: Hai Pham, hpham42@gmail.com
# Last changed: Feb 19, 2003
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


package Gnomish::Genome;

use strict;
use Carp;

use Gnomish::Trait;
use Gnomish::Marker;

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
	else 
	{
		$self->{Name} = "";	
	}
	
	my @recomb;
	$self->{_RecombinationDists} = \@recomb;
	
	my %traits;
	$self->{_Traits} = \%traits;
	
	my %markers;
	$self->{_Markers} = \%markers;
	
	$self->{Nloci} = 0;
	
	$self->{_GenomeID} = 0;

	# if not specified, the values '0', '1' are assumed
	$self->{_AlleleSet} = "0,1";	
	
	return $self;
}

# initialize Genome object with provided parameters.
# Input:
# - Name:String (optional) 
# - Nloci: Int
# - RecombinationDists => "ref_to_array_of_Nloci-2_recombination_distances"
sub makeGenome
{
	my $self = shift;
	my %params = @_;
	my $tmp;

	if ($DEBUG)
	{
		print "-- Gnomish::Genome::makeGenome -- \n";
		print "Creating GenomeName: $params{GenomeName}\n";
	}

	croak "'Name' missing" unless exists $params{Name};
	$self->{Name} = $params{Name};
	
	croak "'Nloci' missing" unless exists $params{Nloci};
	my $nloci = $params{Nloci};
	$self->{Nloci} = $nloci;

	# now make a local copy of the RecombinationDists array.
	croak "'RecombinationDists' missing" unless exists $params{RecombinationDists};
				my $rdists = $params{RecombinationDists};
	$rdists =~ s/[\s\n\r\t]//g; # remove all spaces & return characters. 
	
	my @recomb_dists = split ',', $rdists;
	$self->{_RecombinationDists} = \@recomb_dists;
	
	# check to make sure the number of combination distatances matches nloci
	$rdists = scalar @recomb_dists;
	if ($rdists != ($nloci - 1))
	{
		die "Recombination Distances '$rdists' does not match nloci '$nloci'";	
	}
	
	
	croak "'AlleleSet' missing" unless exists $params{AlleleSet};
	# set the Genome's allele character set.
	# the set must be specified in a comma seperated string
        my $set = $params{AlleleSet};
	$set =~ s/[\s\n\r\t]//g; 	# strip out all spaces from the string
        $self->{_AlleleSet} = $set;
}

# initialize the genome with uniform recombination distances.
# Inputs:
#  - NlociPerChromosome:Int
#  - NChromosomes:Int
#  - RecombinationDist:Real
#  - Name:String (optional)
sub createUniformGenome
{
	my $self = shift;
	my %params = @_;

	croak "'NlociPerChromosome' missing" unless exists $params{NlociPerChromosome};
	croak "'NChromosomes' missing" unless exists $params{NChromosomes};
	croak "'RecombinationDist' missing" unless exists $params{RecombinationDist};
	croak "'AlleleSet' missing" unless exists $params{AlleleSet};	
	croak "'Name' missing" unless exists $params{Name};
	
	my $nloci_per_chr = $params{NlociPerChromosome};
	my $n_chromosomes = $params{NChromosomes};
	my $rdist         = $params{RecombinationDist};
	my $a_set         = $params{AlleleSet};
	my $g_name        = $params{Name};
	
	# create a new RecombinationDists array and fill it with the given rdist
	my ($i, $c, $nloci);
	my @recomb_dists;
	$nloci = 0;
	for ($c = 0; $c < $n_chromosomes; $c++)
	{
		for ($i = 0; $i < ($nloci_per_chr - 1); $i++)
		{
			push @recomb_dists, $rdist;
			$nloci++;
		}
		push @recomb_dists, 0.5;
		$nloci++;
	}
	# remove the last 0.5 rdist from the recombination dists array, since it's unecessary
	pop @recomb_dists;
	
	my $r_dists = join ',', @recomb_dists;
	
	# create the genome
	$self->makeGenome (Nloci => $nloci,
				Name => $g_name,
				RecombinationDists => $r_dists,
				AlleleSet => $a_set);
}


# initialize genome with a random RecombinationDist array for debugging/testing purposes
# Inputs:
# - Nloci:Int
# - MaxChromosomeLen:Int
sub createRandomGenome
{
	my $self = shift();
	my %params = @_;
	
	croak "'createRandomGenome' is currently disabled.";
	
	my ($nloci, $genome_name, $max_chromosome_len);
	my ($i, $chromosome_len, $rdist, @recomb_dists);

	croak "'Nloci' missing" unless exists $params{Nloci};
	croak "'MaxChromosomeLen' missing" unless exists $params{MaxChromosomeLen};
	croak "'AlleleSet' missing" unless exists $params{AlleleSet};
	croak "'Name' missing" unless exists $params{Name};

	my ($g_name) = $params{Name};
	my ($a_set)  = $params{AlleleSet};
			      
	$nloci       = $params{Nloci};
	$max_chromosome_len   = $params{MaxChromosomeLen};			      

	for ($i = 0; $i < ($nloci - 1); $i++)
	{
		$rdist = rand(0.6);
		if ($rdist >= 0.5)
		{
			$rdist = 0.5;
			$chromosome_len = 0;
		}
		else
		{
			$chromosome_len++;
			if ($chromosome_len > $max_chromosome_len)
			{
				$rdist = 0.5;
				$chromosome_len = 0;
			}
		}
		
		$recomb_dists[$i] = $rdist;
	}
		
	$self->makeGenome (Nloci => $nloci, RecombinationDists => \@recomb_dists,
				Name => $g_name,
				AlleleSet => $a_set);
}


# load the named genome from the persistant database.
sub loadGenome
{
	my $self = shift();
	my %params = @_;

	die "Sorry, database functionality has not been implemented yet!\n";
	#$self->makeGenome(GenomeName => $params{GenomeName});	
}


sub getName
{
	my $self = shift();
	return $self->{Name};	
}

sub setName
{
	my $self = shift;
	my %params = @_;
	
	croak "'Name' missing" unless exists $params{Name};
	
	$self->{Name} = $params{Name};
}



sub getAlleleSet
{
	my $self = shift;

	return $self->{_AlleleSet};
}

# check and see if the specified Allele exists in the AlleleSet of this Genome
sub validAllele
{
	my $self = shift;
	my %params = @_;
	
	croak "'Allele' missing" unless exists $params{Allele};

	my $allele_set = $self->{_AlleleSet};
	
	if (index ($allele_set, $params{Allele}) > -1)
	{
		return 1;
	}
		
	return 0;
}


# Return the number of loci in the Genome.
sub getNloci
{
	my $self = shift();
	return $self->{Nloci};	
}


# return a reference to the recombination distances array
sub getRecombinationDistancesRef
{
	my $self = shift();
	return $self->{_RecombinationDists};	
}

# return a copy of the recombination distances array
sub getRecombinationDistances
{
	my $self = shift;

	my ($dist, @rdists);
	my $recomb_dists = $self->{_RecombinationDists};
	foreach $dist (@$recomb_dists)
	{
		push @rdists, $dist;	
	}
	
	return @rdists;
}


##### routines for handling traits in the genome

# add a new trait definition to the Genome.
sub addTrait
{
	my $self = shift;
	my %params = @_;
	
	croak "'Name' missing" unless exists $params{Name};
	croak "'TraitMean' missing" unless exists $params{TraitMean};
	croak "'TraitVariance' missing" unless exists $params{TraitVariance};
	croak "'TraitAlleleEffects' missing" unless exists $params{TraitAlleleEffects};
	
	croak "'AlleleSet' undefined" unless defined $self->{_AlleleSet};
	
	my $tae = $params{TraitAlleleEffects};
	if (not ref $tae)
	{
		croak "'TraitAlleleEffects' must be passed by reference";	
	}
	
	# check and see if the trait already exists in the table
	my $trait_name = $params{Name};
	my $traits = $self->{_Traits};
	
	if ( exists $$traits{$trait_name} )
	{
		croak "A trait with the name '$trait_name' already exists";	
	}
	
	my $trait = new Gnomish::Trait;

	$trait->makeTrait (Name => $params{Name},
				TraitMean => $params{TraitMean}, 
				TraitVariance => $params{TraitVariance},
				TraitAlleleEffects => $tae,
				AlleleSet => $self->{_AlleleSet} );
				
	$$traits{$trait_name} = $trait;
}

sub addMarker
{
	my $self = shift;
	my %params = @_;

	croak "'Name' missing" unless exists $params{Name};
	croak "'MarkerCoding' missing" unless exists $params{MarkerCoding};
	
	# check and see if the marker already exists in the table
	my $marker_name = $params{Name};
	my $markers = $self->{_Markers};
	
	if ( exists $$markers{$marker_name} )
	{
		croak "A maker with the name '$marker_name' already exists";	
	}
	
	my $marker = new Gnomish::Marker;
	$marker->makeMarker (Name => $params{Name}, 
		MarkerCoding => $params{MarkerCoding},
		AlleleSet => $self->{_AlleleSet});
	
	$$markers{$marker_name} = $marker;
}

# return a trait from the genome
sub getTrait
{
	my $self = shift;
	my %params = @_;

	croak "'TraitName' missing" unless exists $params{TraitName};
	my $trait_name = $params{TraitName};
	
	my $traits = $self->{_Traits};
	
	if ( !exists $$traits{$trait_name} )
	{
		croak "Specified trait doesn't exist";	
	}
	
	return $$traits{$trait_name};
}

# return a marker from the genome
sub getMarker
{
	my $self = shift;
	my %params = @_;

	croak "'MarkerName' missing" unless exists $params{MarkerName};
	my $marker_name = $params{MarkerName};

	my $markers = $self->{_Markers};

	if ( !exists $$markers{$marker_name} )
	{
		croak "Specified marker doesn't exist";	
	}
	
	return $$markers{$marker_name};
}

# return the reference to the specified trait.
sub getTraitRef
{
	my $self = shift;
	my %params = @_;

	croak "'Name' missing" unless exists $params{Name};
	my $trait_name = $params{Name};
	
	my $traits = $self->{_Traits};
	
	if ( !exists $$traits{$trait_name} )
	{
		croak "Specified trait doesn't exist";	
	}
	
	return \$$traits{$trait_name};
}

sub getMarkerRef
{
	my $self = shift;

	my %params = @_;

	croak "'MarkerName' missing" unless exists $params{MarkerName};
	my $marker_name = $params{MarkerName};

	my $markers = $self->{_Markers};

	if ( !exists $$markers{$marker_name} )
	{
		croak "Specified marker doesn't exist";	
	}
	
	return \$$markers{$marker_name};
}

# return the number of traits that exists in the genome
sub getTraitNames
{
	my $self = shift;

	my $traits = $self->{_Traits};

	return keys %$traits;	
}

sub getMarkerNames
{
	my $self = shift;

	my $markers = $self->{_Markers};
	
	return keys %$markers;	
}


##### misc methods


# Display a representation of the genome object in ASCII
sub drawInAscii
{
	my $self = shift();

	my ($i, $j, $dist);
	my ($rdist);
	my ($nloci, $genome_name, $recomb_dists);
	
	$genome_name  = $self->{Name};
	$nloci        = $self->{Nloci};
	$recomb_dists = $self->{_RecombinationDists};
	
	print "Genome Name: $genome_name\n";
	print "Nloci      : $nloci\n";
	
	if ($nloci < 1)
	{
		return;	
	}
	
	print "Recombination distances: ";
	foreach $rdist (@$recomb_dists)
	{
		print "$rdist, ";
	}
	print "\n";
	
	my $new_gene = 0;
	print "|"; # get things started
	for ($i = 0; $i < ($nloci - 1); $i++)
	{
		if ($$recomb_dists[$i] == 0.5)
		{
			$new_gene = 1;
			print "\n|";
			next;
		}
		
		# -- Draw the loci --
		# first we compute the actual distance by converting it into an integer
		$dist = int ($$recomb_dists[$i] * 100.0 + 0.5);
		print "$dist";
		for ($j = 0; $j < $dist; $j++)
		{
			print "-";	
		}
		print "|";
	}
	print "\n";
}


1;
