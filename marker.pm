# Marker object.  Markers and Traits are very similar objects, however for the sake of
# clarity, we treat them as seperate entities.
#
# File: Marker.pm
# Author: Hai Pham, hpham42@gmail.com
# Date: Jan 31, 2003
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



#Markers and Traits are very similar object, however for the sake of clarity, we treat
#them as seperate entities.
package Gnomish::Marker;

use strict;
use Carp;

use Gnomish::Stdlib;

# Name => $name, optional name of the marker
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


# Create a new marker list
# Parameters:
#  Name => $name
#  MarkerCoding => %marker_coding_table
#  AlleleSet => $allele_set
#
# Note: 'AlleleSet' is a comma seperated string that specifies the set of valid allele labels
# that belong to the Genome.
sub makeMarker
{
	my $self = shift;
	my %params = @_;

	croak "'Name' missing" unless exists $params{Name};
	croak "'MarkerCoding' missing" unless exists $params{MarkerCoding};
	croak "'AlleleSet' missing" unless exists $params{AlleleSet}; 
	
	$self->{Name} = $params{Name};
	
	my $mc = $params{MarkerCoding};
	if ( not ref $mc )
	{
		croak "'MarkerCoding' must be passed by reference";	
	}
	$self->_setMarkerCoding ($mc);
	
	$self->_setAlleleSet ($params{AlleleSet});
}


# The expected reference should be a hash of hashes in the form of
# $table_name{$locus}{$allele_state} = $marker_code 
sub _setMarkerCoding
{
	my $self = shift;
	my $mc   = shift;
	
	# make a local copy of the hash
	my (%codes_table, $locus);
	foreach $locus (keys %$mc)
	{
		my %marker_codes = %{ $$mc{$locus} };
		$codes_table{$locus} = \%marker_codes;
	}
	
	$self->{_MarkerCoding} = \%codes_table;
}

# this is the same code from the Traits.pm object
sub _setAlleleSet
{
        my $self = shift;
        my $a_set = shift;
        
        # make sure that there are no spaces in the allele set string
        $a_set =~ s/[\s\n\r]//g;
        
        my @allele_set = split ',', $a_set;
        
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

# return the marker coding for a parcular allele state
sub getMarkerCoding
{
	my $self = shift;
	my %params = @_;

        croak "'Locus' missing" unless exists $params{Locus};
        croak "'AlleleState' missing" unless exists $params{AlleleState};
        
        my $locus = $params{Locus};
        my $effects_table = $self->{_MarkerCoding};
        my $state_ord = $self->_getAlleleStateOrd ($params{AlleleState});

        return $effects_table->{$locus}{$state_ord};
}

# return a hash of the marker codings for a given locus
sub getMarkerCodingForLocus
{
	my $self = shift;
	my %params = @_;
	
	croak "'Locus' missing" unless exists $params{Locus};
	
	my $locus = $params{Locus};
	my $coding_table = $self->{_MarkerCoding};
	my %marker_codings = %{ $$coding_table{$locus} };
	
	return %marker_codings;
}

# return the list of loci that belong to the marker
sub getMarkerLoci
{
	my $self = shift;
	
	my $mt = $self->{_MarkerCoding};
	my @loci = keys %$mt;
	
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

sub printMarkerCoding
{
	my $self = shift;

	my $mc = $self->{_MarkerCoding};
	
	my ($locus, $allele);
	foreach $locus (sort keys %$mc)
	{
		print "locus: $locus ";
		my %coding = %{$$mc{$locus}};
		foreach $allele (sort keys %coding)
		{
			print "$allele: ", $coding{$allele}, " ";	
		}
		print "\n";
	}
}

1;

