# Database interfaces for the various Genomish classes
#
# File: DB.pm
# Author: Hai Pham, hpham42@gmail.com
# Last Modified: Feb 26, 2003
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

# NOTE:
# The routines in this module are still in flux, and are inconsistent
# with the interface conventions established by the other modules.
# For example, most of the subroutines here do not have labled
# parameters (eg. GenomeName => $genome_name).
# I'm also not quite sure about about sorts of functionalilty
# people would find useful in the database layer.... so this
# stuff is highly experimental.


package Gnomish::DB;

use strict;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(DB_setParams DB_init 
			DB_saveGenome DB_savePopulation
			DB_loadGenome DB_loadPopulation
			DB_deleteGenome DB_deletePopulation
			DB_renameMarker DB_renameTrait
			DB_renameGenome DB_renamePopulation
			DB_listMarkerNames DB_listTraitNames
			DB_listPopulationNames DB_listGenomeNames);

use DBI;

my ($DB_Type, $DB_Name, $DB_User, $DB_Passwd);

# Keep track of the primary key in each table in the database.  This is
# used by the subroutine '_check_existence'
my %Table_Primary_Keys = (
	"Genome" => "GenomeID",
	"Genotype" => "GenotypeID",
	"Population" => "PopulationID",
	"Traits" => "TraitID",
	"Markers" => "MarkerID"
);



# Set the database parameters.
# Parameters:
# DB_Type: database type ('ODBC' or 'mysql')
# DB_Name: name of the database, or in the case of ODBC, the name of the database handle.
# DB_User: user name to use for database authentication
# DB_Passwd: database password
#
# NOTE: This must be called before any access to the database can be made!
sub DB_setParams (%)
{
	my %params = @_;
	
	croak "'DB_Type' missing" unless $params{DB_Type};
	croak "'DB_Name' missing" unless $params{DB_Name};
	croak "'DB_User' missing" unless $params{DB_User};
	croak "'DB_Passwd' missing" unless $params{DB_Passwd};
	
	$DB_Type = $params{DB_Type};
	$DB_Name = $params{DB_Name};
	$DB_User = $params{DB_User};
	$DB_Passwd = $params{DB_Passwd};
}


# Initialize the database, and return database handle.
# This doesn't necessarily need to be called, but for the sake of efficiency,
# it should be called after a call to 'DB_setParams'.
#
# For example, a typical script with database calls can look like:
#
# DB_setParams (...);
# my $genome = DB_loadGenome ( .... );
# my $population = DB_loadPopulation ( .... );
#
# However the above code gets in efficient if there are lots of
# calls to database code because it forces the code to establish a new
# connection to the database each time there is a database call.
# Code with lots of database access would be more efficent coded like this:
#
# DB_serParams ( ... );
# my $dbh = DB_init();
# my $genome = DB_loadGenome ( ..., $dbh);
# my $population = DB_loadPopulation ( ...., $dbh);
#
# This way, a single database handle (connection) is reused over and over,
# avoiding the need to re-establish a new connection each time a database
# access is made.
sub DB_init ()
{
	my $dbh;
	
	# ACCESS defaults to a String read len that is too short
	# for the AlleleStrings and Recombination Distances arrays
	if ($DB_Type eq "ODBC")
	{
		$dbh = DBI->connect ("DBI:$DB_Type:$DB_Name", $DB_User, $DB_Passwd,
					{LongReadLen => 60000})
			or croak "Couldn't connect to database: " . DBI->errstr;

	}	
	else
	{
		$dbh = DBI->connect ("DBI:$DB_Type:$DB_Name", $DB_User, $DB_Passwd)
			or croak "Couldn't connect to database: " . DBI->errstr;
	}

	return $dbh;		
}


# Internal subroutine to retrieve the row ID of the last insert statement.
# Due to differences in the way that ODBC and mysql works, both the database handle
# and the statment handle (from the DBI driver prepare command) has to be
# passed to this subroutine.
sub _getInsertID ($$)
{
	my $sth = shift;
	my $dbh = shift;
	
	my $id;
        if ($DB_Type eq "mysql")
        {
                $id = $sth->{mysql_insertid};    
        }
        elsif ($DB_Type eq "ODBC")
        {
                $id = $dbh->selectrow_array ("select \@\@IDENTITY");
        }
        else
        {
                croak "Database not supported yet";     
        }

	return $id;
}

# check and see if a particular entry already exists in the database
# column names & values should be specified as
# $params{ColumnName} = 'ColumnValue'; 
sub check_existence ($$%)
{
	my $dbh     = shift;
	my $table   = shift;
	my %columns = @_;

	
	die "The primary key for the table '$table' is undefined" 
		unless exists $Table_Primary_Keys{$table};
	my $tbl_primary_key = $Table_Primary_Keys{$table};
	
	my $q = "select $tbl_primary_key from $table where ";
	my ($col_name, $col_data);
	foreach $col_name (keys %columns)
	{
		$col_data = $columns{$col_name};
		$q .= "$col_name='$col_data',";
	}
	chop $q; # remove the trailing ,
	
	my $sth = $dbh->prepare ($q) or croak $dbh->errstr;
	$sth->execute() or croak $sth->errstr;
	
	my $href = $sth->fetchrow_hashref();
	if (defined $href)
	{
		my $tbl_primary_key = $Table_Primary_Keys{$table};
		my $key_val = $href->{$tbl_primary_key};
		
		return $key_val;
	}
	
	return 0;
}

# check to see if a particular entry already exists in the database.
# if it exists, delete it.  Return 1 if a record was deleted, otherwise return 0
sub check_and_delete ($$%)
{
	my $dbh     = shift;
	my $table   = shift;
	my %columns = @_;
	
	my $key_val = check_existence ($dbh, $table, %columns);
	if ($key_val)
	{
		my $tbl_primary_key = $Table_Primary_Keys{$table};
		$dbh->do ("delete from $table where $tbl_primary_key=$key_val")
			or croak $dbh->errstr;
			
		return 1;
	}
	
	return 0;
}

# get the GenomeID of the named genome.
sub _getGenomeID ($$)
{
	my $genome_name = shift;
	my $dbh = shift;

	my $genome_id = $dbh->selectrow_array ("select GenomeID 
				from Genome where GenomeName='$genome_name'")
			or croak $dbh->errstr;
			
	return $genome_id;
}

# save the specified Genome to the DB
sub DB_saveGenome
{
	my $genome    = shift;
	my $dbh       = shift;

	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	if (ref $genome eq "Gnomish::Genome")
	{
		croak "Genome object must be passed by reference";	
	}
	
	my $genome_name = $$genome->getName();
	my $allele_set  = $$genome->getAlleleSet();
	my $nloci       = $$genome->getNloci();
	my @r_dists     = $$genome->getRecombinationDistances();
	my $recombination_dists = join ',', @r_dists;

	if (not defined $genome_name)
	{
		croak "Genome name must be defined to be saved into database";
	}
	
	my $gid = check_existence ($dbh, "Genome", GenomeName => $genome_name);
	if ($gid)
	{
		croak "Genome '$genome_name' already exists in database";	
	}
	
	# first save the general genome characteristics
	my $sth = $dbh->prepare ("insert into Genome
					(GenomeName,Nloci,RecombinationDists,AlleleSet)
					values
					('$genome_name',$nloci,'$recombination_dists','$allele_set')")
			or croak $dbh->errstr;
	
	#my $sth = $dbh->prepare ("insert into Genome set
	#			GenomeName='$genome_name',
	#			Nloci='$nloci',
	#			RecombinationDists='$recombination_dists',
	#			AlleleSet='$allele_set'")
	#		or croak $dbh->errstr;
	$sth->execute() or croak $sth->errstr;
	
	my $genome_id = _getInsertID ($sth,$dbh);
	
	# save Traits to the database
	my @trait_names = $$genome->getTraitNames();
	my ($tname, $trait_id);
	foreach $tname (@trait_names)
	{
		my $trait = $$genome->getTraitRef (Name => $tname);
		
		$trait_id = DB_saveTrait ($trait, $dbh);
		
		# update the pivot table to link traits & genomes
		#$dbh->do ("insert into GenomeTraits set
		#			GenomeID='$genome_id',
		#			TraitID='$trait_id'")
		#	or croak $dbh->errstr;
		$dbh->do ("insert into GenomeTraits
				(GenomeID,TraitID) values
				($genome_id,$trait_id)")
			or croak $dbh->errstr;
	}
	
	# save Markers to the database
	my @marker_names = $$genome->getMarkerNames();
	my ($marker_name, $marker_id);
	foreach $marker_name (@marker_names)
	{
		my $marker = $$genome->getMarkerRef (MarkerName => $marker_name);
		
		$marker_id = DB_saveMarker ($marker, $dbh);
		
		#$dbh->do ("insert into GenomeMarkers set
		#		GenomeID='$genome_id',
		#		MarkerID='$marker_id'")
		#	or croak $dbh->errstr;
		$dbh->do ("insert into GenomeMarkers (GenomeID,MarkerID) values
				($genome_id,$marker_id)")
			or croak $dbh->errstr;
	}
}

# save the marker to the DB, return the MarkerID of the saved marker
sub DB_saveMarker
{
	my $marker = shift;
	my $dbh = shift;

	if (ref $marker eq "Gnomish::Marker")
	{
		croak "Marker object must be passed by reference";
	}
	
	my $marker_name = $$marker->getName();
	
	my $id = check_existence ($dbh, "Markers", MarkerName => $marker_name);
	if ($id)
	{
		croak "The marker '$marker_name' already exists in the database";
	}
	
	#my $sth = $dbh->prepare ("insert into Markers set
	#				MarkerName='$marker_name'")
	#		or croak $dbh->errstr;
	my $sth = $dbh->prepare ("insert into Markers (MarkerName) values
					('$marker_name')")
			or croak $dbh->errstr;
	$sth->execute() or croak $sth->errstr;
	
	my $marker_id = _getInsertID ($sth,$dbh);
	
	my @loci = $$marker->getMarkerLoci();
	my ($locus);
	foreach $locus (@loci)
	{
		my %locus_coding = $$marker->getMarkerCodingForLocus( Locus => $locus );
		my ($allele_state, $marker_label);
		foreach $allele_state (keys %locus_coding)
		{
			$marker_label = $locus_coding{$allele_state};
			#$dbh->do ("insert into MarkerAlleleCodes set
			#		MarkerID='$marker_id',
			#		Locus='$locus',
			#		AlleleState='$allele_state',
			#		AlleleCode='$marker_label'")
			#	or croak $dbh->errstr;
			$dbh->do ("insert into MarkerAlleleCodes
					(MarkerID,Locus,AlleleState,AlleleCode) values
					($marker_id,$locus,'$allele_state','$marker_label')")
				or croak $dbh->errstr;
		}
	}

	return $marker_id;
}

# save the specified Trait to the DB, return the TraitID
sub DB_saveTrait
{
	# reference to trait object
	my $trait = shift;
	my $overwrite = shift;
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	if (ref $trait eq "Gnomish::Trait")
	{
		croak "Trait object must be passed by reference";	
	}
	
	my $trait_name = $$trait->getName();
	my $trait_variance = $$trait->getTraitVariance();
	my $trait_mean = $$trait->getTraitMean();
	
	my $tid = check_existence ($dbh, "Traits", TraitName => $trait_name);
	if ( $tid )
	{
		croak "Trait '$trait_name' already exists in database";	
	}

	#my $sth = $dbh->prepare ("insert into Traits set
	#		TraitName='$trait_name', 
	#		TraitVariance='$trait_variance',
	#		TraitMean='$trait_mean'") 
	#		or croak $dbh->errstr;
	my $sth = $dbh->prepare ("insert into Traits
					(TraitName,TraitVariance,TraitMean) values
					('$trait_name',$trait_variance,$trait_mean)")
			or croak $dbh->errstr;
	$sth->execute() or croak $sth->errstr;
	
	my $trait_id = _getInsertID($sth,$dbh);
	
	# now save the allele effects			
	my @active_loci = $$trait->getActiveLoci();
	my (%allele_effects, $locus);
	foreach $locus (@active_loci)
	{
		my %locus_effects = $$trait->getLocusEffects(Locus => $locus);
		my ($allele_state, $effect_val);
		foreach $allele_state (keys %locus_effects)
		{
			$effect_val = $locus_effects{$allele_state};
			#$dbh->do ("insert into TraitAlleleEffects set
			#		TraitID='$trait_id',
			#		Locus='$locus',
			#		AlleleState='$allele_state', 
			#		EffectValue='$effect_val'")
			#		or croak $dbh->errstr;
			$dbh->do ("insert into TraitAlleleEffects
					(TraitID,Locus,AlleleState,EffectValue) values
					($trait_id,$locus,'$allele_state',$effect_val)")
				or croak $dbh->errstr;
		}
			
	}
	
	return $trait_id;	
}


# save a Population to the database
sub DB_savePopulation
{
	my $pop = shift;
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();	
	}
	
	if (ref $pop eq "Gnomish::Population")
	{
		croak "Population object must be passed by reference";	
	}
	
	my $pop_name = $$pop->getName();
	if ( (length $pop_name) == 0 )
	{
		croak "The Population must have a name to be saved in the database"; 	
	}
	
	my $gid = check_existence ($dbh, "Population", PopulationName => $pop_name);
        if ($gid)
        {
                croak "Genome '$pop_name' already exists in database";       
        }
	
	# get the genomeID of the genome that belongs to the population
	my $nmembers = $$pop->nMembers();
	
	my ($genome_id, $m);
	if ($nmembers < 1)
	{
		croak "The population has no genotypes defined!";
	}
	
	$m = $$pop->getMemberRef (Member => 0);
	my $genome_name = $$m->getGenomeName();
	$genome_id = _getGenomeID ($genome_name, $dbh); 
	
	# save the population name & retrieve it's ID on the DB
	#my $sth = $dbh->prepare ("insert into Population set
	#				PopulationName='$pop_name',
	#				GenomeID=$genome_id")
	#		or croak $dbh->errstr;
	my $sth = $dbh->prepare ("insert into Population (PopulationName,GenomeID)
					values ('$pop_name',$genome_id)")
			or croak $dbh->errstr;
	$sth->execute() or croak $sth->errstr;
	
	my $pop_id = _getInsertID ($sth, $dbh);
	
	# save the genotypes
	#$sth = $dbh->prepare ("insert into PopulationGenotypes set
	#			PopulationID=?,
	#			GenotypeId=?")
	#	or croak $dbh->errstr;
	$sth = $dbh->prepare ("insert into PopulationGenotypes
				(PopulationID,GenotypeID) values (?,?)")
		or croak $dbh->errstr;
	for (my $n = 0; $n < $nmembers; $n++)
	{
		$m = $$pop->getMemberRef(Member => $n);
		my $genotype_id = DB_saveGenotype ($m, $genome_id, $dbh);
		
		$sth->execute ($pop_id, $genotype_id); 
	}
}

# We should make this into a more generalized interface!!
# should be able to cope with no GenomeID, and no DBH!
sub DB_saveGenotype ($$$)
{
	my $genotype = shift;
	my $genome_id = shift;
	my $dbh = shift;

	if (not defined $dbh)
	{
		$dbh = DB_init();	
	}
	
	if (ref $genotype eq "Gnomish::Genotype")
	{
		croak "Genotype object must be passed by reference";
	}
	
	my $genotype_name = $$genotype->getName();
	if (not defined $genotype_name)
	{
		$genotype_name = 'NULL';
	}
	
	my $allele_string0 = $$genotype->getGamete (Gamete => 0);
	my $allele_string1 = $$genotype->getGamete (Gamete => 1);
	
	# we don't check uniqueness of genotypes, because genotypes do not
	# have to be unique...
	#my $sth = $dbh->prepare ("insert into Genotype set
	#				GenomeID='$genome_id',
	#				GenotypeName='$genotype_name',
	#				AlleleString0='$allele_string0',
	#				AlleleString1='$allele_string1'")
	#	or croak $dbh->errstr;
	my $sth = $dbh->prepare ("insert into Genotype
					(GenomeID,GenotypeName,AlleleString0,AlleleString1)
					values
					($genome_id,'$genotype_name','$allele_string0','$allele_string1')")
			or croak $dbh->errstr;
	$sth->execute or croak $sth->errstr;

	my $genotype_id = _getInsertID ($sth, $dbh);
	
	### now save genotypic & phenotypic values ###
	
	# first we need to get the phenotypic & genotypic value hashes
	my @trait_names = $$genotype->getTraitNames();
	my $trait_name;
	
	my ($pv, $gv);
	#$sth = $dbh->prepare ("insert into GenotypeValues set
	#				GenotypeID=?,
	#				TraitID=?,
	#				GenotypicValue=?,
	#				PhenotypicValue=?")
	#		or croak $dbh->errstr;
	$sth = $dbh->prepare ("insert into GenotypeValues
				(GenotypeID,TraitID,GenotypicValue,PhenotypicValue) 
				values (?,?,?,?)")
		or croak $dbh->errstr;
	foreach $trait_name (@trait_names)
	{
		$pv = $$genotype->getPhenotypicValueForTrait (TraitName => "$trait_name");
		$gv = $$genotype->getGenotypicValueForTrait (TraitName => "$trait_name");
		
		# insert them into the GenotypeValues table, but first
		# we need to get the TraitID from the database
		my $trait_id = $dbh->selectrow_array ("select TraitID
					from Traits where TraitName='$trait_name'")
				or croak $dbh->errstr;
		
		$sth->execute ($genotype_id, $trait_id, $gv, $pv)
			or croak $sth->errstr;
	}
	
	return $genotype_id;
}

####################### LOAD ROUTINES ########################
sub DB_loadGenome ($$)
{
	my $genome_name = shift;
	my $dbh    = shift;

	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	my $genome_id = check_existence ($dbh, "Genome", GenomeName => $genome_name);
	if ($genome_id == 0)
	{
		croak "The Genome '$genome_name' doesn't seem to exist in the database";			
	}
	
	# first load the genome's various attributes available in the immediate table
	my $sth = $dbh->prepare ("select * from Genome where GenomeID=$genome_id")
				or croak $dbh->errstr;
	$sth->execute or croak $sth->errstr;
	
	my $g = $sth->fetchrow_hashref();
	my $nloci      = $g->{Nloci};
	my $allele_set = $g->{AlleleSet};
	my $recombination_dists = $g->{RecombinationDists};

	# create a new genome object
	my $genome = new Gnomish::Genome;
	$genome->makeGenome ( Name => $genome_name,
				Nloci => $nloci,
				RecombinationDists => $recombination_dists,
				AlleleSet => $allele_set );
		
	# Now load up the traits.
	# First, get the ID of traits associated with this Genome
	$sth = $dbh->prepare ("select TraitID from GenomeTraits where GenomeID=$genome_id")
			or croak $dbh->errstr;
	$sth->execute or croak $sth->errstr;
	my (@trait_ids, $t);
	while ( $t = $sth->fetchrow_arrayref() )
	{
		push @trait_ids, $t->[0];	
	}
	
	my $t_id;
	foreach $t_id (@trait_ids)
	{
		$sth = $dbh->prepare ("select * from Traits where TraitID=$t_id")
				or croak $dbh->errstr;
		$sth->execute() or croak $sth->errstr;

		$t = $sth->fetchrow_hashref();
		my $TraitName = $t->{TraitName};
		my $TraitMean = $t->{TraitMean};
		my $TraitVariance = $t->{TraitVariance};

		my %effects_table = _loadTraitAlleleEffectsTable ($dbh, $t_id);

		$genome->addTrait (Name => $TraitName,
					TraitMean => $TraitMean,
					TraitVariance => $TraitVariance,
					TraitAlleleEffects => \%effects_table);
	}
	
	# now load up the markers...
	$sth = $dbh->prepare ("select MarkerID from GenomeMarkers where GenomeID=$genome_id")
		or croak $dbh->errstr;
	$sth->execute() or croak $sth->errstr;
	my (@marker_ids, $m);
	while ($m = $sth->fetchrow_arrayref())
	{
		push @marker_ids, $m->[0];
	}
	
	my $m_id;
	foreach $m_id (@marker_ids)
	{
		$sth = $dbh->prepare ("select * from Markers where MarkerID=$m_id")
			or croak $dbh->errstr;
		$sth->execute() or croak $sth->errstr;
		
		$m = $sth->fetchrow_hashref();
		my $marker_name = $m->{MarkerName};
		
		my %marker_coding = _loadMarkerCodingTable ($dbh, $m_id);
		
		$genome->addMarker (Name => $marker_name,
					MarkerCoding => \%marker_coding);
	}
	
	return $genome;
}

sub _loadTraitAlleleEffectsTable ($$)
{
	my ($dbh, $trait_id) = @_;
	
	my $sth = $dbh->prepare ("select * from TraitAlleleEffects where TraitID=$trait_id")
			or croak $dbh->errstr;
	$sth->execute() or croak $sth->errstr;
	
	my ($ae, $locus, $allele_state, $effect_value, %effects_table);
	while ( $ae = $sth->fetchrow_hashref() )
	{
		$locus        = $ae->{Locus};
		$allele_state = $ae->{AlleleState};
		$effect_value = $ae->{EffectValue};
		
		$effects_table{$locus}{$allele_state} = $effect_value;
	}
	
	return %effects_table;
}

sub _loadMarkerCodingTable ($$)
{
	my ($dbh, $marker_id) = @_;

	my $sth = $dbh->prepare ("select * from MarkerAlleleCodes where MarkerID=$marker_id")
			or croak $dbh->errstr;
	$sth->execute() or croak $sth->errstr;
	
	my ($locus, $allele_state, $allele_label, %coding_table, $mc);
	while ( $mc = $sth->fetchrow_hashref() )
	{
		$locus = $mc->{Locus};
		$allele_state = $mc->{AlleleState};
		$allele_label = $mc->{AlleleCode};
		
		$coding_table{$locus}{$allele_state} = $allele_label;
	}
	
	return %coding_table;
}


sub DB_loadPopulation
{
	my ($pop_name, $dbh) = @_;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	my ($pop_id, $genome_id) = $dbh->selectrow_array ("select PopulationID,GenomeID
					from Population where PopulationName='$pop_name'")
					or croak $dbh->errstr;
	
	my $population = new Gnomish::Population (Name => $pop_name);
	
	# load up the genome object, but first get the genome name...
	my $genome_name = $dbh->selectrow_array ("select GenomeName from
				Genome where GenomeID=$genome_id")
				or croak $dbh->errstr;
	
	my $genome = DB_loadGenome ($genome_name, $dbh);
	
	# now get the list of GenotypeID's the belong to the population
	my $sth = $dbh->prepare ("select GenotypeID from PopulationGenotypes
					where PopulationID=$pop_id")
			or croak $dbh->errstr;
	$sth->execute() or croak $sth->errstr;
	my $array_ref = $sth->fetchall_arrayref();
	
	# load up the genotypes
	foreach my $row (@$array_ref)
	{
		my ($genotype_id) = @$row;
		
		my $genotype = _loadGenotype($genotype_id, \$genome, $dbh);
		
		$population->addMember (Genotype => \$genotype);
	}
	
	return $population;
}

sub _loadGenotype ($$$)
{
	my ($genotype_id, $genome, $dbh) = @_;

	# load the genotype from the database
	my ($genome_id, $genotype_name, $allele_str0, $allele_str1) =
		$dbh->selectrow_array ("select GenomeID,GenotypeName,AlleleString0,AlleleString1
			from Genotype where GenotypeID=$genotype_id")
		or croak $dbh->errstr;
	
	# load up the phenotypic & genotypic values
	my $sth = $dbh->prepare ("select * from GenotypeValues where GenotypeID=$genotype_id")
		or croak $dbh->errstr;
	$sth->execute or croak $sth->errstr;

	my (%trait_names, %genotypic_values, %phenotypic_values);	
	my $array_ref = $sth->fetchall_arrayref();
	foreach my $row (@$array_ref)
	{
		my ($g_id, $trait_id, $genotypic_value, $phenotypic_value) = @$row;
		
		if (not exists $trait_names{$trait_id})
		{
			my ($t_name) = $dbh->selectrow_array 
						("select TraitName from Traits where TraitID=$trait_id")
						or croak $dbh->errstr;
			$trait_names{$trait_id} = $t_name;
		}
		
		my $t_name = $trait_names{$trait_id};
		
		$genotypic_values{$t_name} = $genotypic_value;
		$phenotypic_values{$t_name} = $phenotypic_value;
	}
	
	my $genotype = new Gnomish::Genotype;
		
	$genotype->makeGenotype (Genome => $genome,
					AlleleString0 => $allele_str0,
					AlleleString1 => $allele_str1,
					GenotypicValues => \%genotypic_values,
					PhenotypicValues => \%phenotypic_values);
					
	return $genotype;
}


############## Utility functions ############

sub _listNames ($$$)
{
	my $dbh = shift;
	my $tbl_name = shift;
	my $col_name = shift;
	
	my $sth = $dbh->prepare ("select $col_name from $tbl_name")
			or croak $dbh->errstr;
	$sth->execute or croak $sth->errstr;
	
	my (@names, $n);
	while ($n = $sth->fetchrow_arrayref())
	{
		push @names, $n->[0];
	}
	
	return @names;
}

sub DB_listGenomeNames
{
	my $dbh = shift;

	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	my @genome_names = _listNames ($dbh, 'Genome', 'GenomeName');
	
	return @genome_names;
}

sub DB_listPopulationNames
{
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	my @pop_names = _listNames ($dbh, 'Population', 'PopulationName');
	
	return @pop_names;
}

sub DB_listTraitNames
{
	my $dbh = shift;

	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	my @trait_names = _listNames ($dbh, 'Traits', 'TraitName'); 
}

sub DB_listMarkerNames
{
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	my @marker_names = _listnames ($dbh, 'Markers', 'MarkerName');
	
	return @marker_names;
}

# Do the actual renaming
# Don't really like this (particularly the %elements table), may want to redesign...
sub _renameElement ($$$$)
{
	my ($dbh, $element_type, $old_name, $new_name) = @_;
	
	my %elements = ( 'Genome'     => 'Genome:GenomeName:GenomeID',
			 'Population' => 'Population:PopulationName:PopulationID',
			 'Trait'      => 'Traits:TraitName:TraitID',
			 'Marker'     => 'Markers:MarkerName:MarkerID');
			 
	die unless exists $elements{$element_type};
	
	my ($table, $col_name, $col_id) = split ':', $elements{$element_type};
	
	# first we have to get the row id of the record
	my $id = check_existence ($dbh, $table, "$col_name" => "$old_name");
	
	if ($id == 0)
	{
		croak "'$old_name' does not seem to exist in the $table table";
	}
	
	$dbh->do ("update $table set $col_name = '$new_name' where $col_id = $id")
		or croak $dbh->errstr;
		
	return;
}

sub DB_renamePopulation
{
	my $old_name = shift;
	my $new_name = shift;
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	_renameElement ($dbh, 'Population', $old_name, $new_name);
}

sub DB_renameGenome
{
	my $old_name = shift;
	my $new_name = shift;
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	_renameElement ($dbh, 'Genome', $old_name, $new_name);
}

sub DB_renameTrait
{
	my $old_name = shift;
	my $new_name = shift;
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	_renameElement ($dbh, 'Trait', $old_name, $new_name);
}

sub DB_renameMarker
{
	my $old_name = shift;
	my $new_name = shift;
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	_renameElement ($dbh, 'Marker', $old_name, $new_name);
}


sub DB_deletePopulation
{
	my $pop_name = shift;
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	my ($pop_id) = check_existence ($dbh, "Population", PopulationName => $pop_name);
	if ($pop_id == 0)
	{
		warn "The Population '$pop_name' doesn't exist in the database";
		return;
	}
	
	# first, we have to retrieve all the references to the genotypes from the db
	my $sth = $dbh->prepare ("select GenotypeID from PopulationGenotypes
					where PopulationID=$pop_id")
			or croak $dbh->errstr;
	$sth->execute or croak $sth->errstr;
	
	# now delete the genotypes, one at a time
	my $array_ref = $sth->fetchall_arrayref();
	foreach my $row (@$array_ref)
	{
		my ($genotype_id) = @$row;
		
		_deleteGenotype ($genotype_id, $dbh);
	}
	
	# finally, delete the population entry in the db
	$dbh->do ("delete from Population where PopulationID=$pop_id")
		or croak $dbh->errstr;
}

sub _deleteGenotype ($$)
{
	my $genotype_id = shift;
	my $dbh = shift;
	
	$dbh->do ("delete from GenotypeValues where GenotypeID=$genotype_id")
		or croak $dbh->errstr;
	$dbh->do ("delete from Genotype where GenotypeID=$genotype_id")
		or croak $dbh->errstr;
	$dbh->do ("delete from PopulationGenotypes where GenotypeID=$genotype_id")
		or croak $dbh->errstr;
}

sub DB_deleteGenome
{
	my $genome_name = shift;
	my $dbh = shift;
	
	if (not defined $dbh)
	{
		$dbh = DB_init();
	}
	
	my ($genome_id) = check_existence ($dbh, "Genome", GenomeName => $genome_name);
	if ($genome_id == 0)
	{
		warn "The Genome '$genome_name' doesn't exist in the database";
		return;
	}
	
	# first check and see if there are still any populations
	# that use this genome
	my ($pid) = check_existence ($dbh, "Population", GenomeID => $genome_id);
	if ($pid > 0)
	{
		my $pop_name = $dbh->selectrow_array ("select PopulationName
					from Population where PopulationID=$pid")
				or croak $dbh->errstr;
		
		warn "Unable to delete '$genome_name', it is being referenced by the population '$pop_name'";
		return;
	}
	
	# first retrieve the markers that belong to the genome
	my $sth = $dbh->prepare ("select MarkerID from GenomeMarkers
					where GenomeID=$genome_id")
			or die $dbh->errstr;
	$sth->execute or die $sth->errstr;
	
	# now delete the markers...
	my $array_ref = $sth->fetchall_arrayref();
	foreach my $row (@$array_ref)
	{
		my ($marker_id) = @$row;
		_deleteMarker ($marker_id, $dbh);
	}
	
	# next retrieve the traits that belong to the genome
	$sth = $dbh->prepare ("select TraitID from GenomeTraits
					where GenomeID=$genome_id")
			or die $dbh->errstr;
	$sth->execute or die $sth->errstr;
	
	# delete those too
	$array_ref = $sth->fetchall_arrayref();
	foreach my $row (@$array_ref)
	{
		my ($trait_id) = @$row;
		_deleteTrait ($trait_id, $dbh);
	}
	
	# finally, delete the genome's record from the db
	$dbh->do ("delete from Genome where GenomeID=$genome_id")
		or die $dbh->errstr;
}

sub _deleteMarker ($$)
{
	my $marker_id = shift;
	my $dbh = shift;
	
	$dbh->do ("delete from MarkerAlleleCodes where MarkerID=$marker_id")
		or die $dbh->errstr;
		
	$dbh->do ("delete from Markers where MarkerID=$marker_id")
		or die $dbh->errstr;
	
	$dbh->do ("delete from GenomeMarkers where MarkerID=$marker_id")
		or die $dbh->errstr;
}

sub _deleteTrait ($$)
{
	my $trait_id = shift;
	my $dbh = shift;
	
	$dbh->do ("delete from TraitAlleleEffects where TraitID=$trait_id")
		or die $dbh->errstr;
		
	$dbh->do ("delete from Traits where TraitID=$trait_id")
		or die $dbh->errstr;
		
	$dbh->do ("delete from GenomeTraits where TraitID=$trait_id")
		or die $dbh->errstr;
}

1;
