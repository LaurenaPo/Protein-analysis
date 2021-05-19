#!/usr/bin/perl

#REFERENCE
#   Please cite:
#   Mansiaux Y, Joseph AP, Gelly J-C, de Brevern AG (2011) Assignment of 
#   PolyProline II conformation and analysis od sequence-structure 
#   relationship. PLoS One 6(3): e18401. doi:10.1371/journal.pone.0018401 
#
#   Chebrek R, Leonard S, de Brevern AG, Gelly J-C (2014) 
#   PolyprOnline: polyproline helix II and secondary structure assignment database.
#   Database ; Nov 7;2014 [pmid:25380779]


# Copyright Jean-Christophe Gelly (Jan 20 2012)
#
# jean-christophe.gelly@univ-paris-diderot.fr
#
# This software is a computer program whose purpose is to 
# assign polyproline type II helix from dsspcmbi program output.
#
# This software is governed by the CeCILL-B license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL-B
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL-B license and that you accept its terms.

use strict;

#my $path_dssp="/home/iena/gelly/SS_ASSIGN/dssp/dsspcmbi"; #Add the absolute path to the dssp program (eg: my $path_dssp="/home/jdoe/bin/dsspcmbi";)
my $path_dssp="/usr/bin/dssp"; # !!!!!! ADD PATH TO DSSP INSIDE THE "" !!!!!!!!!

my $ref_tab_out_dssp;

my $pathfile=shift;

my $flag_horiz=0;
if($#ARGV == 0)
{
	if($ARGV[0] eq '-horiz')
	{
		$flag_horiz=1;
	}
}

if ($path_dssp eq "")
{
    print STDERR "Error: No path to dssp program !\n";
    print STDERR "Please open this program ($0) and modify the line 50 and give the complete path to the program dssp for \$path_dssp variable\n";
    exit 1;
}

unless(-e "$path_dssp") 
{
    print STDERR "Error: Path to dssp program ($path_dssp) is not correct!\nPlease check it.\n";
    exit 1;
}

if ($pathfile =~ /^-h$/ and $pathfile =~/^--h$/ and $pathfile =~/^--help$/ and $pathfile =~/^-help$/)
{
    help();
}



unless(-e "$pathfile") 
{
    print STDERR "Error: Not found input PDB file $pathfile!\n\n";
    help();
    
}

$ref_tab_out_dssp = &dssp_ppii_assignment("$path_dssp",$pathfile);


if($flag_horiz == 0)
{
	foreach my $line (@{$ref_tab_out_dssp})
	{
    		print "$line";
	}
}
elsif ($flag_horiz == 1)
{
	my $ref_tab_out_dssp_horiz=Parsing_Dssp($ref_tab_out_dssp);
	foreach my $line (@{$ref_tab_out_dssp_horiz})
	{
    		print "$line";
	}
}

exit 1;


############
# SUBs
#############

sub dssp_ppii_assignment{

    my $bin_dssp = shift;
    my $PDB      = shift;

    my $aa     = "";
    my $ss     = "";
    my $assign = 0;
    my $i      = 0;
    my $idx    =0;

    my %hash_chain;

    my $chain="";

    my $epsilon        = 29;
    my $cannonical_phi =-75;
    my $cannonical_psi = 145;

    my $sup_phi=$cannonical_phi+$epsilon;
    my $inf_phi=$cannonical_phi-$epsilon;

    my $sup_psi=$cannonical_psi+$epsilon;
    my $inf_psi=$cannonical_psi-$epsilon;

    my @tab_new_dssp;

    #Launch DSSP
    my $run_dssp="$bin_dssp $PDB";
    my @tab_output=`$run_dssp`;


    #Parsing DSSP to detect Polyproline II Helix
    WHILE1:while (my $line=$tab_output[$i])
    {
        if($line =~ "#  RESIDUE AA STRUCTURE")
        {
            $assign = 1;

        }
        elsif($assign == 1)
        {
            $ss          = substr($line,16 ,1);
            $aa          = substr($line,13 ,1);

            if($aa =~ /\!/)
            {
                ${ $hash_chain{$chain} }[$idx] = 0;
                $idx++;
                $i++;
                next WHILE1;
            }

            my $index    = substr($line,0,5);
            my $position = substr($line,5,5);

            $chain    = substr($line,11,1);

            my $phi      = substr($line,103,6);
            my $psi      = substr($line,109,6);

            $index=~s/\s//g;
            $phi  =~s/\s//g;
            $psi  =~s/\s//g;

            if ($chain eq ' ') {$chain="_";}

            if (not exists $hash_chain{$chain})
            {
                $idx=0;
                $hash_chain{$chain}=();
            }

            if ($phi <= $sup_phi and $phi >= $inf_phi and $psi <= $sup_psi and $psi >= $inf_psi and $ss eq ' ')
            {
                ${ $hash_chain{$chain} }[$idx]=1;
                if ($idx != 0) 
                { 
                    if ( ${ $hash_chain{$chain} }[$idx - 1] >= 1 )
                    {
                        ${ $hash_chain{$chain} }[$idx]       = 2;
                        ${ $hash_chain{$chain} }[$idx - 1]   = 2;
                    }
                }
                else
                {
                    ${ $hash_chain{$chain} }[$idx]=0;
                }
            }
            else
            {
                ${ $hash_chain{$chain} }[$idx]=0;
            }

            $idx++;
        }
        $i++;

    }

    #PRINT MODIFED OUTPUT

    $i=0;
    $idx=0;
    $assign=0;
    WHILE2:while (my $line=$tab_output[$i])
	{
		if($line =~ "#  RESIDUE AA STRUCTURE")
		{
			$assign = 1;
            push @tab_new_dssp,$line;
            $idx=0;
		}
        elsif($assign == 1)
		{
            $aa          = substr($line,13 ,1);
            my $chain    = substr($line,11,1);

            if($aa =~ /\!/)
            {
                $idx++;
                push @tab_new_dssp,$line;
                $i++;
                next WHILE2;
            }

            if ( ${ $hash_chain{$chain} }[$idx] >= 2 )
            {
                substr($line,16 ,1 ,"P");
            }
            push @tab_new_dssp,$line;
            $idx++;
        }
        else
        {
            push @tab_new_dssp,$line;
        }
        $i++;
    }

    # return ref of a tab
    return(\@tab_new_dssp);
}
sub Parsing_Dssp{

    my $ref_tab_out_dssp=shift;
    my $aa = "";
    my @AA = ();

	my @tab_out;
	my $ss = "";
	my @SS = ();
	my $buff = 1;
	my $assign = 0;
	my $i = 0;
	my $j = 0;
my %dsspto3=
(
    "H" => "H",
    "G" => "G",
    "I" => "I",
    "E" => "E",
    "B" => "B",
    "C" => "C",
    "S" => "S",
    "T" => "T",
    " " => "-",
    "P" => "P"
);


	push(@tab_out,sprintf(">SEQ\n"));

	foreach my $line (@{$ref_tab_out_dssp})
	{
		if($line =~ "#  RESIDUE AA STRUCTURE")
		{
			$assign = 1;

		}
        elsif($assign == 1)
		{
			$ss = substr($line,16,1);
            		$aa = substr($line,13,1);
			my $chain=substr($line,11,1);
			#print "AA:$aa SS:$ss \n$line\n";
			if($aa !~ /\!/)
			{

				#if( ($buff+1) < substr($line,7,))
				#{			
				#	push(@SS,"*");
				#}
				push(@SS,$dsspto3{$ss});
                		push(@AA,$aa);
				$buff = substr($line,7,3);
			}

		}
	}

    push(@tab_out,sprintf(join('',@AA)."\n"));
    push(@tab_out,sprintf(">DSSPPII\n"));
    push(@tab_out,sprintf(join('',@SS)."\n"));

	return(\@tab_out);
}



##############
# Sub help
##############
sub help 
{
    print STDERR "
USAGE
    $0 <PDB_File> -> Read PDB_File and print DSSP output (with PPII helix) on stdout

OPTIONS
    -horiz : Give output DSSP in 1D sequence fashion 

ADDITIONAL INFOS
    To use this program you must open it to modify the line 40
    and give the complete path to the program dssp (\$path_dssp variable).

LICENCE
    Copyright Jean-Christophe Gelly (Jan 27 2012)

    jean-christophe.gelly\@univ-paris-diderot.fr

    This software is a computer program whose purpose is to 
    assign polyproline type II helix from dsspcmbi program output.

    This software is governed by the CeCILL-B license under French law and
    abiding by the rules of distribution of free software.  You can  use, 
    modify and/ or redistribute the software under the terms of the CeCILL-B
    license as circulated by CEA, CNRS and INRIA at the following URL
    \"http://www.cecill.info\". 

    As a counterpart to the access to the source code and  rights to copy,
    modify and redistribute granted by the license, users are provided only
    with a limited warranty  and the software's author,  the holder of the
    economic rights,  and the successive licensors  have only  limited
    liability. 

    In this respect, the user's attention is drawn to the risks associated
    with loading,  using,  modifying and/or developing or reproducing the
    software by the user in light of its specific status of free software,
    that may mean  that it is complicated to manipulate,  and  that  also
    therefore means  that it is reserved for developers  and  experienced
    professionals having in-depth computer knowledge. Users are therefore
    encouraged to load and test the software's suitability as regards their
    requirements in conditions enabling the security of their systems and/or 
    data to be ensured and,  more generally, to use and operate it in the 
    same conditions as regards security. 

    The fact that you are presently reading this means that you have had
    knowledge of the CeCILL-B license and that you accept its terms.

REFERENCES
   If you find this programm useful please cite:

   Mansiaux Y, Joseph AP, Gelly J-C, de Brevern AG (2011) Assignment of 
   PolyProline II conformation and analysis od sequence-structure 
   relationship. PLoS One 6(3): e18401. doi:10.1371/journal.pone.0018401 

   Chebrek R, Leonard S, de Brevern AG, Gelly J-C (2014) 
   PolyprOnline: polyproline helix II and secondary structure assignment database.
   Database ; Nov 7;2014 [pmid:25380779]


";
exit 1;
}


