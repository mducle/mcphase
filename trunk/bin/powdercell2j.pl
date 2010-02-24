#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}

#function make_input_mcphas
#% Lucian G Pascut
#% pascutlucian@yahoo.com

unless ($#ARGV >=0)
{print STDERR << "EOF";

    $0:

    program $0  used to create mcphas.j type file from output of powdercell

    usage: $0 *.*

    *.* .. filenname of powdercell file
    output is written to mcphas.j

 In case of any error please check the format of the input file *.*

 Example:
  No   name       crystal coordinates          cartesian coordinates
                 x        y        z           x        y        z
  ------------------------------------------------------------------
  1     Sr1    0.3644   0.0000   0.2500     1.0962  -4.1497  -2.7991
  ...

EOF
 exit 0;
}
    $p=0;
  foreach (@ARGV)
  {   $file=$_;
 #   %
 #   %begin:"find how many different types of atoms do we have in the structure"
 #   %
 #   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    unless (open (Fin, $file)){die "\n error:unable to open $file\n";}
    while($line=<Fin>)
     {
      if ($line=~/^\s*\d/) {@c=split(" ",$line);
                            $atomtype[$p]=$c[1];
                            $x[$p]=$c[2];$y[$p]=$c[3];$z[$p]=$c[4];
                            $p++;
                           }
     }
    $pp=0;
    for($i=0;$i<$p;++$i)
    { $ppp=1;
      foreach(@atomtypetempp)
       {if($_=~$atomtype[$i]){$ppp=0;}
       }
      if($ppp==1){$atomtypetempp[$pp]=$atomtype[$i];$pp+=$ppp;}
    }
    print "\n   There are ".($pp)." different atom types:\n";
#    %find the number of atoms in the same group
    $noofatom=0;$i=0;
      foreach(@atomtypetempp)
       {foreach(@atomtype){
        if($_=~$atomtypetempp[$i])        {++$nofatoms[$i];}
        }
        ++$i;
       }
    print "@atomtypetempp\n";
    print "    Number of different atoms for each type:\n";
    print "@nofatoms\n";
  }
  print "Give the lattice constants a(A) and then press enter\n";
  $alattice= <STDIN>;$alattice=~s/\n//;
  print "Give the lattice constants b(A) and then press enter\n";
  $blattice=<STDIN>;$blattice=~s/\n//;
  print "Give the lattice constants c(A) and then press enter\n";
  $clattice=<STDIN>;$clattice=~s/\n//;
  print "Give the angle alpha(deg) and then press enter\n";
  $alpha=<STDIN>;$alpha=~s/\n//;
  print "Give the angle beta(deg) and then press enter\n";
  $beta=<STDIN>;$beta=~s/\n//;
  print "Give the angle gamma(deg)  and then press enter\n";
  $gamma=<STDIN>;$gamma=~s/\n//;

 
  foreach(@atomtypetempp)
  {  open (Fout1, ">$_.sipf");
     print Fout1 "#IONTYPE=$_\n";
     print "For $_ atom the oxidation state is? \nPlease give the charge in units of elementary charge |e|; for ex. -2 or +1 or +3:\n";
     $charge=<STDIN>;$charge=~s/\n//;
     print Fout1 "CHARGE=$charge\n";
     close Fout1;
  }
    open (Fout, ">mcphas.j");
  print "\n\n\n\n\n\n";
  print Fout "#<!--mcphase.mcphas.j-->\n";
  print Fout "#***************************************************************\n";
  print Fout "# Lattice and Exchange Parameter file for\n";
  print Fout "# mcphas version 3.0\n";
  print Fout "# - program to calculate static magnetic properties\n";
  print Fout "# reference: M. Rotter JMMM 272-276 (2004) 481\n";
  print Fout "# mcdisp version 3.0\n";
  print Fout "# - program to calculate the dispersion of magnetic excitations\n";
  print Fout "# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n";
  print Fout "#***************************************************************\n";
  print Fout "#\n";
  print Fout "# ";
  $i=0;foreach(@atomtypetempp)
  {print Fout " $_ (".($nofatoms[$i]/min(@nofatoms)).")";
   print " $_ (".($nofatoms[$i]/min(@nofatoms)).")";++$i;
  }

print STDOUT << "EOF";


 If this is not the formula of your system
 please check the format of the input file

Example:
 No   name       crystal coordinates          cartesian coordinates
                x        y        z           x        y        z
 ------------------------------------------------------------------
 1     Sr1    0.3644   0.0000   0.2500     1.0962  -4.1497  -2.7991
 ...
EOF
 my $total = 0;
($total+=$_) for @nofatoms;
  print Fout "\n";
  print Fout "#\n";
  print Fout "# Lattice Constants (A)\n";
  print Fout "#\n";
  print Fout "#! a= $alattice b= $blattice c= $clattice  alpha= $alpha beta= $beta gamma= $gamma\n";
  print Fout "#\n";
  print Fout "#! r1a=   1 r2a= 0 r3a=  0\n";
  print Fout "#! r1b=   0 r2b= 1 r3b=  0   primitive lattice vectors [a][b][c]\n";
  print Fout "#! r1c=   0 r2c= 0 r3c=  1\n";
  print Fout "#\n";
  print Fout "#! nofatoms= $total  nofcomponents=6  number of atoms in primitive unit cell/number of components of each spin\n";

$i=0;$nr=1;
  foreach(@atomtypetempp)
   {$ii=0;
    foreach(@atomtype)
      {  if($_=~$atomtypetempp[$i])
           {
           print Fout "#********************************************************************* \n";
           print Fout "#ATOM TYPE $_ ; number of the atom in the UNIT CELL = $nr ; number of the atom within this type = $nofatoms[$i]\n";
           print Fout "#! da= $x[$ii] [a] db= $y[$ii] [b] dc= $z[$ii] [c] nofneighbours=0 diagonalexchange=1 gJ= 0 cffilename= $_.sipf\n";
           ++$nr; }
       ++$ii;
      }
    ++$i;
   }
close Fout;

 
print "\n\n";
print "The folowing file have been created by this program:\n\n";
print "mcphas.j - the input file for McPhase program\n";
print "\n";
foreach(@atomtypetempp) {print "$_.sipf - contain the oxidation state of $_\n";}

print STDOUT << "EOF";

 running \"makenn R\" command in McPhase will create a number of files 
 \"makenn.aN.pc\" equal with the number of the atoms in the unit cell.
 in this file you find the neighbours of the N atom in the unit cell.
 It is important to know what number corresponds to each atom type!

EOF


sub min
{
    my (@numbers);

    @numbers = @_;

    my ($min);

    $min = $numbers[0];

    foreach my $i (@numbers)
    {
        if ($i < $min)
        {
            $min = $i;
        }
    }

    return ($min);
}

sub sum
{
    my (@numbers);

    @numbers = @_;

    my ($sum);



    foreach my $i (@numbers)
    { $sum+=$i;
    }

    return ($min);
}