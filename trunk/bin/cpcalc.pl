#!/usr/bin/perl

$Tmin=$ARGV[0];
$Tmax=$ARGV[1];
$deltaT=$ARGV[2];
unless($#ARGV>=2)
  {print "Program cpcalc - calculates specific heat from output file results/cfield.out\n";
   print "             use as :  cpcalc tmin tmax deltat [-option]\n"; 
   print "             alternatively \n";
   print "             use as: cpcalc col1 col2 datafile [-option]\n";
   print "                     (take cp-data from datafile and calculate\n";
   print "                      standard deviation)\n";
   print "     output is written to stdout, energies have to be given in meV,\n";
   print "     temperatures tmin tmax deltat in Kelvin\n";
   print "     Options: -s   .... calculate entropy  (J/molK) instead of cp\n";
   print "              -f   .... calculate free energy (J/mol) instead of cp\n";
   print "              -u   .... calculate magnetic energy (J/mol) instead of cp\n";
   print "              -z   .... calculate partition sum instead of cp\n";
   exit(1);
  }
      $cptext="cp(J/molK)";
      if($ARGV[3]=~/-f/){$cptext="f(J/mol)";}
       if($ARGV[3]=~/-z/){$cptext="z";}
       if($ARGV[3]=~/-u/){$cptext="u(J/mol)";}
       if($ARGV[3]=~/-s/){$cptext="s(J/molK)";}

unless (open (Fin,"results/cfield.out")) {unless (open (Fin,"cfield.out")){print "ERROR cpcalc: file results/cfield.out not found\n";exit(1);}else{print "#reading cfield.out\n"; }}
else {print "#reading results/cfield.out\n";}

$noflevels=1; # initialize noflevels
# read energies from cfield.out
  while($line=<Fin>)
  {if($line=~/^.*\QEnergy Eigenvalues are in\E/){unless ($line=~/^.*\QEnergy Eigenvalues are in  meV\E/) 
                                                        {print "ERROR cpcalc: energies in cfield.out must be in meV ! - ".$line."\n";exit(1);}
                                                }
   if($line=~/^.*\QNumber of different energy levels\E\s*:/) # read noflevels
      {($noflevels)=($line=~m|\QNumber of different energy levels\E\s*:\s*([\d.eEdD\Q-\E\Q+\E]+)|);
       print "#Number of different energy levels: $noflevels\n";
      } 
   if($line=~/^.*\QEnergy shift\E\s*\Q(Eshift)\E\s*:/) # read energyshift
      {($energyshift)=($line=~m|\QEnergy shift\E\s*\Q(Eshift)\E\s*:\s*([\d.eEdD\Q-\E\Q+\E]+)|);
       print "#Energy shift: $energyshift meV\n";
      } 
   for($i=1;$i<=$noflevels;++$i)
     {if($line=~/^.*\QE( \E$i\Q)\E\s*=/) 
      {($E[$i])=($line=~m|\QE( \E$i\Q)\E\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);
       ($deg[$i])=($line=~m|.*\QDegeneracy\E\s*:\s*([\d.eEdD]+)|);
      print "#E($i)=".$E[$i]."meV degeneracy ".$deg[$i]."-fold\n";
      } 
     }
  }     
  close Fin;

$dT=0.1; # this is fixed value for calculation of derivative of energy

if (open (Fin,$deltaT))
{
$colT=$Tmin;
$colcp=$Tmax;
$ii=0;$sta=0;
# read temperatures and calculate cp + comp to experiment
  print  "# T(K)  $cptext  exp-$cpctext\n";
  while($line=<Fin>)
  {
       if ($line=~/^\s*#/) {}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);
           	$T0=$numbers[$colT-1];
		$cpexp=$numbers[$colcp-1];
            ($cpclc)=cp(); # CALCULATE cp !!!!
	 ++$ii;
	 $sta+=($cpclc-$cpexp)*($cpclc-$cpexp);
         print   "$T $cpclc $cpexp\n";
        }
  }
  $sta/=$ii;
  print  "#sta=$sta\n";       
  close Fin;
}
else
{print  "# T(K)  $cptext\n";
 for($T0=$Tmin;$T0<$Tmax;$T0+=$deltaT)
 {($cpclc)=cp();
  print "$T0  $cpclc\n";
 }
}
exit(0);

sub cp {
$T=$T0+$dT/2;
            $Zp=0;$Up=0;
            for ($i=1;$i<=$noflevels;++$i)
                {$x=$E[$i]/$T/0.0862;
                 $Zp+=$deg[$i]*exp(-$x);
                 $Up+=$deg[$i]*$E[$i]*exp(-$x);
                 }
            $Up/=$Zp*0.0862; # magnetic energy per ion in K
            $Up+=$energyshift/0.0862; # shift energy 
            $Up*=1.38066*6.023; # magnetic energy in J per mol
	    
	    $T=$T0-$dT/2;
            $Zm=0;$Um=0;
            for ($i=1;$i<=$noflevels;++$i)
                {$x=$E[$i]/$T/0.0862;
                 $Zm+=$deg[$i]*exp(-$x);
                 $Um+=$deg[$i]*$E[$i]*exp(-$x);
                 }
            $Um/=$Zm*0.0862; # magnetic energy per ion in K
            $Um+=$energyshift/0.0862; # shift energy 
            $Um*=1.38066*6.023; # 0.0862meV/K*1.602e-22J/meV*6.023e23ion/mol 
                                # magnetic energy in J per mol
         
	 $cpc=($Up-$Um)/$dT; # specific heat 	    
       if($ARGV[3]=~/-f/){$cpc=($energyshift/0.0862-$T0*log(0.5*($Zm+$Zp)))*1.38066*6.023;} # helmholtz function
       if($ARGV[3]=~/-z/){$cpc=0.5*($Zm+$Zp)*exp(-$energyshift/$T0/0.0862);} # partition sum
       if($ARGV[3]=~/-u/){$cpc=0.5*($Um+$Up);} # magnetic energy
       if($ARGV[3]=~/-s/){$cpc=0.5*($Um+$Up)/$T0-($energyshift/$T0/0.0862-log(0.5*($Zm+$Zp)))*1.38066*6.023;} # entropy

return $cpc;
}