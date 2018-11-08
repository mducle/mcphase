#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}



unless ($#ARGV >=0) 

{print "\n program formfactor used to calculate the neutron formfactor for an ion \n in dipole approximation - parameters taken from single ion paramater file *.*\n";

 print "\n usage: formfactor  *.*  \n *.* .. filenname\n\n";

   print "Note the formfactor in dipole approximation is different for rare earth and transition metals\n";

   print "FRE(Q)= <j0(Q)>+<j2(Q)> (2/GJ - 1)\n";

   print "FTM(Q)= <j0(Q)>-<j2(Q)> (2/GJ - 1)\n";

 exit 0;}

 

  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

  while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {}

       else{

           if($line=~/^.*IONTYPE\s*=/) {($IONTYPE)=($line=~m|IONTYPE\s*=\s*(.+)|);} 

           if($line=~/^.*GJ\s*=/) {($gJ)=($line=~m|GJ\s*=\s*(.+)|);} 

           if($line=~/^.*FFj0A\s*=/) {($FFj0A)=($line=~m|FFj0A\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj0a\s*=/) {($FFj0a)=($line=~m|FFj0a\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj0B\s*=/) {($FFj0B)=($line=~m|FFj0B\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj0b\s*=/) {($FFj0b)=($line=~m|FFj0b\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj0C\s*=/) {($FFj0C)=($line=~m|FFj0C\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj0c\s*=/) {($FFj0c)=($line=~m|FFj0c\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj0D\s*=/) {($FFj0D)=($line=~m|FFj0D\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj2A\s*=/) {($FFj2A)=($line=~m|FFj2A\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj2a\s*=/) {($FFj2a)=($line=~m|FFj2a\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj2B\s*=/) {($FFj2B)=($line=~m|FFj2B\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj2b\s*=/) {($FFj2b)=($line=~m|FFj2b\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj2C\s*=/) {($FFj2C)=($line=~m|FFj2C\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj2c\s*=/) {($FFj2c)=($line=~m|FFj2c\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           if($line=~/^.*FFj2D\s*=/) {($FFj2D)=($line=~m|FFj2D\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 

           }

      }

      close Fin;

   print ">\n";

   print "IONTYPE=$IONTYPE\n";

   if($gJ==0){ print "Lande factor GJ not found in single ion property file - trying to deduce GJ from IONTYPE ...\n";

              if ($IONTYPE=~/.*Ce3+/){$gJ=6.0/7;}

              if ($IONTYPE=~/.*Pr3+/){$gJ=4.0/5;}

              if ($IONTYPE=~/.*Nd3+/){$gJ=8.0/11;}

              if ($IONTYPE=~/.*Pm3+/){$gJ=3.0/5;}

              if ($IONTYPE=~/.*Sm3+/){$gJ=2.0/7;}

              if ($IONTYPE=~/.*Eu3+/){$gJ=0.0;}

              if ($IONTYPE=~/.*Gd3+/){$gJ=2.0;}

              if ($IONTYPE=~/.*Tb3+/){$gJ=3.0/2;}

              if ($IONTYPE=~/.*Dy3+/){$gJ=4.0/3;}

              if ($IONTYPE=~/.*Ho3+/){$gJ=5.0/4;}

              if ($IONTYPE=~/.*Er3+/){$gJ=6.0/5;}

              if ($IONTYPE=~/.*Tm3+/){$gJ=7.0/6;}

              if ($IONTYPE=~/.*Yb3+/){$gJ=8.0/7;}

              if ($IONTYPE=~/.*U4+/){$gJ=4.0/5;}

              if ($IONTYPE=~/.*U3+/){$gJ=8.0/11;}

              if ($IONTYPE=~/.*Np4+/){$gJ=8.0/11;}

              if ($IONTYPE=~/.*Nd2+/){$gJ=3.0/5;}

              if ($IONTYPE=~/.*Sm2+/){$gJ=0.0;}

              if ($IONTYPE=~/.*Eu2+/){$gJ=2.0;}

              if ($IONTYPE=~/.*Gd2+/){$gJ=3.0/2;}

              if ($IONTYPE=~/.*Tb2+/){$gJ=4.0/3;}

              if ($IONTYPE=~/.*Dy2+/){$gJ=5.0/4;}

              if ($IONTYPE=~/.*Ho2+/){$gJ=6.0/5;}

              if ($IONTYPE=~/.*Er2+/){$gJ=7.0/6;}

              if ($IONTYPE=~/.*Tm2+/){$gJ=8.0/7;}

              

                 if ($IONTYPE=~/.*Co/){$gJ=2;}

                 if ($IONTYPE=~/.*Co+/){$gJ=2;}

                 if ($IONTYPE=~/.*Co2+/){$gJ=2;}

                 if ($IONTYPE=~/.*Co3+/){$gJ=2;}

                 if ($IONTYPE=~/.*Co4+/){$gJ=2;}

                 if ($IONTYPE=~/.*V2+/){$gJ=2;}

                 if ($IONTYPE=~/.*V3+/){$gJ=2;}

                 if ($gJ==0) {print "ERROR: ion not recognized\n exiting ...\n"; exit(1);}

            }

   print "#Note the formfactor in dipole approximation is different for rare earth and transition metals\n";

   print "#FRE(Q)= <j0(Q)>+<j2(Q)> (2/GJ - 1)\n";

   print "#FTM(Q)= <j0(Q)>-<j2(Q)> (2/GJ - 1)\n";

   print "#|Q|(1/A)  FRE(Q) |FRE(Q)|^2 FTM(Q) |FTM(Q)|^2\n";

  for($Q=0;$Q<10;$Q+=0.1)

   {$s=$Q/4/3.1415;

   $j0 = $FFj0A * exp(-$FFj0a * $s * $s) + $FFj0B  * exp(-$FFj0b * $s * $s);

   $j0 = $j0 + $FFj0C * exp(-$FFj0c * $s * $s) + $FFj0D;

     

   $j2 = $FFj2A*$s*$s*exp(-$FFj2a*$s*$s)+$FFj2B*$s*$s*exp(-$FFj2b*$s*$s);

   $j2 = $j2 + $FFj2C * $s * $s * exp(-$FFj2c * $s * $s) + $s * $s * $FFj2D;

   $F=$j0 + $j2 * (2 / $gJ - 1);

   $FF=$F*$F;

   $FTM=$j0 - $j2 * (2 / $gJ - 1);

   $FFTM=$FTM*$FTM;

   print "$Q $F $FF $FTM $FFTM\n";

   }

  }

