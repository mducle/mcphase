#!/usr/bin/perl

unless ($#ARGV >0) 
{print " program simannfit used to perform simulating annealing\n";
 print " usage: simannfit  10 [-t 100][-s 1000] * \n  10 .. initial temperature \n * .. filename(s) of paramter file(s)\n";
 print " ... there must exist a program calcsta, the output of which\n";
 print " contains 'sta = 249' - this program is called from this output\n";
 print " at every iteration step and  the standard sta deviation is minimzed\n";
 print " simannfit gives as a parameter a maximum number for sta - if during the\n";
 print " calculation of sta in calcsta this number is exceeded calcsta can exit (saves time)\n";
 print " at every iteration step and  the standard sta deviation is minimzed\n";
 print " option -t sets time limit until program end (in seconds)\n";
 print " option -s gives maximal number of interation steps to be done\n";
 exit 0;}


#********************************************************************************
 #load from files * parameters, stepwidths, interval and copy original files to *.fit
#format STDOUT_TOP =
#parameter   [ value,     min,      max,     err,       stp     ]   
#.
format STDOUT =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sprintf("%s [%e,%e,%e,%e,%e]",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i])
.

 @parnam=();@par=();@parmin=();@parmax=();@parerr=();@parstp=();@parav=();
				 @parhisto=();@parhistostart=();$hh=0;
  $stattemp=$ARGV[0]; shift @ARGV;
  $starttime=time;$maxtim=1e10;$maxstep=1e24;
  if ($ARGV[0]=~"-t") {shift @ARGV; $maxtim=$ARGV[0]; shift @ARGV;}
  if ($ARGV[0]=~"-s") {shift @ARGV; $maxstep=$ARGV[0]; shift @ARGV;}
  
 foreach (@ARGV)
 {$file=$_; system ("cp -f ".$file." ".$file.".par"); open (Fin, $file);
   while($line=<Fin>)
 {while ($line=~/^.*par\w+\s*\Q[\E/) {++$#par;#load another parameter
				 ($parnam[$#par])=($line=~m|(par\w+)\s*\Q[\E|);
				 ($par[$#par])=($line=~m|par\w+\s*\Q[\E\s*([^,]+)|);
				 ($parmin[$#par])=($line=~m|par\w+\s*\Q[\E\s*[^,]+\s*,\s*([^,]+)|);
				 ($parmax[$#par])=($line=~m|par\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)|);
				 ($parerr[$#par])=($line=~m|par\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)|);
				 ($parstp[$#par])=($line=~m|par\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^\Q]\E]+)|);
                                 $i=$#par;write STDOUT;
                         #check if parmin<=parmax
                          if ($parmin[$#par]>$parmax[$#par]) 
                            {print "ERROR simannfit reading parameterrange: parmin > parmax\n";
                             exit 1;
                             }
				 $line=~s|par\w+\s*\Q[\E||;                     
				 $parhistostart[$i]=$hh;
				 $hh+=int(($parmax[$i]-$parmin[$i])/$parstp[$i])+2;
				 $#parhisto=$hh;
				  for($hx=0;$hx<=int(($parmax[$i]-$parmin[$i])/$parstp[$i])+1;++$hx)
                                   {$parhisto[($hx+$parhistostart[$i])]=0;}
                                  
				 }
     } close Fin;
 } print "starting fit ...\n";
#*******************************************************************************
# fitting loop
$rnd=1;$stasave=1e20;
 ($sta)=sta();$stps=1;$noofupdates=0;$stepnumber=0;
if($sta>0)
{
$SIG{INT} = \&catch_zap;
while($sta>0)
{  $stasave=$sta;
 # modify parameters
 print "\n ...next fitting loop ...\n";
 @parsav=@par;$i=0;
 foreach(@par){$rnd=rand;$par[$i]+=($rnd-0.5)*$parstp[$i]*$stps;
               if ($par[$i]<$parmin[$i]) {$par[$i]=$parmin[$i];}
	       if ($par[$i]>$parmax[$i]) {$par[$i]=$parmax[$i];}
               write STDOUT;
	       ++$i;}
 $rnd=rand;
   ($sta)=sta(); # CALCULATE sta !!!!
   ++$stepnumber;
   print " ...  current sta=$sta, statistical T=$stattemp, step ratio=$stps\nsta of stored parameters=$stasave\n";

 if (time-$starttime>$maxtim){$sta=0;}
 if ($stepnumber>$maxstep){$sta=0;}
 if ($sta==0) {#recover ol pars
      @par=@parsav; 
      ($sta)=sta(); # CALCULATE sta !!!!
      }
 last if ($sta==0);
 if ($sta>$stasave)
  {if($rnd>exp(-($sta-$stasave)/$stattemp))
     {#recover ol pars
      @par=@parsav;#foreach(@parstp){$_*=0.9995;}
      $stps*=0.999;$sta=$stasave; 
      if ($stps<0.01){$stps=10;}# if stepwidth decreased too much make large steps again to get out of side minimum !!!
     }
   else
     {$stattemp=$stattemp*0.999;
     }     
  }  
 else
  {$stattemp=$stattemp*0.995;
   #update errors
   $i=0;
   print "STA DECREASED in the last LOOP:\n";
   foreach(@par){$p=$_;
   $parav[$i]=($parav[$i]*$noofupdates + $p)/($noofupdates+1);
   $parerr[$i]=sqrt($parerr[$i]*$parerr[$i]*$noofupdates+
                   ($p-$parav[$i])*($p-$parav[$i]))/($noofupdates+1);     
   $hx=int(($p-$parmin[$i])/$parstp[$i]);
   ++$parhisto[($hx+$parhistostart[$i])];
    open(Fout,">".$parnam[$i].".hst");
    print Fout "#{Histogram of parameter ".$parnam[$i]."\n# value vs. number of  occurrences in good solutions (sta decreased)}\n";
    for($hx=0;$hx<=int(($parmax[$i]-$parmin[$i])/$parstp[$i])+1;++$hx)
     {print Fout (($hx+0.5)*$parstp[$i]+$parmin[$i])."   ".($parhisto[($hx+$parhistostart[$i])])."\n"; 
     } close Fout;

   ++$i;} ++$noofupdates;
   #printout current parameters
   $i=0;foreach(@par){write STDOUT;++$i;}
   print "   sta of stored parameters=$sta\n";
  }	       
 print "*";
 }	       

   print "best fit:\n";
   $i=0;foreach(@par){write STDOUT;++$i;}
   ($sta)=sta(); # CALCULATE sta !!!!
}
   print "      sta=$sta\n";
# move files
 foreach (@ARGV)
 {$file=$_; system ("mv ".$file." ".$file.".fit");
            system ("mv ".$file.".par ".$file);}
exit 0;
# END OF MAIN PROGRAM
#****************************************************************************** 

sub catch_zap {
 my $signame = shift;
 foreach (@ARGV)
 {$file=$_; system ("mv ".$file." ".$file.".fit");
            system ("mv ".$file.".par ".$file);}
 die "SIG$signame stopping fit\n";
}
#****************************************************************************** 

#use Carp;
# sub swrite {
#             croak "usage: swrite PICTURE ARGS" unless @_;
#                my $format = shift;
#               $^A = "";
#                  formline($format,@_);
#             return $^A;
#            }          

sub sta {local $SIG{INT}='IGNORE'; 
 #print "#write modified parameterset to files *\n";
 foreach (@ARGV)
 {$file=$_; open (Fin, $file.".par");open (Fout1, ">".$file);open (Fout2,">simannfit.par");
   while($line=<Fin>)
     {$modline=$line;
      if ($line=~/^.*par/) {#here write modified parameter set to line
                            while ($line=~/^.*par\w+\s*\Q[\E/) 
			          {$i=0;#insert a parameter
				   foreach (@parnam)
				    {$pnam=$_;
				     if ($line=~/^.*$pnam\s*\Q[\E/)
                                        {$line=~s|$pnam\s*\Q[\E[^\Q]\E]*\Q]\E|$par[$i] |;
#                                   $dd=swrite("[@##.######,@##.####,@##.####,@##.######,@##.######]",
#                                   ;
                                   $dd=sprintf("%s [%e,%e,%e,%e,%e]",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i]);
                                   $modline=~s|$pnam\s*\Q[\E[^\Q]\E]+\Q]\E|$dd|;
					}
				     ++$i;
				    }
				  }
#                             while ($line=~/^.*function\s*\Q[\E/) 
#			          {$i=0;#insert a function parameter
#				   foreach (@parnam)
#				    {$pnam=$_;
#				     if ($line=~/^.*function\s*\Q[\E\s*$pnam/)
#                                        {$line=~s|function\s*\Q[\E\s*$pnam\s*\Q]\E|$par[$i] |;
#					}
#				     ++$i;
#				    }
#				  }

#new method of treating functions
                             while ($line=~/^.*?function\s*\Q[\E(.*?)\Q]\E/) 
			          {$expression=$1;
				   #substitute into the expression the paramter values
				   $i=0;#insert a function parameter
				   foreach (@parnam)
				    {$pnam=$_;
				     while ($expression=~/.*$pnam/)
                                        {$expression=~s|$pnam|$par[$i]|;
					}
				     ++$i;
				    }
				    # calculate the expression by a little perl program
				    open (Foutcc, ">./ccccccc.ccc");
				    printf Foutcc "#!/usr/bin/perl\nprint ".$expression.";\n";
				    close Foutcc;
                            system "chmod 755 ./ccccccc.ccc"; 
				    system "./ccccccc.ccc > ./cccccc1.ccc";
				    open (Fincc,"./cccccc1.ccc");
				    $data=<Fincc>; close Fincc;
				    system "rm ./ccccccc.ccc ./cccccc1.ccc";
				    # $data contains now the result of the mathematical expression
                                    $line=~s|function\s*\Q[\E.*?\Q]\E|$data |;
				   }



                            } print Fout1 $line;print Fout2 $modline;
     } close Fin;close Fout1;close Fout2;
     unless (rename "simannfit.par",$file.".par")
     {die "\n error:perhaps (perl rename cannot cross filesystems) \n";}
 }
# print "#call routine calcsta to calculate standard deviation\n";
$staboundary=$stasave-log($rnd)*$stattemp;
 system ("./calcsta $staboundary > simannfit.sta");
 open (Fin,"simannfit.sta"); 
 while($line=<Fin>){
           if($line=~/^.*sta\s*=/) {($sta)=($line=~m|sta\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 
                   }
 close Fin;
 system ("rm simannfit.sta");
 return $sta;
}  

    
