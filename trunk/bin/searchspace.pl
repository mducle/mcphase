#!/usr/bin/perl
BEGIN{@ARGV=map{glob($_)}@ARGV}
use File::Copy;

# copy parameter files to specific locations
# remember a few useful perl commands 
# copy("alter_name","neuer_name");
# rename("alter_name","neuer_name");
# chdir('../.');

unless ($#ARGV >0) 
{print " program searchspace used to scan parameter space\n";
 print " usage: searchspace i  * \n";
 print " i ... level of parameter search (i=0 for first scan, for i>0 a\n";
 print " file searchspace.i with parameter sets must exist)\n";
 print " * .. filename(s) of paramter file(s)\n";
 print " \n... there must exist a program calcsta, the output of which\n";
 print " contains 'sta = 249' - this program is called from this output\n";
 print " at every iteration step and  the standard sta deviation is recorded in file\n";
 print " searchspace.i  i=1,2,3, .... \n";
 print " Note: at each level searchspace covers parameterspace with a more dense net of\n";
 print " points. The minimum distance of these points is taken from the parameter-stepwidths\n";
 print " <Press enter to close>";$in=<STDIN>;
 exit 0;}

#in dos version
# mv -> rename
# cp -> copy
# .\calcsta -> calcsta.bat
# rm -> del


#********************************************************************************
 #load from files * parameters, stepwidths, interval and copy original files to *.fit
#format STDOUT_TOP =
#parameter   [ value,     min,      max,     err,       stp     ]   
#.
format STDOUT =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sprintf("%s [%+e,%+e,%+e,%+e,%+e]",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i])
.
format Fout =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sprintf ("%s [%+e,%+e,%+e,%+e,%+e]",$parnam[$ii],$par[$ii],$parmin[$ii],$parmax[$ii],$parerr[$ii],$parstp[$ii])
.

 @parnam=();@par=();@parmin=();@parmax=();@parerr=();@parstp=();@parav=();
  $searchlevel=$ARGV[0]; shift @ARGV;
 while(!open(Fout,">results/searchspace.status")){print "Error opening file results/searchspace.status\n";<STDIN>;}
  print Fout "parameter[value,      min,           max,           (not used)   ,minimum meshwidth]\n";
 foreach (@ARGV)
 {$file=$_; if(mycopy ($file,$file.".bak")){print "\n error copying $file \n";print " <Press enter to close>";$in=<STDIN>;exit 1;} 
   unless (open (Fin, $file.".forfit")){print "\n error:unable to open $file\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}   
   while($line=<Fin>)
{while ($line=~/^(#!|[^#])*?\bpar\w+\s*\Q[\E/) {++$#par;#load another parameter
				 ($parname)=($line=~m/(?:#!|[^#])*?\b(par\w+)\s*\Q[\E/);
                                 foreach(@parnam){if (/$parname/){print "ERROR searchspace: parameter $parname occurs more than one time in input files\n"; print " <Press enter to close>";$in=<STDIN>;exit 1;}}
                                 $parnam[$#par]=$parname;
				 ($par[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*([^,]+)/);
				 ($parmin[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*([^,]+)/);
				 ($parmax[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)/);
				 ($parerr[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^,]+)/);
				 ($parstp[$#par])=($line=~m/(?:#!|[^#])*?\bpar\w+\s*\Q[\E\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*[^,]+\s*,\s*([^\Q]\E]+)/);
                                 $i=$#par;write STDOUT;$ii=$i;write Fout;
				  #check if parmin<=parmax
                          if ($parmin[$#par]>$parmax[$#par]) 
                            {print "ERROR searchspace reading parameterrange: parmin > parmax\n";
                              print " <Press enter to close>";$in=<STDIN>;exit 1;
                             }

                         $line=~s/(?:#!|[^#])*?\bpar\w+\s*\Q[\E//;                     
				 }
     } close Fin;
 }  
    if ($#par<0) {print "Error searchspace: no parameters found in input files @ARGV\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}
    print ($#par+1);print " parameters found - starting search of level $searchlevel ...\n";close Fout;
#*******************************************************************************
$minnumber=1;$pointcounter=0;
$starttime=time;$stamin=1e10;
# open and initialize output file
print "storing points in file searchspace.$searchlevel\n";
open (Foutlevel,">results/searchspace.$searchlevel");
print Foutlevel "#"; 
foreach(@parnam){print Foutlevel $_." ";} 
print Foutlevel "sta\n";

if ($searchlevel==0)
{$pointstocalculate=($#par+1)*2;
            $i=0;foreach(@par){$par[$i]=$parmin[$i]+($parmax[$i]-$parmin[$i])/2;++$i}
            ($sta)=sta(); #sta at input parameter position
       $staorigin=$sta;$i=0;$minimum=1;
       foreach(@par)
        {$dpar=($parmax[$i]-$parmin[$i])/2**($searchlevel+2);
         if (abs($dpar)>abs($parstp[$i]))
         {
          $par[$i]+=$dpar;
          ($sta)=sta();
           read_write_statusfile();
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                        mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
			}

          if ($sta<=$staorigin){$minimum=0;}
          $ii=0;foreach(@par){$dd=sprintf("%e ",$par[$ii]);print Foutlevel $dd;++$ii} print Foutlevel $sta."\n";
	    ++$pointcounter; 
          $par[$i]-=2*$dpar;
          ($sta)=sta();
          read_write_statusfile();
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                        mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
			}
          if ($sta<=$staorigin){$minimum=0;}
          $ii=0;foreach(@par){$dd=sprintf("%e ",$par[$ii]);print Foutlevel $dd;++$ii} print Foutlevel $sta."\n";
          $par[$i]+=$dpar;
  	    ++$pointcounter;
         }
         last if ($sta==0);
        ++$i;
       }
       if ($pointcounter!=0)
       { if ($minimum==1&&$sta!=0)
        {#the input parameters are probably a (local) minimum - so save a file
          ($sta)=sta();
#        last if ($sta==0);
       if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                        mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
			}
          foreach (@ARGV){$file=$_; mycopy ($file.".forfit",$file.".$searchlevel.$minnumber");}
	  ++$minnumber;
       }
        }
}
else
{print "reading points from file searchspace.".($searchlevel-1)."\n";
 $pointstocalculate=0;open (Fin1,"results/searchspace.".($searchlevel-1));
 while($line=<Fin1>){if ($line=~/^\s*#/) {;} else {++$pointstocalculate;}}close Fin1;
 $pointstocalculate*=2*($#par+1);
 open (Fin1,"results/searchspace.".($searchlevel-1));
 while($line=<Fin1>)
  {
   if ($line=~/^\s*#/) {;}
   else{$line=~s/D/E/g;@numbers=split(" ",$line);
            $i=0;foreach(@par){$par[$i]=$numbers[$i];++$i}
            $sta=$numbers[$i]; #sta at input parameter position
       $staorigin=$sta;$i=0;$minimum=1;
       foreach(@par)
        {$dpar=($parmax[$i]-$parmin[$i])/2**($searchlevel+2);
         if (abs($dpar)>abs($parstp[$i]))
         {
          $par[$i]+=$dpar;
          ($sta)=sta();
          open(Fin,"results/searchspace.status");$line=<Fin>;
          if ($line=~/exiting searchspace/){$sta=0;close Fin;}
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                        mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
			}

          if ($sta<=$staorigin){$minimum=0;}
          $ii=0;foreach(@par){$dd=sprintf("%e ",$par[$ii]);print Foutlevel $dd;++$ii} print Foutlevel $sta."\n";
	   ++$pointcounter; 
          $par[$i]-=2*$dpar;
          ($sta)=sta();
          read_write_statusfile();
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                        mycopy ($file.".forfit ",$file.".forfit.min.$searchlevel");}
			}
          if ($sta<=$staorigin){$minimum=0;}
          $ii=0;foreach(@par){$dd=sprintf("%e ",$par[$ii]);print Foutlevel $dd;++$ii} print Foutlevel $sta."\n";
          $par[$i]+=$dpar;
	   ++$pointcounter;
          }
         last if ($sta==0);
        ++$i;
       }
         last if ($sta==0);
         last if ($pointcounter==0);
       if ($minimum==1&&$sta!=0)
        {#the input parameters are probably a (local) minimum - so save a file
          ($sta)=sta();if ($sta>$staorigin*1.00000001){$warning="#WARNING - calculation of sta gave another result at some point than stored in file searchspace.".($searchlevel-1)."\n";}
       if ($sta<$stamin&&$sta!=0){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; mycopy ($file,$file.".min.".$searchlevel);
                                        mycopy ($file.".forfit",$file.".forfit.min.$searchlevel");}
			}
          foreach (@ARGV){$file=$_; mycopy ($file.".forfit",$file.".$searchlevel.$minnumber");}
	  ++$minnumber;
        }
       }
    
  }
  close Fin1;
}
$hours=(time-$starttime)/3600;
$estimate=sprintf("%6.2f",$hours*2*($#par+1));
print "$pointcounter points calculated in  $hours h.\n Time estimate for next level ".($searchlevel+1).": $estimate h\n";
print Foutlevel "#$pointcounter points calculated in  $hours h.\n# Time estimate for next level ".($searchlevel+1).": $estimate h\n";
close Foutlevel;
 foreach (@ARGV) {$file=$_; mycopy($file.".min.$searchlevel",$file);}
     open(Fout,">results/searchspace.status");
print Fout ".......... searchspace stopped...............\n"; 
print Fout ".............................................\n";
print Fout "$pointcounter points calculated in  $hours h.\n";
print Fout ".............................................\n";
print Fout "Mininimal found sta=$stamin\n";
print Fout ".............................................\n";
print Fout "Time estimate for next level ".($searchlevel+1).": $estimate h\n";
print Fout ".............................................\n";
close Fout;
exit 0;

# END OF MAIN PROGRAM
#****************************************************************************** 


sub sta {local $SIG{INT}='IGNORE'; 
 #print "#write modified parameterset to files *\n";
 foreach (@ARGV)
 {$file=$_; open (Fin, $file.".forfit");open (Fout1, ">".$file);open (Fout2,">results/searchspace.par");
   while($line=<Fin>)
     {$modline=$line;
      if ($line=~/^(#!|[^#])*?\bpar/) {#here write modified parameter set to line
                            while ($line=~/^(#!|[^#])*?\bpar\w+\s*\Q[\E/) 
			          {my $i=0;#insert a parameter
				   foreach (@parnam)
				    {$pnam=$_;
				     if ($line=~/^(#!|[^#])*?\b$pnam\s*\Q[\E/)
                                        {$line=~s|$pnam\s*\Q[\E[^\Q]\E]*\Q]\E|$par[$i] |;
                                   $dd=sprintf("%s [%e,%e,%e,%e,%e]",$parnam[$i],$par[$i],$parmin[$i],$parmax[$i],$parerr[$i],$parstp[$i]);
                                   $modline=~s|$pnam\s*\Q[\E[^\Q]\E]+\Q]\E|$dd|;
					}
				     ++$i;
				    }
				  }

                             while ($line=~/^(?:#!|[^#])*?\bfunction\s*\Q[\E(.*?)\Q]\E/) 
			          {$expression=" ".$1." ";
				   #substitute into the expression the paramter values
				   my $i=0;#insert a function parameter
				   foreach (@parnam)
				    {$pnam=$_;
				     while ($expression=~/.*\W$pnam\W/)
                                        {$expression=~s|$pnam|$par[$i]|;
					}
				     ++$i;
				    }
				    # calculate the expression by a little perl program
				    open (Foutcc, ">./results/ccccccc.ccc");
				    printf Foutcc "#!/usr/bin/perl\nprint ".$expression.";\n";
				    close Foutcc;$systemcall="perl ./results/ccccccc.ccc > ./results/cccccc1.ccc";
                                    if ($^O=~/MSWin/){$systemcall=~s|\/|\\|g;}
				    if(system $systemcall){print "error evaluating expression in results/ccccc*";print " <Press enter to close>";$in=<STDIN>;exit 1;}
				    open (Fincc,"results/cccccc1.ccc");
				    $data=<Fincc>; close Fincc;
				    mydel("./results/ccccccc.ccc");mydel("./results/cccccc1.ccc");
				    # $data contains now the result of the mathematical expression
                                    $line=~s|function\s*\Q[\E.*?\Q]\E|$data |;
				   }



                            } print Fout1 $line;print Fout2 $modline;
     } close Fin;close Fout1;close Fout2;
     if (mycopy("results/searchspace.par",$file.".forfit"))
     {print "\n error copying results/searchspace.par to  $file.forfit\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}
 }
 # print "#call routine calcsta to calculate standard deviation\n";
 if ($^O=~/MSWin/){
                   if(system ("calcsta.bat > results\\searchspace.sta")){print "\n error executing calcsta.bat\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}
                  }
 else
                  {
                   if(system ("./calcsta > results/searchspace.sta")){print "\n error executing calcsta.bat\n";print " <Press enter to close>";$in=<STDIN>;exit 1;}
                  }

 open (Fin,"./results/searchspace.sta"); 
 while($line=<Fin>){
           if($line=~/^(#!|[^#])*?\bsta\s*=/) {($sta)=($line=~m/(?:#!|[^#])*?\bsta\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)/);} 
                   }
 close Fin;
 mydel("./results/searchspace.sta");
 return $sta;
}  

sub mycopy { my ($file1,$file2)=@_;

             if ($^O=~/MSWin/){$file1=~s|\/|\\|g;$file2=~s|\/|\\|g;
                               return system("copy ".$file1." ".$file2);
                              }
                 else
                              {return system("cp -f ".$file1." ".$file2);
                              }

           }

sub mydel  { my ($file1)=@_;

             if ($^O=~/MSWin/){$file1=~s|\/|\\|g;
                               return system("del ".$file1);
                              }
                 else
                              {return system("rm ".$file1);
                              }

           }
    
sub read_write_statusfile {
          open(Fin,"./results/searchspace.status");$line=<Fin>;close Fin;
          if ($line=~/exiting searchspace/){$sta=0;}
          else
     {open(Fout,">./results/searchspace.status");$ii=0;
     $est=sprintf("%6.2f",($pointstocalculate-$pointcounter)*(time-$starttime)/3600/($pointcounter+1));
     $passed=sprintf("%6.2f",(time-$starttime)/3600);
     print Fout "Searchlevel $searchlevel   Current sta=$sta Minimum sta=$stamin\n";
     print Fout "--------------------------------------------------------------------\n";
     print Fout "Time since Start: $passed h Estimated time to complete: $est h\n";
     print Fout "--------------------------------------------------------------------\n";
     print Fout "parameter[value,      min,           max,           (not used)   ,minimum meshwidth]\n";
     foreach(@par){write Fout;++$ii;}
     print Fourt $warning;
     close Fout;}

                          }