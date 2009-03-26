#!/usr/bin/perl

unless ($#ARGV >0) 
{print " program searchspace used to scan parameter space\n";
 print " usage: searchspace i  * \n i ... level of parameter search (i=0 for first scan, for i>0 a file searchspace.i with parameter sets must exist)\n * .. filename(s) of paramter file(s)\n";
 print " ... there must exist a program calcsta, the output of which\n";
 print " contains 'sta = 249' - this program is called from this output\n";
 print " at every iteration step and  the standard sta deviation is recorded in file\n";
 print " searchspace.i  i=1,2,3, .... \n";
 print " Note: at each level searchspace covers parameterspace with a more dense net of\n";
 print " points. The minimum distance of these points is taken from the parameter-stepwidths\n";
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
format Fout =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sprintf ("%s [%e,%e,%e,%e,%e]",$parnam[$ii],$par[$ii],$parmin[$ii],$parmax[$ii],$parerr[$ii],$parstp[$ii])
.

 @parnam=();@par=();@parmin=();@parmax=();@parerr=();@parstp=();@parav=();
  $searchlevel=$ARGV[0]; shift @ARGV;
 open(Fout,">results/searchspace.status");
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
                                 $i=$#par;write STDOUT; $ii=$i;write Fout;
				  #check if parmin<=parmax
                          if ($parmin[$#par]>$parmax[$#par]) 
                            {print "ERROR searchspace reading parameterrange: parmin > parmax\n";
                             exit 1;
                             }
                          $line=~s|par\w+\s*\Q[\E||;                     
				 }
     } close Fin;
 } print "starting search of level $searchlevel ...\n";close Fout;
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
{
            $i=0;foreach(@par){$par[$i]=$parmin[$i]+($parmax[$i]-$parmin[$i])/2;++$i}
            ($sta)=sta(); #sta at input parameter position
       $staorigin=$sta;$i=0;$minimum=1;
       foreach(@par)
        {$dpar=($parmax[$i]-$parmin[$i])/2**($searchlevel+2);
         if (abs($dpar)>abs($parstp[$i]))
         {
          $par[$i]+=$dpar;
          ($sta)=sta();
          open(Fin,"results/searchspace.status");$line=<Fin>;
          if ($line=~/exit searchspace/){$sta=0;close Fin;}
          else
     {open(Fout,">results/searchspace.status");$ii=0;
     foreach(@par){write Fout;++$ii;}
     print Fout " ...  current sta=$sta\n";
     close Fout;}
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; system ("cp -f ".$file." ".$file.".min.".$searchlevel);
                                        system ("cp -f ".$file.".par ".$file.".parmin.$searchlevel");}
			}

          if ($sta<=$staorigin){$minimum=0;}
          $ii=0;foreach(@par){$dd=sprintf("%e ",$par[$ii]);print Foutlevel $dd;++$ii} print Foutlevel $sta."\n";
	   ++$pointcounter; 
          $par[$i]-=2*$dpar;
          ($sta)=sta();
          open(Fin,"results/searchspace.status");$line=<Fin>;
          if ($line=~/exit searchspace/){$sta=0;close Fin;}
          else
     {open(Fout,">results/searchspace.status");$ii=0;
     foreach(@par){write Fout;++$ii;}
     print Fout " ...  current sta=$sta\n";
     close Fout;}
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; system ("cp -f ".$file." ".$file.".min.".$searchlevel);
                                        system ("cp -f ".$file.".par ".$file.".parmin.$searchlevel");}
			}
          if ($sta<=$staorigin){$minimum=0;}
          $ii=0;foreach(@par){$dd=sprintf("%e ",$par[$ii]);print Foutlevel $dd;++$ii} print Foutlevel $sta."\n";
          $par[$i]+=$dpar;
	   ++$pointcounter;
          }
         last if ($sta==0);
         ++$i;
       }
       last if($pointcounter==0);
       if ($minimum==1&&$sta!=0)
        {#the input parameters are probably a (local) minimum - so save a file
          ($sta)=sta();
          open(Fin,"results/searchspace.status");$line=<Fin>;
          if ($line=~/exit searchspace/){$sta=0;close Fin;}
          else
     {open(Fout,">results/searchspace.status");$ii=0;
     foreach(@par){write Fout;++$ii;}
     print Fout " ...  current sta=$sta\n";
     close Fout;}
        last if ($sta==0);
       if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; system ("cp -f ".$file." ".$file.".min.".$searchlevel);
                                        system ("cp -f ".$file.".par ".$file.".parmin.$searchlevel");}
			}
          foreach (@ARGV){$file=$_; system ("cp ".$file.".par ".$file.".$searchlevel.$minnumber");}
	  ++$minnumber;
        }
}
else
{print "reading point from file searchspace.".($searchlevel-1)."\n";
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
         {$par[$i]+=$dpar;
          ($sta)=sta();
          open(Fin,"results/searchspace.status");$line=<Fin>;
          if ($line=~/exit searchspace/){$sta=0;close Fin;}
          else
     {open(Fout,">results/searchspace.status");$ii=0;
     foreach(@par){write Fout;++$ii;}
     print Fout " ...  current sta=$sta\n";
     close Fout;}
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; system ("cp -f ".$file." ".$file.".min.".$searchlevel);
                                        system ("cp -f ".$file.".par ".$file.".parmin.$searchlevel");}
			}

          if ($sta<=$staorigin){$minimum=0;}
          $ii=0;foreach(@par){$dd=sprintf("%e ",$par[$ii]);print Foutlevel $dd;++$ii} print Foutlevel $sta."\n";
	   ++$pointcounter; 
          $par[$i]-=2*$dpar;
          ($sta)=sta();
          open(Fin,"results/searchspace.status");$line=<Fin>;
          if ($line=~/exit searchspace/){$sta=0;close Fin;}
          else
     {open(Fout,">results/searchspace.status");$ii=0;
     foreach(@par){write Fout;++$ii;}
     print Fout " ...  current sta=$sta\n";
     close Fout;}
          last if ($sta==0);
          if ($sta<$stamin){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; system ("cp -f ".$file." ".$file.".min.".$searchlevel);
                                        system ("cp -f ".$file.".par ".$file.".parmin.$searchlevel");}
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
         last if($pointcounter==0);
       if ($minimum==1&&$sta!=0)
        {#the input parameters are probably a (local) minimum - so save a file
          ($sta)=sta();if ($sta>$staorigin*1.00000001){print "ERROR - calculation of sta gives another result than stored in file searchspace.".($searchlevel-1)."\n";exit 1;}
       if ($sta<$stamin&&$sta!=0){$stamin=$sta; foreach (@ARGV)
                                       {$file=$_; system ("cp -f ".$file." ".$file.".min.".$searchlevel);
                                        system ("cp -f ".$file.".par ".$file.".parmin.$searchlevel");}
			}
          foreach (@ARGV){$file=$_; system ("cp ".$file.".par ".$file.".$searchlevel.$minnumber");}
	  ++$minnumber;
        }
       }
    
  }
  close Fin1;
}
$hours=(time-$starttime)/3600;
$estimate=$hours*2*$#par;
print "$pointcounter points calculated in  $hours h.\n Time estimate for next level ".($searchlevel+1).": $estimate h\n";
print Foutlevel "#$pointcounter points calculated in  $hours h.\n# Time estimate for next level ".($searchlevel+1).": $estimate h\n";
close Foutlevel;
 foreach (@ARGV) {$file=$_; system ("cp ".$file.".parmin.$searchlevel  ".$file);}
     open(Fout,">results/searchspace.status");print Fout " ... searchspace stopped\n"; close Fout;

exit 0;

# END OF MAIN PROGRAM
#****************************************************************************** 


sub sta {local $SIG{INT}='IGNORE'; 
 #print "#write modified parameterset to files *\n";
 foreach (@ARGV)
 {$file=$_; open (Fin, $file.".par");open (Fout1, ">".$file);open (Fout2,">results/searchspace.par");
   while($line=<Fin>)
     {$modline=$line;
      if ($line=~/^.*par/) {#here write modified parameter set to line
                            while ($line=~/^.*par\w+\s*\Q[\E/) 
			          {my $i=0;#insert a parameter
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
#			          {my $i=0;#insert a function parameter
#				   foreach (@parnam)
#				    {$pnam=$_;
#				     if ($line=~/^.*function\s*\Q[\E\s*$pnam/)
#                                        {$line=~s|function\s*\Q[\E\s*$pnam\s*\Q]\E|$par[$i] |;
#					}
#				     ++$i;
#				    }
#				  }

                             while ($line=~/^.*?function\s*\Q[\E(.*?)\Q]\E/) 
			          {$expression=$1;
				   #substitute into the expression the paramter values
				   my $i=0;#insert a function parameter
				   foreach (@parnam)
				    {$pnam=$_;
				     while ($expression=~/.*$pnam/)
                                        {$expression=~s|$pnam|$par[$i]|;
					}
				     ++$i;
				    }
				    # calculate the expression by a little perl program
				    open (Foutcc, ">./results/ccccccc.ccc");
				    printf Foutcc "#!/usr/bin/perl\nprint ".$expression.";\n";
				    close Foutcc;
				    system "./results/ccccccc.ccc > ./results/cccccc1.ccc";
                            system "chmod 755 ./results/ccccccc.ccc";
				    open (Fincc,"./results/cccccc1.ccc");
				    $data=<Fincc>; close Fincc;
				    system "rm ./results/ccccccc.ccc ./results/cccccc1.ccc";
				    # $data contains now the result of the mathematical expression
                                    $line=~s|function\s*\Q[\E.*?\Q]\E|$data |;
				   }



                            } print Fout1 $line;print Fout2 $modline;
     } close Fin;close Fout1;close Fout2;
     unless (rename "searchspace.par",$file.".par")
     {die "\n error:perhaps (perl rename cannot cross filesystems) \n";}
 }
# print "#call routine calcsta to calculate standard deviation\n";
 system ("./calcsta > results/searchspace.sta");
 open (Fin,"results/searchspace.sta"); 
 while($line=<Fin>){
           if($line=~/^.*sta\s*=/) {($sta)=($line=~m|sta\s*=\s*([\d.eEdD\Q-\E\Q+\E]+)|);} 
                   }close Fin;
 system ("rm results/searchspace.sta");
 return $sta;
}  

    
