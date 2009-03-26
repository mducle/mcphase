#!/usr/bin/perl
#\begin{verbatim}

unless ($#ARGV >1) 
{print " program dif  used to differentiate d(columny)/d(columnx) in data file, the result goes to coly\n";
 print " differentiation is done point by point\n";
 print " usage: dif colx coly  *.*   \n colx=columnx, coly=columny \n *.* .. filenname\n";
 exit 0;}
 
$colx=$ARGV[0];shift @ARGV;
$coly=$ARGV[0];shift @ARGV;

  foreach (@ARGV)
  {
   $file=$_;
   print "<".$file;
   open (Fin, $file);
   open (Fout, ">range.out");
   while($line=<Fin>)
     {
       if ($line=~/^\s*#/) {print Fout $line;}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);
            if (!@numbers1){@numbers1=@numbers;}
                 unless(0==($numbers[$colx-1]-$numbers1[$colx-1]))
		  {@numout=@numbers;
		   $numout[$coly-1]=($numbers[$coly-1]-$numbers1[$coly-1])/($numbers[$colx-1]-$numbers1[$colx-1]);                       
		   $i=0;
		   foreach (@numout)
		   {unless($i==$coly-1){$numout[$i]=($numbers[$i]+$numbers1[$i])/2;}
                    print Fout $numout[$i]." ";++$i;}     
                    print Fout "\n";
		  } 
	    @numbers1=@numbers;
           }
      }
      close Fin;
      close Fout;
       unless (rename "range.out",$file)
          {unless(open (Fout, ">$file"))     
      {die "\n error:unable to write to $file\n";}
      open (Fin, "range.out");
      while($line=<Fin>){ print Fout $line;}
      close Fin;
      close Fout;
      system "del range.out"; 
     }

   print ">\n";
   }

#\end{verbatim} 