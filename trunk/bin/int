#!/usr/bin/perl
#\begin{verbatim}

unless ($#ARGV >1) 
{print " program int  used to integrate columnx vs columny in data file, the result is displayed\n";
 print " integration is done point by point\n";
 print " usage: dif colx coly  *.*   \n colx=columnx, coly=columny \n *.* .. filenname\n";
 exit 0;}
 
$colx=$ARGV[0];shift @ARGV;
$coly=$ARGV[0];shift @ARGV;

  foreach (@ARGV)
  {
   $file=$_;
   print "<".$file;
   open (Fin, $file);
   $integral=0;
   while($line=<Fin>)
     {
       if ($line=~/^\s*#/) {print Fout $line;}
       else{$line=~s/D/E/g;@numbers=split(" ",$line);
            if (!@numbers1){@numbers1=@numbers;}
                 unless(0==($numbers[$colx-1]-$numbers1[$colx-1]))
		  {@numout=@numbers;
                   $integral+=($numbers[$coly-1]+$numbers1[$coly-1])/2*($numbers[$colx-1]-$numbers1[$colx-1]);
		  } 
	    @numbers1=@numbers;
           }
      }
      close Fin;
   print " x=col".$colx." y=f(x)=col".$coly." integral=".$integral.">\n";
   }

#\end{verbatim} 