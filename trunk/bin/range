#!/usr/bin/perl
#\begin{verbatim}

unless ($#ARGV >2) 
{print " program range used to delete (comment out) data out of specified range from a file\n";
 print " usage: range  col min max *.*   \n (col=column,[min,max]..range \n *.* .. filenname\n";
 exit 0;}
 
$column=$ARGV[0];shift @ARGV;
$min=$ARGV[0];   shift @ARGV;
$max=$ARGV[0];   shift @ARGV;

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
#          print $numbers[$column-1]."\n";
          if ($min<$numbers[$column-1]&&$numbers[$column-1]<$max)
           	  {print Fout $line;}  
          else   
                  {print Fout "#:".$line;}  
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