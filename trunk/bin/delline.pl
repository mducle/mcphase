#!/usr/bin/perl

#\begin{verbatim}



unless ($#ARGV >1) 

{print " program delline  used to delete lines from line1 to line2  \n";

 print " usage: comment line1 line2  *.*   \n  *.* .. filenname\n";

 exit 0;}

 

$col1=$ARGV[0];shift @ARGV;

$col2=$ARGV[0];shift @ARGV;



  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;

   open (Fout, ">range.out");

   $i=0;

   while($line=<Fin>)

     {++$i;

       if ($i<$col1||$i>$col2) {print Fout $line;}

       else{print "deleted line ".$i." :".$line;}

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