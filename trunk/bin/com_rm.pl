#!/usr/bin/perl

# batch to remove every comment line 





unless ($#ARGV >=0) 

{print " program com_rm  used to remove all comment lines (starting with #) in a file\n";

 print " usage: com_rm *.*  \n *.* .. filenname\n";

 exit 0;}



#$command="#addc";foreach $d(@ARGV){$command.= " ".$d;}; $command.="\n";



  foreach (@ARGV)

  {

   $file=$_;

   unless (open (Fin, $file)){die "\n error:unable to open $file\n";}   
   print "<".$file;


   open (Fout, ">range.out");

   while($line=<Fin>)

     {

       if ($line=~/^\s*#/) {}

       else{print Fout $line;}

       

      }

      close Fin;

#      print Fout $command;

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

