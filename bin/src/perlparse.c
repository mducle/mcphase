//#include<cstdio>
#include<stdio.h>
#include<stdlib.h>

#include "martin.h"
#include "myev.h"
#include "perlparse.h"

void PrintComplexMatrix(FILE * fout, char * name, ComplexMatrix & mat)
{fprintf(fout,"my $%s=new Math::Matrix (\n",name);
 for (int i=mat.Rlo();i<=mat.Rhi();++i){fprintf(fout,"[");
   for (int j=mat.Clo();j<=mat.Chi();++j){fprintf(fout,"%g %+g *i",real(mat(i,j)),imag(mat(i,j)));
                                     if(j<mat.Chi())fprintf(fout,",");}
                                   fprintf(fout,"]");if(i<mat.Rhi())fprintf(fout,",\n");
                                   }fprintf(fout,");\n");
}

int perlparse(char*sipffilename
              ,double ** numbers,char ** numbernames    // number 
              ,char **  strings,char ** stringnames    // string variables
              ,ComplexMatrix **  operators,char ** operatornames    // operators
              )
{FILE * fin,*fout; 
 char  command[MAXNOFCHARINLINE],perlfile[MAXNOFCHARINLINE],instr[MAXNOFCHARINLINE];

 fin=fopen(sipffilename,"rb");
 while(feof(fin)==0){if(fgets(instr,MAXNOFCHARINLINE,fin)!=NULL){if(strncmp(instr,"#!noperl",8)==0)
                     {fclose(fin);return myparse(sipffilename,numbers,numbernames,strings,stringnames,operators,operatornames);}
                    }}
  
 sprintf(perlfile,"results/_%s.pl",sipffilename);
 fout=fopen(perlfile,"w");
  // output header to perl file
  fprintf(fout,"#--- start of autogenerated prolog ---\n"
               "use lib $ENV{MCPHASE_DIR}.\"/bin\";"
               "use warnings;\n"
               "use Math::Complex;\n"
               "use Matrix;\n"
               "$Math::Matrix::Precision = 5;\n");
  // output initial values for variables to perl file
    for(int i=0;numbernames[i]!=NULL;++i)fprintf(fout,"my $%s=%g;\n",numbernames[i],(*numbers[i]));
    for(int i=0;stringnames[i]!=NULL;++i)fprintf(fout,"my $%s=\"%s\";\n",stringnames[i],strings[i]);
    for(int i=0;operatornames[i]!=NULL;++i)PrintComplexMatrix(fout,operatornames[i],(*operators[i]));

 fprintf(fout,"#---- end autogenerated prolog ----\n\n\n");
  // output sipffile to perl program
  fseek(fin,0,SEEK_SET);
  while(feof(fin)==0){if(fgets(instr,MAXNOFCHARINLINE,fin)!=NULL)fprintf(fout,"%s",instr);}
  fclose(fin);


  // output printing statements to perl program  
  fprintf(fout,"\n\n\n#--- start of autogenerated epilog ---\n");
  fprintf(fout,"open Fout, \">results/_%s.pl.out\";\n",sipffilename);
  for(int i=0;numbernames[i]!=NULL;++i)fprintf(fout,"print Fout \"%s=\".$%s.\"\\n\";\n",numbernames[i],numbernames[i]);
  for(int i=0;stringnames[i]!=NULL;++i)fprintf(fout,"print Fout \"%s=\".$%s.\"\\n\";\n",stringnames[i],stringnames[i]);
  for(int i=0;operatornames[i]!=NULL;++i)fprintf(fout,"print Fout \"%s=\\n\",$%s->as_blocks();\n",operatornames[i],operatornames[i]);
  fprintf(fout,"close Fout;\n");
  fprintf(fout,"#--- end of autogenerated epilog ---\n");

 fclose(fout);

// system call to perl
 sprintf(command,"perl results/_%s.pl\n",sipffilename);
  
if(system(command)){fprintf(stderr,"Error parsing sipffile through perl\n");return false; }

// read perl output and fill variables with values
 double dummy;
 sprintf(perlfile,"results/_%s.pl.out",sipffilename);
  fin=fopen(perlfile,"rb");
  while(feof(fin)==0){
   if(fgets(instr,MAXNOFCHARINLINE,fin)!=NULL)
   {for(int i=0;numbernames[i]!=NULL;++i)extract(instr,numbernames[i],(*numbers[i]));
    for(int i=0;stringnames[i]!=NULL;++i)extract(instr,stringnames[i],strings[i],(size_t)MAXNOFCHARINLINE);
    for(int i=0;operatornames[i]!=NULL;++i)if(extract(instr,operatornames[i],dummy)==0)
                                            {if(myReadComplexMatrix (fin, (*operators[i]))==false)
                                              {fprintf(stderr,"Error parsing sipffile through perl - reading matrix from output\n");return false; }
                                             }
    }
                     }
 fclose(fin); 

return true;
}


// this substitutes the perl parsing in order to avoid big memory problems
// however the syntax of the commands to parse must be much simpler ...
int myparse(char*sipffilename
              ,double ** numbers,char ** numbernames    // number 
              ,char **  strings,char ** stringnames    // string variables
              ,ComplexMatrix **  operators,char ** operatornames    // operators
              ,int **sq2, double **csn, char **statements, int *iocs, int* oes)
{FILE * fin; 
 char  instr[MAXNOFCHARINLINE];
 printf("# using noperl parsing\n");
 fin=fopen(sipffilename,"rb");

 char line[MAXNOFCHARINLINE], *t0, *t1, *t2, *tokvar, *tokeq, *tokop, statement[MAXNOFCHARINLINE], *varpos[99], *oppos[99];
 char operr[]="Error parsing statement: %s\nOnly up to 99 %s allowed.\n";
 char sterr[]="Error. Inconsistent parsing of statement: %s\n";
 double constval[99];
 int idvar[99], idop[99], vc, oc, cc, seq[999], seq2[999], iline=0;

 while(feof(fin)==0)
 if(fgets(instr,MAXNOFCHARINLINE,fin)!=NULL && instr[strspn(instr," \t")]!='#' && strlen(instr)>3)
 {
    int lc=0, ops=0, eqs=0, ndl=0, var=0;

    memset(line,0,MAXNOFCHARINLINE-1);
    // First run, remove all spaces, counts number of operators, variables and delimiters
    t1=instr; t2 = strstr(instr,"\n");
    while(t1<t2) { 
       if( *t1!=' ' && *t1!='\t' )             line[lc++]=*t1;
       if( *t1=='+' || *t1=='-' || *t1=='*' )  ops++; else
       if( *t1=='$' )                          var++; else
       if( *t1=='=' )                          eqs++; else
       if( *t1==';' )                          ndl++;
       t1++;
    } line[lc]='\0';

    // Some sanity checks
    if(ndl==0)   { fprintf(stderr,"Error parsing line: %s\nEach statement must end with a ';'\n",instr); return false; }
    if(eqs!=ndl) { fprintf(stderr,"Error parsing line: %s\nEach statement should have a '=' and end with ';'\n",instr); return false; }
    if(var<eqs)  { fprintf(stderr,"Error parsing line: %s\nEach statement must involve a variable beginning with '$'\n",instr); return false; }

    // Splits the line into statements, delimited by ';' - and parses them
    t1 = line;
    for(int ist=0; ist<ndl; ist++)
    {
       vc=0; t2=strstr(t1,";"); tokeq=strstr(t1,"="); tokvar=strstr(t1,"$"); memset(statement,0,98); strncpy(statement,t1,(t2-t1)+1);
       if(tokeq ==NULL||tokeq >t2)    { fprintf(stderr,"Error parsing statement: %s\nStatement must have an '=' sign.\n",statement); return false; }
       if(tokvar==NULL||tokvar>tokeq) { fprintf(stderr,"Error parsing statement: %s\nYou must specify a variable name before the '='.\n",statement); return false; }
       while(tokvar!=NULL && tokvar<t2)
       {
          unsigned int opnamelen=0;
          for(int iop=0; operatornames[iop]!=NULL; iop++) 
             if(strncmp(tokvar+1,operatornames[iop],strlen(operatornames[iop]))==0) { 
                if(strlen(operatornames[iop])>opnamelen) { idvar[vc]=iop; opnamelen=strlen(operatornames[iop]); } }
          varpos[vc++]=tokvar;
          tokvar=strchr(tokvar+1,'$');
          if(vc>0 && tokvar<tokeq && tokvar!=NULL) { 
             fprintf(stderr,"Error parsing statement: %s\nOnly one variable may be LHS of '=' in each statement\n",statement); return false; }
          if(vc>98) { fprintf(stderr,"Error parsing statement: %s\nOnly up to 99 variables allowed.\n",statement); return false; }
       }
       // The commands in each statement is stored in the sequence array seq[], where each element may take one of these values:
       //    0 - end of statement
       //    1 - variable (currently only operator names allowed) 
       //    2 - constant (double precision floating point number)
       //    3 - + addition operator
       //    4 - - subtraction operator
       //    5 - * left scaling operator (multiplies a matrix by a constant from the left, e.g. a.M)
       //    6 - * matrix multiplication operator
       //    7 - * right scaling operator (multiplies a matrix by a constant from the right, e.g. a.M)
       tokop=t1; oc=0; cc=0;
       while(tokop<t2)
       {
          if(*tokop=='+') { oppos[oc]=tokop; idop[oc++]=3; } else
          if(*tokop=='-') { oppos[oc]=tokop; idop[oc++]=4; } else
          if(*tokop=='*') { 
             oppos[oc]=tokop;
             bool isrmat = *(tokop+1)=='$';                                       // True if right hand side of operator is a matrix
             t0 = tokop-1; while(*t0!='='&&*t0!='+'&&*t0!='-'&&*t0!='*') t0--;    // Search backwards in string statement until a delimiter is seen.
             t0++;
             if(t0[strspn(t0," \t.0123456789")]=='*' && isrmat){ idop[oc]=5;      // Scalar * Matrix
                constval[cc++]=strtod(t0,NULL); } else
             if( strchr(t0,'$')<tokop                && isrmat)  idop[oc]=6; else // Matrix * Matrix
             if( strchr(t0,'$')<tokop                && !isrmat){idop[oc]=7;      // Matrix * Scalar
                constval[cc++]=strtod(tokop+1,NULL); }
             oc++;
          }
          tokop++;
          if(oc>98) { fprintf(stderr,operr,statement,"operators"); return false; }
          if(oc>98) { fprintf(stderr,operr,statement,"constants"); return false; }
       }
       // Goes through the list of variables and operators to construct the sequence
       int ioc=0, ivc=1, iseq=1, oe=0;
       seq[0]=1;                                   // Assume the sequence starts with the operator to be assigned.
       if(oc>0 && oppos[0]<tokeq) {                // Operator-Assignment type, e.g. +=, -=, *=
          if(*(oppos[0]+1)=='=') oe=1;
          else { fprintf(stderr,"Error parsing statement: %s\nOperator '%c' before '=' but not of form '%c='.\n",statement,*oppos[0],*oppos[0]); return false; }
       }
       while(ioc<oc && ivc<vc)
       {
          if(oppos[ioc]<varpos[ivc]) {
             if(idop[ioc]==5) seq[iseq++]=2;       // Scalar * Matrix - set scalar position first
             seq[iseq++]=idop[ioc++];
             if(idop[ioc-1]==7) seq[iseq++]=2;     // Matrix * Scalar - set scalar position after
          }
          else {
             seq[iseq++]=1; ivc++; }
       }
       while(ivc<vc) { seq[iseq++]=1; ivc++; }
       while(ioc<oc) { if(idop[ioc]==5) seq[iseq++]=2; seq[iseq++]=idop[ioc++]; if(idop[ioc-1]==7) seq[iseq++]=2; }
       seq[iseq]=0;
//     printf("%s\t[",statement); for(int ii=0; ii<=iseq; ii++) printf("%i ",seq[ii]); printf("]\n");
       t1=t2+1;

       // Rearranges the sequence in terms of operations and variables: [operation,leftside,rightside,...]
       int i0=1, i1=1, pmflag=0; ioc=0; ivc=0;
       seq2[0] = idvar[ivc++]; // Assume seq[0] is the id of the operator to be assigned
       if(oe&&seq[1]>2) { seq2[i1++]=seq[1]; i0++; }                    // First operation is assignment-operator, +=,-=,*=
       if(seq[i0]==1&&seq[i0+1]<5) { seq2[i1++]=-idvar[ivc++]; i0++; }  // First operation is simple assignment
       while(seq[i0]!=0) {
          if(seq[i0]<3) i0++; else {
             switch(seq[i0]) {
                case 3:                        // + operation, 
                   // if next in sequence is a matrix and following operation not matrix*scalar, then assign,
                   //     else set flag to assign the product of the multiplication after next operation.
                   if(seq[i0+1]==1&&seq[i0+2]!=7) { seq2[i1++]=3; seq2[i1++]=-idvar[ivc++]; } else pmflag=3; i0++;
                   break;
                case 4:                        // - operation
                   if(seq[i0+1]==1&&seq[i0+2]!=7) { seq2[i1++]=4; seq2[i1++]=-idvar[ivc++]; } else pmflag=4; i0++;
                   break;
                case 5:                        // Scalar * Matrix, put matrix id into seq2 array
                   seq2[i1++]=seq[i0];
                   if(seq[i0-1]==2) seq2[i1++]=2;             else { fprintf(stderr,sterr,statement); return false; }
                   if(seq[i0+1]==1) seq2[i1++]=-idvar[ivc++]; else { fprintf(stderr,sterr,statement); return false; }
                   i0+=2; if(pmflag!=0) { seq2[i1++]=pmflag; pmflag=0; } 
                   break;
                case 6:                        // Matrix * Matrix, put both matrix id into seq2 array
                   seq2[i1++]=seq[i0];
                   if(seq[i0-1]==1) seq2[i1++]=-idvar[ivc++]; else { fprintf(stderr,sterr,statement); return false; }
                   if(seq[i0+1]==1) seq2[i1++]=-idvar[ivc++]; else { fprintf(stderr,sterr,statement); return false; }
                   i0+=2; if(pmflag!=0) { seq2[i1++]=pmflag; pmflag=0; } 
                   break;
                case 7:                        // Matrix * Scalar, put matrix id into seq2 array
                   seq2[i1++]=seq[i0];
                   if(seq[i0+1]==2) seq2[i1++]=2;             else { fprintf(stderr,sterr,statement); return false; }
                   if(seq[i0-1]==1) seq2[i1++]=-idvar[ivc++]; else { fprintf(stderr,sterr,statement); return false; }
                   i0+=2; if(pmflag!=0) { seq2[i1++]=pmflag; pmflag=0; } 
                   break;
             }
          }
       }
       seq2[i1]=0;
//     printf("%s\t[[",statement); for(int ii=0; ii<=iseq; ii++) printf("%i ",seq2[ii]); printf("]]\n");

       if(sq2!=NULL && csn!=NULL && statements!=NULL && iocs!=NULL && oes!=NULL)  // Output only the sequence and return
       {
          i0=0; int *sqv = new int[999]; while(seq2[i0]!=0&&i0<999) { sqv[i0]=seq2[i0]; i0++; } sq2[iline]=sqv;
          double *vcs = new double[99]; for(i0=0; i0<99; i0++) vcs[i0]=constval[i0]; csn[iline]=vcs;
          char *stam = new char[MAXNOFCHARINLINE]; strcpy(stam,statement); statements[iline]=stam;
          iocs[iline]=ioc; oes[iline++]=oe;
          if(iline>999) { fprintf(stderr,"noperl: Sorry I can only handle up to 999 statements for now.\n"); return false; }
       }
       else 
          if(myparse_execute(operators, operatornames, seq2, constval, statement, ioc, oe)==false) return false;
    }
 }
/*
  while(feof(fin)==0)
  {if(fgets(instr,MAXNOFCHARINLINE,fin)!=NULL)
   {if(instr[strspn(instr," \t")]!='#'&&strlen(instr)>3){
      // do something it is a command line
        // get first variable index
        char *token,*t1;int i1=-1,i0=strspn(instr," \t");
        if(instr[i0]!='$')return false; // must start with $ sign
t1=strstr(instr,";"); t1++; while(*t1!='\0') { *t1=' '; t1++; } printf("%s\t//\t",instr);
        for(int i=0;operatornames[i]!=NULL;++i)if ((token = strstr (instr, operatornames[i]))==instr+i0+1)
                                               {t1=token+strlen(operatornames[i]);i0=strspn(t1," \t");if(t1[i0]=='='||t1[i0]=='+')i1=i;
//printf("%s\t%s\n",operatornames[i],t1);
} 
        if (i1==-1)return false; // check if a parameter and a = after it is found - if not return false
        int pe=0;
//printf("i0=%i, %c %c\n",i0,t1[i0],t1[i0+1]);
        if(t1[i0]=='+'&&t1[i0+1]=='='){++t1;pe=1;}else if (t1[i0]=='+')return false;
        t1+=i0+1; // now i1 is the index of the ComplexMatrix to be assigned ...
       // now read the phrase after the = sign
       double cst=1.0;
       if(t1[strspn(instr," \t")]!='$'){//read cst to be multiplied
                      cst= strtod (t1, NULL);
//printf("%s\n",t1);
                      if((t1=strstr(t1,"*"))==NULL)return false;
                      ++t1;
                     }
       // check if $ ... i.e. variable name
//printf("%s\n",t1);
        i0=strspn(t1," \t");
        if(t1[i0]=='$')
          {// lets look which matrix is to be assigned
           int i2=-1;char * t2=t1;
           for(int i=0;operatornames[i]!=NULL;++i)
              if ((token = strstr (t1, operatornames[i]))==t1+i0+1) {
                 t2 = token+strlen(operatornames[i]); 
                 if(strchr(" \t;+-",t2[strspn(t2," \t")])!=NULL) i2=i; break;
//printf("%s\t%s\n",operatornames[i],t2);
} 
         //if (i2==-1) return false; // check if a parameter and a ; after it is found - if not return false
           if(t2[strspn(t2," \t")]!=';')
           {
              // Do the first operation
printf("op(%i)[%s] %s %f * op(%i)[%s] ",i1,operatornames[i1],(pe?"+=":"="),cst,i2,operatornames[i2]); 
              if(pe) if(cst==1.0) (*operators[i1])+=(*operators[i2]); else (*operators[i1])+=cst*(*operators[i2]);
              else   if(cst==1.0) (*operators[i1]) =(*operators[i2]); else (*operators[i1]) =cst*(*operators[i2]);
              // Loop through other operations
              while(t2[strspn(t2," \t")]!=';')
              {
                 cst=1.; int pm=1; i0=strspn(t2," \t"); 
                 if(t2[i0]=='+') pm=1; else if(t2[i0]=='-') pm=0; else return false; t2+=i0+1; i0=strspn(t2," \t"); 
                 if(t2[i0]!='$') { 
                    cst = strtod(t2,NULL); if((t2=strstr(t2,"*"))==NULL) return false; t2++; i0=strspn(t2," \t"); }
                 if(t2[i0]!='$') return false;
                 for(int i=0; operatornames[i]!=NULL; ++i)
                    if ( (token=strstr(t2, operatornames[i])) == t2+i0+1 ) { t2=token+strlen(operatornames[i]); i2=i; break; }
                 if(pm==1) // +
{ printf(" + %f * op(%i)[%s] ",cst,i2,operatornames[i2]);
                    if(cst==1) (*operators[i1]) += (*operators[i2]); else (*operators[i1]) += cst*(*operators[i2]);
}                else      // -
{ printf(" - %f * op(%i)[%s] ",cst,i2,operatornames[i2]);
                    if(cst==1) (*operators[i1]) -= (*operators[i2]); else (*operators[i1]) -= cst*(*operators[i2]);
}             }
printf("\n");
           } else { 
           if(pe){
           //  += assignment
printf("op(%i)[%s] += %f * op(%i)[%s]\n",i1,operatornames[i1],cst,i2,operatornames[i2]);
            if(cst==1.0){(*operators[i1])+=(*operators[i2]);//printf("%i=%i\n",i1,i2);
                       } else {(*operators[i1])+=cst*(*operators[i2]);}
           }else{
           // = asssign
printf("op(%i)[%s] = %f * op(%i)[%s]\n",i1,operatornames[i1],cst,i2,operatornames[i2]);
           if(cst==1.0){(*operators[i1])=(*operators[i2]);//printf("%i=%i\n",i1,i2);
                       } else {(*operators[i1])=cst*(*operators[i2]);}
           }

         }} else{return false;} // constants other than matrix names are not yet implemented
   }}
  }
*/
 fclose(fin);
 if(sq2!=NULL && csn!=NULL && statements!=NULL && iocs!=NULL && oes!=NULL) return iline;
 return true;
}

int myparse_execute(ComplexMatrix **operators, char **operatornames, int *seq2, double *constval, char *statement, int ioc, int oe)
{
    char *t0; //line[MAXNOFCHARINLINE], *t0, *t1, *t2, *tokvar, *tokeq, *tokop, statement[MAXNOFCHARINLINE], *varpos[99], *oppos[99];
    char dummystr[MAXNOFCHARINLINE];
    char sterr[]="Error. Inconsistent parsing of statement: %s\n";
    int i0=0;
    ComplexMatrix dummy, dummy2;

       // Prints out what we think the input should be for debugging purposes
       printf("%s\t==>\t",statement);

       // Now do the matrix manipulations, do the first operation outside the loop
       int is5=0;
       t0=statement+strlen(operatornames[seq2[0]])+1;
       if(*(t0+1)=='=') {                          // Operator-assignment [+=,-=,*=]
          if(!oe||seq2[1]<3) { fprintf(stderr,sterr,statement); return false; }
          if(seq2[2]<0) {
                dummy = (*operators[-seq2[2]]); i0=3; } else {              // Simple assignment
          switch(seq2[2]) {
             case 5:          sprintf(dummystr," %f * op(%i)[%s] ",constval[ioc],-seq2[4],operatornames[-seq2[4]]); is5=1;
             case 7: if(!is5){sprintf(dummystr," op(%i)[%s] * %f ",-seq2[4],operatornames[-seq2[4]],constval[ioc]);}is5=0;
                dummy =    constval[ioc++]     * (*operators[-seq2[4]]);    // Scalar*Matrix or Matrix*Scalar
                break;
             case 6: sprintf(dummystr," op(%i)[%s] * op(%i)[%s] ",-seq2[3],operatornames[-seq2[3]],-seq2[4],operatornames[-seq2[4]]);
                dummy = (*operators[-seq2[3]]) * (*operators[-seq2[4]]);    // Matrix*Matrix
                break;
             default: fprintf(stderr,"Error carrying out matrix operations in statement: %s\n",statement); return false; 
          } i0 = 5; }
          switch(*t0) {
             case '+': if(seq2[1]!=3) { fprintf(stderr,sterr,statement); return false; }
                (*operators[seq2[0]]) += dummy; printf(" op(%i)[%s] += %s ",seq2[0],operatornames[seq2[0]],dummystr); break;
             case '-': if(seq2[1]!=4) { fprintf(stderr,sterr,statement); return false; }
                (*operators[seq2[0]]) -= dummy; printf(" op(%i)[%s] -= %s ",seq2[0],operatornames[seq2[0]],dummystr); break;
             case '*': if(seq2[1]!=7) { fprintf(stderr,sterr,statement); return false; }
                                                printf(" op(%i)[%s] *= %s ",seq2[0],operatornames[seq2[0]],dummystr); break;
                dummy2 = (*operators[seq2[0]]) * dummy; (*operators[seq2[0]]) = dummy2;
          }
       }
       else {                                      // Assignment
          if(seq2[1]<0) {     sprintf(dummystr," op(%i)[%s] ",-seq2[1],operatornames[-seq2[1]]);
                (*operators[seq2[0]]) = (*operators[-seq2[1]]); i0=2; } else {              // Simple assignment
          switch(seq2[1]) {
             case 5:          sprintf(dummystr," %f * op(%i)[%s] ",constval[ioc],-seq2[3],operatornames[-seq2[3]]); is5=1;
             case 7: if(!is5){sprintf(dummystr," op(%i)[%s] * %f ",-seq2[3],operatornames[-seq2[3]],constval[ioc]);}is5=0;
                (*operators[seq2[0]]) =    constval[ioc++]     * (*operators[-seq2[3]]);    // Scalar*Matrix or Matrix*Scalar
                break;
             case 6: sprintf(dummystr," op(%i)[%s] * op(%i)[%s] ",-seq2[2],operatornames[-seq2[2]],-seq2[3],operatornames[-seq2[3]]);
                (*operators[seq2[0]]) = (*operators[-seq2[2]]) * (*operators[-seq2[3]]);    // Matrix*Matrix
                break;
             default: fprintf(stderr,"Error carrying out matrix operations in statement: %s\n",statement); return false; 
          } i0 = 4; }
          printf(" op(%i)[%s] = %s ",seq2[0],operatornames[seq2[0]],dummystr);
       }

       int pm;
       while(seq2[i0]!=0) {
          switch(seq2[i0]) {
             case 3:
             case 4:          sprintf(dummystr," op(%i)[%s] ",-seq2[i0+1],operatornames[-seq2[i0+1]]);
                dummy = (*operators[-seq2[i0+1]]); 
                pm = seq2[i0];   i0 += 2;
                break;
             case 5:          sprintf(dummystr," %f * op(%i)[%s] ",constval[ioc],-seq2[i0+2],operatornames[-seq2[i0+2]]); is5=1;
             case 7: if(!is5){sprintf(dummystr," op(%i)[%s] * %f ",-seq2[i0+2],operatornames[-seq2[i0+2]],constval[ioc]);}is5=0;
                dummy =      constval[ioc++]      * (*operators[-seq2[i0+2]]); // Scalar*Matrix or Matrix*Scalar
                pm = seq2[i0+3]; i0 += 4;
                break;
             case 6: sprintf(dummystr," op(%i)[%s] * op(%i)[%s] ",-seq2[i0+1],operatornames[-seq2[i0+1]],-seq2[i0+2],operatornames[-seq2[i0+2]]);
                dummy = (*operators[-seq2[i0+1]]) * (*operators[-seq2[i0+2]]); // Matrix*Matrix
                pm = seq2[i0+3]; i0 += 4;
                break;
             default: fprintf(stderr,"Error carrying out matrix operations in statement: %s\n",statement); return false; 
          }
          printf(" %c %s ",(pm==3?'+':'-'),dummystr);
          if(pm==3) (*operators[seq2[0]]) += dummy; else (*operators[seq2[0]]) -= dummy;
       }
       printf("\n");
  return true;
}

/* for test 
int main(int argc,char **argv)
{
//char sipffilename[]="test.sipf";
char numnam[]="alpha\0beta \0gamma\0";
double numbers[]={2.3,2.5,3.7};
char * numbernames[4];
for (int i=0;i<3;++i)numbernames[i]=numnam+6*i;numbernames[3]=NULL;

char * strings[3];char * stringnames[3];
char strnam[]="IONTYPE\0TITLE\0";   
stringnames[0]=strnam;stringnames[1]=strnam+8;stringnames[2]=NULL;
strings[0]=new char [MAXNOFCHARINLINE];
strings[1]=new char [MAXNOFCHARINLINE];
strcpy(strings[0],"Ce3p");
strcpy(strings[1],"Tomomi");

ComplexMatrix * operators[3];char * operatornames[3];
char opnam[]="I1 \0I2 \0";   
for (int i=0;i<2;++i)operatornames[i]=opnam+4*i;operatornames[2]=NULL;
operators[0]=new ComplexMatrix(1,3,1,3);(*operators[0])=0;
operators[1]=new ComplexMatrix(1,3,1,3);(*operators[1])=0;

if(perlparse(argv[1],numbers,numbernames,strings,stringnames,operators,operatornames)==false){printf("Error perl parsing sipf file\n");}

// print out modified number set
printf("%s = %g\n",numbernames[0],numbers[0]);
printf("%s = %g\n",numbernames[1],numbers[1]);
printf("%s = %g\n",numbernames[2],numbers[2]);

printf("%s = %s\n",stringnames[0],strings[0]);
printf("%s = %s\n",stringnames[1],strings[1]);

printf("%s =\n",operatornames[0]);myPrintComplexMatrix(stdout,(*operators[0]));
printf("%s =\n",operatornames[1]);myPrintComplexMatrix(stdout,(*operators[1]));
}

*/
