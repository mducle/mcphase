//File: basicio.cpp
//$Id: basicio.cpp,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $"$
//$Log: basicio.cpp,v $
//Revision 1.1  2001/10/08 15:00:22  xausw
//Initial revision
//
//Revision 1.2  1999/03/15 09:08:37  herbie
//*** empty log message ***
//

//Basic Version from C.Dusek
/* pty.c */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <unistd.h>
#include <signal.h>
#include <sys/ioctl.h>
#include <sys/param.h>
#include <sysexits.h>
#include <stdarg.h>
#include <sys/wait.h>

#include "stdfunc.h"
#include "basicio.h"

static char pty_name[] ="/dev/ptyXY";


int pty_master(char **name);
int pty_slave(void);
int fd_popen(...);

extern int errno;
// *****************************************************************
int pty_master(char **name)
/*
 *    Allocate master terminal.
 * Arguments:
 *    name: Unless NULL, ptr to slave terminal name is inserted.
 * Return Value:
 *    File descriptor of master terminal, -1 on failure.
 */
{
   int         i, master_fd;
   char        *ptr;
   struct stat statbuff;
   static char ptychar[] =  "pqrs";             /* X */
   static char hexdigit[] = "0123456789abcdef"; /* Y */

   for(ptr=ptychar;*ptr;ptr++) {
      pty_name[5] = 'p';
      pty_name[8] = *ptr;  /* X */
      pty_name[9] = '0';   /* Y */

      if(stat(pty_name,&statbuff)<0) break;

      for(i=0;i<16;i++) {
         pty_name[9]=hexdigit[i]; /* 0-15 -> 0-9a-f */
         if((master_fd=open(pty_name,O_RDWR))>=0) {
            int j;
            /* try opening slave */
            pty_name[5]='t';
            if((j=open(pty_name,O_RDWR))<0) pty_name[5]='p';
            else {
               close(j);
               if(name!=NULL) *name=pty_name;
               return master_fd; }}}}
   return -1;
}
// ******************************************************
int pty_slave(void)
/*
 *    Open slave terminal belonging to previously allocated master.
 * Return Value:
 *    File descriptor of slave terminal, -1 on failure.
 */
{
   pty_name[5]='t';
   return open(pty_name,O_RDWR); 
}
// ********************************************************
int fd_popen(const char *path,char *const argv[])
{ // returns a file dsecriptor to a pipe to communicate
  // with the process specified in the parameterlist
  // argv must be terminated wit NULL;

 int filedes[2];   // two loneley pipes to appease the god of concurrency ...
 //int one=1;

 int fd=pty_master(NULL);
// ioctl(fd,FIONBIO,&one);

 struct sigaction act;
 sigfillset(&act.sa_mask);
 act.sa_handler=SIG_IGN;
 act.sa_flags=SA_RESTART;
 sigaction(SIGCLD,&act,NULL);

 if(pipe(filedes)) exit(1);

 int pid=fork();

 if(pid==-1){PRINT_DEBUG("Uuups, can not create child\n") 
             exit(EXIT_FAILURE);
             }
 if(pid){
    int dummy;
    close(filedes[1]);
    read(filedes[0],&dummy,sizeof(dummy));
    close(filedes[0]);

    return fd;  //parent
 }

 //child
 else
   {close (fd);
    setsid();
    fd=open("/dev/tty",O_RDWR);
    if(fd>=0)
      {ioctl(fd,TIOCNOTTY,NULL);
       close(fd);
      }
    for(fd=0;fd<NOFILE;fd++) if(fd!=filedes[1]) close(fd);

    fd=pty_slave();
    if(fd!=0) dup2(fd,0);
    if(fd!=1) dup2(fd,1);
    if(fd!=2) dup2(fd,2);
    if(fd>2) close(fd);

    fd=open("/dev/tty",O_RDWR);
    close(fd);

    { int dummy=1;
      write(filedes[1],&dummy,sizeof(dummy));
      close(filedes[1]); }

    execvp(path,argv);
    PRINT_DEBUG("Error executing %s\n",path)
    exit(EXIT_FAILURE);
   } 
}
// ******************************************************
int ExecWait(const char *path, char * const argv[], char * const env[])
{ 
  // Starts subprocess specified in the parameterlist
  // Waits until child terminates
  // argv & env must be terminated with NULL;
  // returns status of child process

 struct sigaction act;
 sigfillset(&act.sa_mask);
 act.sa_handler=SIG_IGN;
 act.sa_flags=SA_RESTART;
 sigaction(SIGCLD,&act,NULL);


 pid_t pid=fork();

 if(pid==-1){PRINT_DEBUG("Uuups, can not create child\n") 
             exit(EXIT_FAILURE);
             }
 if(pid)// parent
   {
   int r;
    pid_t pi=wait(&r);
    if(pi!=pid)PRINT_DEBUG("Uuups, programmers error when executing %s\n",path)
    if(!WIFEXITED(r))PRINT_DEBUG("Error executing %s\n",path)
    return r;    

   }

 //child
 else
   {//  FILE *f=fopen("bug.tmp","w");
//      fprintf(f,"CHILD: %s arg[0]:%s\n",path,argv[0]);
//      while(argv[i])
//           {fprintf(f,"CHILD: Arg[%d]: %s\n",i,argv[i]);
//            i++;}
//       fclose(f);
     int i=execve(path,argv,env);
     PRINT_DEBUG("ERROR %d executing execve(%s)",i,path); 
     perror("Error: ");
     _exit(EX_OSERR); 
     return -1; //not reached
   } 

	
}
// ******************************************************
static int pids[FOPEN_MAX]={0};
FILE *spopen(const char *mode, const char *path,...)
/*
 *    Secure popen.
 * Arguments:
 *    file:  Path to executable
 *    mode:  "r|w/e"
 * Return Value:
 *    I/O Stream.
 */
{
   int filedes[2];
   int std;  
   int pid;

   if(pipe(filedes)) 
     {PRINT_DEBUG("ERORR at pipe():");
      exit(EXIT_FAILURE);
     } 
   /*   
    * Mode:    r/e       w
    * par: pipe[0]  pipe[1]
    * cld: pipe[1]  pipe[0]
    */

   while((pid=fork())<0) {
      PRINT_DEBUG("WARNING: fork():");
      sleep(1); }
   if(pid) 
     { /* parent */
      std=*mode=='w';
      /* mode:  r/e  w
       * std:   0    1
       */
      close(filedes[!std]); /* close unused descriptor */
//      fprintf(stdout,"PARENT: mode: %c, closing fd[%d]:%d\n",
//                      *mode,!std,filedes[!std]);
      pids[filedes[std]]=pid;
//      fprintf(stdout,"PARENT: mode: %c, std:%d fd[0]:%d, fd[1]:%d\n",
//                      *mode,std,filedes[0],filedes[1]);
//      fprintf(stdout,"PARENT: mode: %c, opening fd[%d]:%d\n",
//                      *mode,std,filedes[std]);
      return fdopen(filedes[std],*mode=='e'?"r":mode);
     }
   
   /* child */
   std=*mode!='w';
   /* mode:  r/e  w
    * std:   1    0
    */
   close(filedes[!std]); /* close unused descriptor */
//   fprintf(stdout,"CHILD: mode: %c, closing fd[%d]:%d\n",
//                      *mode,!std,filedes[!std]);

   if(filedes[std]!=(*mode=='e'?2:std))
     {/* filedes[std] could be already in right position after 
       * a previous call to dae_detatch()                      */
//      fprintf(stdout,"CHILD: fd[%d]:%d != %d\n",std,filedes[std],(*mode=='e'?2:std));
//      fprintf(stdout,"CHILD: closing %d\n",(*mode=='e'?2:std));
      close(*mode=='e'?2:std);

//      int id2;
      if(dup2(filedes[std],(*mode=='e'?2:std))==-1)
        {PRINT_DEBUG("dup2(%d,%d): failed ",filedes[std],*mode=='e'?2:std);
         _exit(EX_OSERR);
	}
//      fprintf(stdout,"CHILD: dup2(fd[%d]:%d,%d)=%d\n",
//                               std,filedes[std],(*mode=='e'?2:std),id2);
//      fprintf(stdout,"CHILD: closing fd[%d]:%d\n",std,filedes[std]);
      close(filedes[std]);
//      fprintf(stdout,"CHILD: fd[%d]:%d closed=%d\n",std,filedes[std],id2);
     }
      
//   fprintf(stdout,"CHILD: fds finished\n");
     char     **argv=NULL;
     int      i=0;
     char     *c;
     va_list  args;

     va_start(args,path);

      do
       {c=va_arg(args,char*);
        argv=(char **) realloc(argv,(++i)*sizeof(char*));
        argv[i-1]=c;
       }
      while(c!=NULL);
     va_end(args);
//    fprintf(stdout,"CHILD: Executing %s %s %s %s\n",
//                    path,argv[0],argv[1],argv[2]);

     i=execvp(path,argv);
     
     PRINT_DEBUG("ERROR %d executing execvp(%s)",i,path); 
     perror("Error: ");
     _exit(EX_OSERR); 
   return NULL; /* not reached */
}
// *******************************************************
int spclose(FILE *stream)
/*
 * secure pclose
 */
{
   int *pid;
   int stat;

   if(!*(pid=pids+fileno(stream))) return -1;
   fflush(stream);
   fclose(stream);
   errno=0;
   if(waitpid(*pid,&stat,0)!=*pid)
      PRINT_DEBUG("WARNING: waitpid(%d)",*pid);
   *pid=0;
   return stat;
}
// *********************************************************
struct FILE2 popen2(const char *path,...)
/*
 *    Secure popen for stdout and stderr.
 * Arguments:
 *    file:  Path to executable
 *
 * Return Value:
 *    2 I/O Streams.
 */
{
   int filedes_o[2];
   int filedes_e[2];
   int std;  
   int pid;

   if(pipe(filedes_o) || pipe(filedes_e) ) 
     {PRINT_DEBUG("ERORR at pipe():");
      exit(EXIT_FAILURE);
     } 
   /*   
    * Mode:    r/e       w
    * par: pipe[0]  pipe[1]
    * cld: pipe[1]  pipe[0]
    */

   while((pid=fork())<0) {
      PRINT_DEBUG("WARNING: fork():");
      sleep(1); }
   if(pid) 
     { /* parent */
       /* mode r/e
          std  0  */
      std=0;	  
      close(filedes_o[!std]); /* close unused descriptor */
      close(filedes_e[!std]); /* close unused descriptor */
      pids[filedes_o[std]]=pid;
      pids[filedes_e[std]]=pid;
//       fprintf(stdout,"PARENT: fd_o[0]:%d, fd_o[1]:%d\n",filedes_o[0],filedes_o[1]);
//       fprintf(stdout,"PARENT: fd_e[0]:%d, fd_e[1]:%d\n",filedes_e[0],filedes_e[1]);

      struct FILE2 f2={0,0};
//       fprintf(stdout," sto: %p  ste: %p\n",f2.sto,f2.ste);
      f2.sto=fdopen(filedes_o[std],"r");
      f2.ste=fdopen(filedes_e[std],"r");
//       fprintf(stdout,"opened sto: %p  ste: %p\n",f2.sto,f2.ste);
      return f2;
     }
   
   /* child */
   /* mode r/e
      std  1 */
   std=1;
   close(filedes_o[!std]);
   close(filedes_e[!std]);
   if(filedes_o[std]!=std)
     {/* filedes[1] could be already in right position after 
       * a previous call to dae_detatch()                      */
//       fprintf(stdout,"CHILD: fd_o[%d]:%d != %d\n",std,filedes_o[std],std);
//       fprintf(stdout,"CHILD: closing %d: \n",std);
      close(std);

      int id2;
      if((id2=dup2(filedes_o[std],std))==-1)
        {PRINT_DEBUG("dup2(%d,%d): failed ",filedes_o[std],std);
         _exit(EX_OSERR);
	}
//      fprintf(stdout,"CHILD: dup2(fd_o[%d]:%d,%d)=%d\n",std,filedes_o[std],std,id2);
      close(filedes_o[std]);
//      fprintf(stdout,"CHILD: fd_o[%d]:%d closed\n",std,filedes_o[std]);
     }

   if(filedes_e[std]!=2)
     {/* filedes[0] could be already in right position after 
       * a previous call to dae_detatch()                      */
//      fprintf(stdout,"CHILD: fd_e[%d]:%d != 2\n",std,filedes_e[std]);
//      fprintf(stdout,"CHILD: closing 2\n");
      close(2);

      int id2;
      if((id2=dup2(filedes_e[std],2))==-1)
        {PRINT_DEBUG("dup2(%d,2): failed ",filedes_e[std]);
         _exit(EX_OSERR);
	}
//      fprintf(stdout,"CHILD: dup2(fd_e[%d]:%d,2)=%d\n",std,filedes_e[std],id2);
//      fprintf(stdout,"CHILD: fd_e[%d]:%d closed\n",std,filedes_e[std]);
      close(filedes_e[std]);
     }
      
//   fprintf(stdout,"CHILD: fds finished\n");
     char     **argv=NULL;
     int      i=0;
     char     *c;
     va_list  args;

     va_start(args,path);

      do
       {c=va_arg(args,char*);
        argv=(char **) realloc(argv,(++i)*sizeof(char*));
        argv[i-1]=c;
       }
      while(c!=NULL);
     va_end(args);
      
     i=0;
//      FILE *f=fopen("bug.tmp","w");
//      fprintf(f,"CHILD: %s arg[0]:%s\n",path,argv[0]);
//      while(argv[i])
//           {fprintf(f,"CHILD: Arg[%d]: %s\n",i,argv[i]);
//            i++;}
//       fclose(f);	   
     i=execvp(path,argv);
     
     PRINT_DEBUG("ERROR %d executing execvp(%s)",i,path); 
     perror("Error: ");
     _exit(EX_OSERR); 
   struct FILE2 n={0,0};
   return n; /* not reached */
}
// *******************************************************
int pclose2(struct FILE2 stream)
/*
 * secure pclose
 */
{
   int *pid;
   int stat;

   if(!*(pid=pids+fileno(stream.sto))) return -1;
   fflush(stream.sto);
   fclose(stream.sto);

   if(!*(pid=pids+fileno(stream.ste))) return -1;
   fflush(stream.ste);
   fclose(stream.ste);
   errno=0;
   if(waitpid(*pid,&stat,0)!=*pid)
      PRINT_DEBUG("WARNING: waitpid(%d)",*pid);
   *pid=0;

   return stat;
}
// *********************************************************


