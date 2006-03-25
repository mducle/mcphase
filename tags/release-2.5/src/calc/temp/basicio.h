// File: basicio.h
// $Id: basicio.h,v 1.1 2001/10/08 15:00:22 xausw Exp xausw $
// $Log: basicio.h,v $
// Revision 1.1  2001/10/08 15:00:22  xausw
// Initial revision
//

#ifndef BASICIO_H
#define BASICIO_H

#include <stdio.h>

struct FILE2 { FILE * sto;
               FILE * ste;
	     };
	       
int pty_master(char **name);
int pty_slave(void);
int fd_popen(const char *path, char *const argv[] );

FILE *spopen(const char *mode, const char *path,...);
int spclose(FILE *stream);

struct FILE2 popen2(const char *path,...);
int pclose2(struct FILE2 stream);

int ExecWait(const char *path, char * const argv[], char * const env[]);
#endif
