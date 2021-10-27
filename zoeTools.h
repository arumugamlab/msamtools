/******************************************************************************\
 zoeTools.h - part of the ZOE library for genomic analysis
 
 Copyright (C) 2002-2005 Ian Korf

\******************************************************************************/

#ifndef ZOE_TOOLS_H
#define ZOE_TOOLS_H

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void   zoeLibInfo (void);
void   zoeSetProgramName (const char*);
char * zoeGetProgramName (void);

void zoeS (FILE *, const char *, ...);
void zoeExit  (const char *, ...);

void * zoeMalloc (size_t);
void * zoeCalloc (size_t, size_t);
void * zoeRealloc (void *, size_t);
void   zoeFree (void *);

struct zoeTVec  {
	char ** elem;
	int     size;
	int     limit;
	char  * last;
};
typedef struct zoeTVec * zoeTVec;
void    zoeDeleteTVec (zoeTVec);
zoeTVec zoeNewTVec (void);
void    zoePushTVec (zoeTVec, const char *);
int zoeTcmp (const void * a, const void * b);

struct zoeVec  {
	void ** elem;
	int     size;
	int     limit;
	void  * last;
};
typedef struct zoeVec * zoeVec;
void   zoeDeleteVec (zoeVec);
zoeVec zoeNewVec (void);
void   zoePushVec (zoeVec, void *);

struct zoeHash  {
	int      level;
	int      slots;
	zoeTVec  keys;
	zoeVec   vals;
	zoeVec * key;
	zoeVec * val;
};
typedef struct zoeHash * zoeHash;
void    zoeDeleteHash (zoeHash);
zoeHash zoeNewHash (void);
void    zoeSetHash (zoeHash, const char *, void *);
void *  zoeGetHash (const zoeHash, const char *);
zoeTVec zoeKeysOfHash (const zoeHash);
zoeVec  zoeValsOfHash (const zoeHash);
void    zoeStatHash (const zoeHash);

#endif
