/******************************************************************************\
 zoeTools.h - part of the ZOE library for genomic analysis
 
Copyright (C) Ian Korf 2002-2013.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

The MIT License (MIT) - opensource.org/licenses/MIT

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
