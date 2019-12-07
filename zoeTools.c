/******************************************************************************\
 zoeTools.c - part of the ZOE library for genomic analysis
 
Copyright (C) 2002-2005 Ian Korf

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

\******************************************************************************/

#ifndef ZOE_TOOLS_C
#define ZOE_TOOLS_C

#include "zoeTools.h"

/******************************************************************************\
 Library Information
\******************************************************************************/

static char zoeVersionNumber[] = "2006-07-28";

void zoeLibInfo (void) {
	zoeS(stderr, "ZOE library version %s\n", zoeVersionNumber);
}


/******************************************************************************\
 Program Name
\******************************************************************************/

static char PROGRAM_NAME[256] = "unnamed program";

void zoeSetProgramName (const char * string) {
	(void)strcpy(PROGRAM_NAME, string);
}

char * zoeGetProgramName (void) {
	return PROGRAM_NAME;
}

/******************************************************************************\
 Printing and Error Messages
\******************************************************************************/

void zoeS (FILE * stream, const char * fmt, ...) {
	va_list args;
	
	va_start(args, fmt);
	(void)vfprintf(stream, fmt, args);
	va_end(args);
	(void)fflush(stream);	
}

void zoeExit (const char * fmt, ...) {
	va_list args;
		
	(void)fflush(stdout);
	(void)fprintf(stderr, "ZOE ERROR (from %s): ", zoeGetProgramName());
	va_start(args, fmt);
	(void)vfprintf(stderr, fmt, args);
	va_end(args);
	(void)fprintf(stderr, "\n");
	zoeLibInfo();
	exit(2);
}


/******************************************************************************\
 Memory Tools
\******************************************************************************/

void * zoeMalloc (size_t size) {
	void * buffer;
	
	if ((buffer = malloc(size)) == NULL) zoeExit("zoeMalloc");
	return buffer;
}

void * zoeCalloc (size_t nobj, size_t size) {
	void * buffer;

	if ((buffer = calloc(nobj, size)) == NULL) zoeExit("zoeCalloc");
	return buffer;  
}

void * zoeRealloc (void * p, size_t size) {
	void * buffer;
	
	if ((buffer = realloc(p, size)) == NULL) zoeExit("zoeRealloc");
	return buffer;  
}

void zoeFree (void * p) {
	free(p);
	p = NULL;
}

/******************************************************************************\
 Text Vector
\******************************************************************************/

void zoeDeleteTVec (zoeTVec vec) {
	int i;
	
	if (vec == NULL) return;
	if (vec->elem) {
		for (i = 0; i < vec->size; i++) zoeFree(vec->elem[i]);
		zoeFree(vec->elem);
		vec->elem = NULL;
	}
	zoeFree(vec);
	vec = NULL;
}

zoeTVec zoeNewTVec (void) {
	zoeTVec vec = zoeMalloc(sizeof(struct zoeTVec));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void zoePushTVec (zoeTVec vec, const char * text) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zoeRealloc(vec->elem, vec->limit * sizeof(char *));
	}
	vec->elem[vec->size] = zoeMalloc(strlen(text) +1);
	(void)strcpy(vec->elem[vec->size], text);
	vec->last = vec->elem[vec->size];
	vec->size++;
}


/******************************************************************************\
 zoeVec
\******************************************************************************/

void zoeDeleteVec (zoeVec vec) {	
	if (vec == NULL) return;
	if (vec->elem) {
		zoeFree(vec->elem);
		vec->elem = NULL;
	}
	zoeFree(vec);
	vec = NULL;
}

zoeVec zoeNewVec (void) {
	zoeVec vec = zoeMalloc(sizeof(struct zoeVec));
	vec->size  = 0;
	vec->limit = 0;
	vec->elem  = NULL;
	return vec;
}

void zoePushVec (zoeVec vec, void * thing) {
	if (vec->limit == vec->size) {
		if (vec->limit == 0) vec->limit  = 1;
		else                 vec->limit *= 2;
		vec->elem = zoeRealloc(vec->elem, vec->limit * sizeof(void *));
	}
	vec->elem[vec->size] = thing;
	vec->last = vec->elem[vec->size];
	vec->size++;
}

int zoeTcmp (const void * a, const void * b) {
	return strcmp( *(char **)a, *(char **)b );
}


/******************************************************************************\
 Generic Hash

The hash function is my own creature. I used the multiplication method as
described  in Cormen et al. but added a 7 periodic component so that for
example, ATG and TGA don't hash to the same index. Because DNA might be hashed,
I thought it would be a good idea to avoid some multipe of 3 and because number
systems are 2 or 10 based, I stayed away from those too. The 7 values chosen
were kind of arbitrary. I have tested this on some rather large text files, and
the hash function separation is quite good even after several levels of
re-hashing. Performance is good too, about 6x faster than the C++ map and 3
times faster than a Perl hash.

\******************************************************************************/

static double zoeHASH_MULTIPLIER[7] = {
	3.1415926536, /* PI */
	2.7182818285, /* e */
	1.6180339887, /* golden mean */
	1.7320508076, /* square root of 3 */
	2.2360679775, /* square root of 5 */
	2.6457513111, /* square root of 7 */
	3.3166247904, /* square root of 11 */
};

static float zoeMAX_HASH_DEPTH = 2.0;         /* The hash will remain between */

static int zoeHashLevelToSlots (int level) {  /* half full and twice-filled */	
	return pow(4, level);                     /* with these values */
}

static int zoeHashFunc (const zoeHash hash, const char * key) {
	size_t i;
	double sum;
	
	sum = 0;
	for (i = 0; i < strlen(key); i++) {
		sum += key[i] * zoeHASH_MULTIPLIER[i % 7];
	}
	
	return (int) (hash->slots * (sum - floor(sum)));
}

static void zoeExpandHash (zoeHash hash) {
	int      i, j;
	char   * key = NULL;
	void   * val = NULL;
	int      oldslots = hash->slots;
	zoeVec * oldkey = hash->key;
	zoeVec * oldval = hash->val;
	zoeVec   kvec;
	zoeVec   vvec;
	zoeTVec  keys;
		
	/* create the new hash */
	hash->level = hash->level +1;
	hash->slots = zoeHashLevelToSlots(hash->level);
	hash->key   = zoeMalloc(hash->slots * sizeof(struct zoeVec));
	hash->val   = zoeMalloc(hash->slots * sizeof(struct zoeVec));
	for (i = 0; i < hash->slots; i++) {
		hash->key[i] = zoeNewVec();
		hash->val[i] = zoeNewVec();
	}
	
	/* brand new hash? */
	if (hash->keys->size == 0) return;

	keys = hash->keys;
	hash->keys = zoeNewTVec();
	
	/* transfer old stuff to new hash */
	for (i = 0; i < oldslots; i++) {
		kvec = oldkey[i];
		vvec = oldval[i];
		for (j = 0; j < kvec->size; j++) {
			key = kvec->elem[j];
			val = vvec->elem[j];
			zoeSetHash(hash, key, val);
		}
	}
	
	/* free old stuff */
	for (i = 0; i < oldslots; i++) {
		kvec = oldkey[i];
		vvec = oldval[i];
		zoeDeleteVec(kvec);
		zoeDeleteVec(vvec);
	}
	zoeFree(oldkey);
	zoeFree(oldval);
	zoeDeleteTVec(keys);
	
}

void zoeDeleteHash (zoeHash hash) {
	int i;
	
	if (hash == NULL) return;
	
	for (i = 0; i < hash->slots; i++) {
		if (hash->key[i]) {
			zoeDeleteVec(hash->key[i]);
			hash->key[i] = NULL;
		}
		if (hash->val[i]) {
			zoeDeleteVec(hash->val[i]);
			hash->val[i] = NULL;
		}
	}
	zoeDeleteTVec(hash->keys);
	hash->keys = NULL;
	zoeDeleteVec(hash->vals);
	hash->vals = NULL;
	zoeFree(hash->key);
	hash->key = NULL;
	zoeFree(hash->val);
	hash->val = NULL;
	zoeFree(hash);
	hash = NULL;
}

zoeHash zoeNewHash (void) {
	zoeHash hash = zoeMalloc(sizeof(struct zoeHash));
	hash->level = 0;
	hash->slots = 0;
	hash->keys  = zoeNewTVec();
	hash->vals  = zoeNewVec();
	hash->key   = NULL;
	hash->val   = NULL;
	zoeExpandHash(hash);
	return hash;
}

void * zoeGetHash (const zoeHash hash, const char * key) {
	int    i, index;
	char * string = NULL;

	index = zoeHashFunc(hash, key);
	/* resolve collisions */
	for (i = 0; i < hash->key[index]->size; i++) {
		string = hash->key[index]->elem[i];
		if (strcmp(key, string) == 0) {
			return hash->val[index]->elem[i];
		}
	}
	return NULL; /* return is NULL if not found */
}

void zoeSetHash (zoeHash hash, const char * key, void * val) {
	int    i, index;
	char * string = NULL;
	int    new_key = 1;
	
	index = zoeHashFunc(hash, key);
	
	/* reassign unless new key */
	for (i = 0; i < hash->key[index]->size; i++) {
		string = hash->key[index]->elem[i];
		if (strcmp(key, string) == 0) {
			hash->val[index]->elem[i] = val;
			new_key = 0;
			return;
		}
	}
	
	if (new_key) {
		zoePushTVec(hash->keys, key);
		zoePushVec(hash->key[index], hash->keys->last);
		zoePushVec(hash->vals, val);
		zoePushVec(hash->val[index], hash->vals->last);
	}
	
	/* check if we have to expand the hash */
	if ((float)hash->keys->size / (float)hash->slots >= zoeMAX_HASH_DEPTH) {
		zoeExpandHash(hash);
	}
}

zoeTVec zoeKeysOfHash (const zoeHash hash) {
	int     i;
	zoeTVec vec = zoeNewTVec();
	
	for (i = 0; i < hash->keys->size; i++) zoePushTVec(vec, hash->keys->elem[i]);
	
	return vec;
}

zoeVec zoeValsOfHash (const zoeHash hash) {
	int    i;
	zoeVec vec = zoeNewVec();

	for (i = 0; i < hash->vals->size; i++) zoePushVec(vec, hash->vals->elem[i]);

	return vec;
}

void zoeStatHash (const zoeHash hash) {
	int i, max, min, total, count;

	max = 0;
	min = 1000000;
	total = 0;
	for (i = 0; i < hash->slots; i++) {
		count = hash->val[i]->size;
		total += count;
		if (count > max) max = count;
		if (count < min) min = count;
	}
	zoeS(stdout, "HashStats: level=%d slots=%d keys=%d min=%d max=%d ave=%f\n",
		hash->level, hash->slots, hash->keys->size, min, max,
		(float)total / (float)hash->slots);
}

#endif
