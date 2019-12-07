#include "mBamVector.h"

/* Generic BAM/SAM utilities */

int32_t bam_highquality_differences_only(const bam1_t *b, uint8_t minqual) {

	/* iterators */

	int32_t  l;
	int32_t curr;

	/* MD */

	uint8_t  *mdz;
	char     *md;

	/* tokeniser for MD */

	ks_tokaux_t aux;
	char s[8];
	char *p;
	char *prev;

	/* Get quality */

	uint8_t *q   = bam1_qual(b);

	/* get MD values */

	mdz = bam_aux_get(b, "MD");

	/* If there is no MD, then we can't proceed, so skip the read! */

	if (!mdz) return 0;

	md  = bam_aux2Z(mdz);

	prev = md;
	curr = 0;
	for (p = kstrtok(md, "^acgtnACGTN", &aux); p; p = kstrtok(0, 0, &aux)) {
		int32_t len;
		char *x;

		/* This loop covers positions where there are mismatches */

		for (x=prev; x<p; x++, curr++) {
			if (q[curr] < minqual) {
				return 0;
			}
		}

		/* If you come across a refdel sequence, skip it */
		if (p[0] == '^') {
			do {
				p = kstrtok(0, 0, &aux);
			} while (index("acgtnACGTN", p[0]) != NULL);
		}

		strncpy(s, p, aux.p - p);
		s[aux.p-p] = '\0';
		len = atoi(s);
		for (l=0; l<len; l++,curr++);
		prev = (char*) aux.p;
	}

	return 1;
}

int32_t bam_cigar2alnlen(const bam1_core_t *c, const uint32_t *cigar) {
	uint32_t k;
	int32_t length = 0;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		if (!(op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP || op == BAM_CREF_SKIP || op == BAM_CPAD)) {
			length += cigar[k] >> BAM_CIGAR_SHIFT;
		}
	}
	return length;
}

/*
alen - alignment length,
qlen - length of the query read including clipped ends,
qclip - length of the clipped bases (from both sides)
*/
void bam_cigar2details(const bam1_core_t *c, const uint32_t *cigar, int32_t *alen, int32_t *qlen, int32_t *qclip) {
	uint32_t k;
	*alen = *qlen = *qclip = 0;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		int w  = cigar[k] >> BAM_CIGAR_SHIFT;
		if (op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP) {
			*qclip += w;
			*qlen  += w;
		} else if (!(op == BAM_CREF_SKIP || op == BAM_CPAD)) {
			*alen += w;
			if (op == BAM_CMATCH || op == BAM_CINS)
				*qlen += w;
		}
	}
}

void bam_get_summary(const bam1_t *b, mAlignmentSummary *summary) {

	bam1_core_t *c     = (bam1_core_t*) &b->core;
	uint32_t    *cigar = bam1_cigar(b);

	/* MD stuff */

	uint8_t  *mdz;
	char     *md;

	/* tokeniser for MD */

	ks_tokaux_t aux;
	char *p;

	int32_t alen = 0, qlen = 0, qclip = 0, match = 0, edit = 0;

	uint32_t k;

	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		int w  = cigar[k] >> BAM_CIGAR_SHIFT;
		switch (op) {
		/* Aligned Position */
			case BAM_CMATCH:
				match += w;
				qlen  += w;
				alen  += w;
				break;

		/* Insertion in reference */
			case BAM_CINS:
				qlen  += w;
				/* FALL-THROUGH */

		/* Deletion in reference */
			case BAM_CDEL:
				edit  += w;
				alen  += w;
				break;

		/* Clipped Read */
			case BAM_CHARD_CLIP:  
			case BAM_CSOFT_CLIP:
				qclip += w;
				qlen  += w;
				break;

		/* Not aligned to reference */
			case BAM_CREF_SKIP:  
			case BAM_CPAD:
			default:
				break;
		}
	}

	/* get MD values */

	mdz = bam_aux_get(b, "MD");
	if (mdz) {
		md  = bam_aux2Z(mdz);

		/*****
		 * Remember: when you use kstrtok,
		 * p gives the next token,
		 * aux.p gives the end of this token,
		 * so aux.p-p gives length of this token
		 *****/

		for (p = kstrtok(md, "^0123456789", &aux); p; p = kstrtok(0, 0, &aux)) {
			char *x;
			if (p > md && p[-1] != '^')
				for (x=p; x<aux.p; x++)
					edit++;
		}
		match -= edit;
	}

	summary->match = match;
	summary->edit = edit;
	summary->query_length = qlen;
	summary->query_clip = qclip;
	summary->length = alen;

/*
	fprintf(stdout, "match=%d, edit=%d, qlen=%d, qclip=%d\n", match, edit, qlen, qclip);
*/

	return;
}

void bam_get_extended_summary(const bam1_t *b, mAlignmentSummary *summary) {

	bam1_core_t *c     = (bam1_core_t*) &b->core;
	uint32_t    *cigar = bam1_cigar(b);

	/* MD stuff */

	uint8_t  *mdz;
	char     *md;

	/* tokeniser for MD */

	ks_tokaux_t aux;
	char *p;

	int32_t alen = 0, qlen = 0, qclip = 0, match = 0, mismatch = 0, gapopen = 0, gapextend = 0;

	uint32_t k;

	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		int w  = cigar[k] >> BAM_CIGAR_SHIFT;
		switch (op) {
			/* Aligned Position */
			case BAM_CMATCH:
				match += w;
				qlen  += w;
				alen  += w;
				break;

			/* WATCH THE RUN THROUGH - THERE'S NO BREAK HERE FOR EFFICIENCY */
			/* Insertion in reference */
			case BAM_CINS:
				qlen  += w;
			/* Deletion in reference */
			case BAM_CDEL:
				gapopen++;
				gapextend += (w-1);
				alen  += w;
				break;

		/* Clipped Read */
			case BAM_CHARD_CLIP:  
			case BAM_CSOFT_CLIP:
				qclip += w;
				qlen  += w;
				break;

		/* Not aligned to reference */
			case BAM_CREF_SKIP:  
			case BAM_CPAD:
			default:
				break;
		}
	}

	/* get MD values */

	mdz = bam_aux_get(b, "MD");
	if (mdz) {
		md  = bam_aux2Z(mdz);

		/*****
		 * Remember: when you use kstrtok,
		 * p gives the next token,
		 * aux.p gives the end of this token,
		 * so aux.p-p gives length of this token
		 *****/

		for (p = kstrtok(md, "^0123456789", &aux); p; p = kstrtok(0, 0, &aux)) {
			char *x;
			if (p > md && p[-1] != '^')
				for (x=p; x<aux.p; x++)
					mismatch++;
		}
		match -= mismatch;
	}

	summary->match        = match;
	summary->mismatch     = mismatch;
	summary->gapopen      = gapopen;
	summary->gapextend    = gapextend;
	summary->query_length = qlen;
	summary->query_clip   = qclip;
	summary->length       = alen;
	summary->edit         = mismatch + qclip + gapopen + gapextend;

/*
	fprintf(stdout, "match=%d, mismatch=%d, gapopen=%d, gapextend=%d\n", match, mismatch, gapopen, gapextend);
*/

	return;
}

/* mBamVector */

void mSortBamVector(mBamVector *vec, int(*compar)(const void *, const void *)) {
	qsort(vec->elem, vec->size, sizeof(bam1_t*), compar);
}

void mInitBamVector(mBamVector *vec, int limit) {
	vec->limit = (limit>0)?limit:1;
	/* Don't want to handle empty vectors. This is faster. Memory doesnt matter */
	vec->elem = (bam1_t**) mMalloc(vec->limit*sizeof(bam1_t*));
	vec->size = 0;
}

void mPushBamVector(mBamVector *vec, bam1_t *item) {
	if (vec->size == vec->limit) {
#ifdef DEBUG
		fprintf(stderr, "Expanding vector to %d\n", vec->limit *2);
#endif
		vec->limit *= 2;
		vec->elem = (bam1_t**) mRealloc(vec->elem, vec->limit*sizeof(bam1_t*));
	}
	vec->elem[vec->size] = item;
	vec->size++;
}

void mFreeBamVector(mBamVector *vec) {
	mFree(vec->elem);
}

void mEmptyBamVector(mBamVector *vec) {
	int i;
	for (i=0; i<vec->size; i++) {
		if (vec->elem[i] != NULL) {
			bam_destroy1(vec->elem[i]);
		}
	}
	vec->size = 0;
}

void mWriteBamVector(samfile_t *stream, mBamVector *bamvector) {
	int i;
	for (i=0; i<bamvector->size; i++) {
		bam1_t *b = bamvector->elem[i];
		if (b != NULL) {
			samwrite(stream, b);
		}
	}
}

/* mBamPool */

void mInitBamPool(mBamPool *pool, int limit) {
	int i;
	bam1_t **elem;
	pool->limit = (limit>0)?limit:1;
	pool->elem = (bam1_t**) mMalloc(pool->limit*sizeof(bam1_t*));
	pool->size = 0;
	elem = pool->elem;
	for (i=0; i<pool->limit; i++)
		elem[i] = bam_init1();
}

void mExpandBamPool(mBamPool *pool) {
	int i;
	int limit = pool->limit;
	bam1_t **elem = pool->elem;
	elem = (bam1_t**) mRealloc(elem, 2*limit*sizeof(bam1_t*));

	/* the following lines should NOT use pool->limit and limit interchangeably, this is the old limit */

	/* init the new bams */
	for (i=limit; i<2*limit; i++)
		elem[i] = bam_init1();

	/* Reset the limit */
	pool->limit *= 2;

	pool->elem  = elem;
}

void mReOriginateBamPool(mBamPool *pool) {
	if (pool->size == 0) {
		return;
	}
	bam_copy1(pool->elem[0], pool->elem[pool->size]);
	pool->size = 0;
}

bam1_t* mAdvanceBamPool(mBamPool *pool) {
	pool->size++;
	if (pool->size == pool->limit) {
		mExpandBamPool(pool);
#ifdef DEBUG
		fprintf(stderr, "BAM pool expanded to %d!\n", pool->limit);
#endif
	}
	return pool->elem[pool->size];
}

void mFreeBamPool(mBamPool *pool) {
	int i;
	for (i=0; i<pool->limit; i++)
		bam_destroy1(pool->elem[i]);
	mFree(pool->elem);
}

void mWriteBamPool(samfile_t *stream, mBamPool *pool) {
	int i;
	bam1_t **elem = pool->elem;
	for (i=0; i<pool->size; i++) 
		samwrite(stream, elem[i]);
}

