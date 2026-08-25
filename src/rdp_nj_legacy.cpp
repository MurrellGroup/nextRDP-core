#include "MathFuncsDll.h"

#include <cmath>
#include <cstdlib>
#include <cstdint>

// This is the literal Clearcut/NJ execution path used by RDP's MakeNJTreesP2.
// It is kept in its own translation unit so its floating-point evaluation can
// be matched without changing the preceding distance/matrix routines.
namespace MathFuncs::rdp_nj_legacy {

/* The original 32-bit RDP build evaluates the Clearcut expressions on an
 * x87 stack.  WASM has no x87 register format, so these helpers retain the
 * 64-bit x87 significand after each arithmetic operation until the source
 * explicitly stores a float. */
static inline long double x87_round(const long double input) {
    if (!std::isfinite(input) || input == 0.0L) return input;
    int exponent = 0;
    const long double mantissa = std::frexp(input, &exponent);
    const long double scaled = std::ldexp(mantissa, 64);
    const long double magnitude = std::fabs(scaled);
    const long double lower = std::floor(magnitude);
    const long double fraction = magnitude - lower;
    long double rounded = lower;
    if (fraction > 0.5L ||
        (fraction == 0.5L && std::fmod(lower, 2.0L) != 0.0L)) {
        rounded += 1.0L;
    }
    if (scaled < 0.0L) rounded = -rounded;
    return std::ldexp(rounded, exponent - 64);
}
static inline long double x87_add(const long double a, const long double b) {
    return x87_round(a + b);
}
static inline long double x87_sub(const long double a, const long double b) {
    return x87_round(a - b);
}
static inline long double x87_div(const long double a, const long double b) {
    return x87_round(a / b);
}

static inline float f32_add(const float a, const float b) {
    volatile float result = a + b;
    return result;
}

static inline float x87_add_store(const float a, const long double b) {
    return static_cast<float>(x87_add(static_cast<long double>(a), b));
}

DMAT *
NJ_parse_distance_matrix(float *dists,int nextno) {

	DMAT *dmat = NULL;

	int state, dmat_type;
	int row;
	int fltcnt;
	int x, y, i, xx, j;
	int numvalread;
	int expectedvalues = -1;
	float val;
	int first_state = 0;


	/* allocate our distance matrix and token structure */
	dmat = (DMAT *)calloc(1, sizeof(DMAT));
	
	dmat->ntaxa = nextno + 1;

	/* set our initial working size according to the # of taxa */
	dmat->size = dmat->ntaxa;

	/* allocate space for the distance matrix values here */
	dmat->val =
		(float *)calloc(NJ_NCELLS(dmat->ntaxa), sizeof(float));
	

	/*  taxa names */
	dmat->taxaname = (char **)calloc(dmat->ntaxa, sizeof(char *));
	
	for (xx = 0; xx <= nextno; xx++) {
		dmat->taxaname[xx] = (char *)calloc(3, sizeof(char));
	}

	/* set the initial state of our state machine */
	dmat_type = NJ_PARSE_SYMMETRIC;
	row = -1;
	fltcnt = 0;
	numvalread = 0;


	//this is where clearcut parses the distance matrix file


	//dmat-> val = dists;

	for (i= 0; i < nextno; i++){
		for (j = i+1; j<=nextno;j++)
			dmat->val[NJ_MAP(i, j, dmat->size)] = dists[i+j*(nextno+1)];
	}




	/* now lets allocate space for the r and r2 columns */
	dmat->r = (float *)calloc(dmat->ntaxa, sizeof(float));
	dmat->r2 = (float *)calloc(dmat->ntaxa, sizeof(float));
	
	/* track some memory addresses */
	dmat->rhandle = dmat->r;
	dmat->r2handle = dmat->r2;
	dmat->valhandle = dmat->val;

	return(dmat);



	/* clean up our partial progress */

}


void
NJ_init_r(DMAT *dmat) {

	std::int32_t i, j, size;
	std::int32_t index;
	float *r, *r2, *val;
	std::int32_t size1;
	long double size2;

	r = dmat->r;
	r2 = dmat->r2;
	val = dmat->val;
	size = dmat->size;
	size1 = size - 1;
	size2 = size - 2;

	index = 0;
	for (i = 0; i<size1; i++) {
		index++;
		for (j = i + 1; j<size; j++) {
			r[i] = f32_add(r[i], val[index]);
			r[j] = f32_add(r[j], val[index]);
			index++;
		}

		r2[i] = static_cast<float>(x87_div(r[i], size2));
	}

	return;
}

NJ_VERTEX *
NJ_init_vertex(DMAT *dmat) {

	std::int32_t i;
	NJ_VERTEX *vertex;

	/* allocate the vertex here */
	vertex = (NJ_VERTEX *)calloc(1, sizeof(NJ_VERTEX));

	/* allocate the nodes in the vertex */
	vertex->nodes = (NJ_TREE **)calloc(dmat->ntaxa, sizeof(NJ_TREE *));
	vertex->nodes_handle = vertex->nodes;

	/* initialize our size and active variables */
	vertex->nactive = dmat->ntaxa;
	vertex->size = dmat->ntaxa;

	/* initialize the nodes themselves */
	for (i = 0; i<dmat->ntaxa; i++) {

		vertex->nodes[i] = (NJ_TREE *)calloc(1, sizeof(NJ_TREE));

		vertex->nodes[i]->left = NULL;
		vertex->nodes[i]->right = NULL;

		vertex->nodes[i]->taxa_index = i;
	}

	return(vertex);
}

float
NJ_min_transform(DMAT *dmat,
	std::int32_t *ret_i,
	std::int32_t *ret_j) {

	std::int32_t i, j;   /* indices used for looping        */
	std::int32_t tmp_i = 0;/* to limit pointer dereferencing  */
	std::int32_t tmp_j = 0;/* to limit pointer dereferencing  */
	long double smallest;  /* track the smallest trans. dist  */
	long double curval;    /* the current trans. dist in loop */

	float *ptr;      /* pointer into distance matrix    */
	float *r2;       /* pointer to r2 matrix for computing transformed dists */

	smallest = (float)HUGE_VAL;

	/* track these here to limit pointer dereferencing in inner loop */
	ptr = dmat->val;
	r2 = dmat->r2;

	/* for every row */
	for (i = 0; i<dmat->size; i++) {
		ptr++;  /* skip diagonal */
		for (j = i + 1; j<dmat->size; j++) {   /* for every column */

											   /* find transformed distance in matrix at i, j */
			curval = x87_sub(static_cast<long double>(*(ptr++)),
				x87_add(static_cast<long double>(r2[i]),
					static_cast<long double>(r2[j])));

			/* if the transformed distanance is less than the known minimum */
			if (curval < smallest) {

				smallest = x87_round(curval);
				tmp_i = i;
				tmp_j = j;
			}
		}
	}

	/* pass back (by reference) the coords of the min. transformed distance */
	*ret_i = tmp_i;
	*ret_j = tmp_j;

	return(static_cast<float>(smallest));  /* return the min transformed distance */
}


NJ_TREE *
NJ_decompose(DMAT *dmat,
	NJ_VERTEX *vertex,
	std::int32_t x,
	std::int32_t y,
	int last_flag) {

	NJ_TREE *new_node;
	long double x2clade, y2clade;

	/* compute the distance from the clade components to the new node */
	if (last_flag) {
		x2clade =
			(dmat->val[NJ_MAP(x, y, dmat->size)]);
	}
	else {
		x2clade = x87_add(
			x87_div(dmat->val[NJ_MAP(x, y, dmat->size)], 2.0L),
			x87_div(x87_sub(dmat->r2[x], dmat->r2[y]), 2.0L));
	}

	vertex->nodes[x]->dist = x2clade;

	if (last_flag) {
		y2clade =
			(dmat->val[NJ_MAP(x, y, dmat->size)]);
	}
	else {
		y2clade = x87_add(
			x87_div(dmat->val[NJ_MAP(x, y, dmat->size)], 2.0L),
			x87_div(x87_sub(dmat->r2[y], dmat->r2[x]), 2.0L));
	}

	vertex->nodes[y]->dist = y2clade;

	/* allocate new node to connect two sub-clades */
	new_node = (NJ_TREE *)calloc(1, sizeof(NJ_TREE));

	new_node->left = vertex->nodes[x];
	new_node->right = vertex->nodes[y];
	new_node->taxa_index = NJ_INTERNAL_NODE;  /* this is not a terminal node, no taxa index */

	if (last_flag) {
		return(new_node);
	}

	vertex->nodes[x] = new_node;
	vertex->nodes[y] = vertex->nodes[0];

	vertex->nodes = &(vertex->nodes[1]);

	vertex->nactive--;

	return(new_node);
}

static inline
void
NJ_compute_r(DMAT *dmat,
	std::int32_t a,
	std::int32_t b) {

	std::int32_t i;         /* a variable used in indexing */
	float *ptrx, *ptry; /* pointers into the distance matrix */

						/* some variables to limit pointer dereferencing in loop */
	std::int32_t size;
	float *r, *val;

	/* to limit pointer dereferencing */
	size = dmat->size;
	val = dmat->val;
	r = dmat->r + a + 1;

	/*
	* Loop through the rows and decrement the stored r values
	* by the distances stored in the rows and columns of the distance
	* matrix which are being removed post-join.
	*
	* We do the rows altogether in order to benefit from cache locality.
	*/
	ptrx = &(val[NJ_MAP(a, a + 1, size)]);
	ptry = &(val[NJ_MAP(b, b + 1, size)]);

	for (i = a + 1; i<size; i++) {
		long double updated = x87_sub(static_cast<long double>(*r),
			static_cast<long double>(*ptrx++));
		if (i>b) {
			updated = x87_sub(updated, static_cast<long double>(*ptry++));
		}
		*r = static_cast<float>(updated);

		r++;
	}

	/* Similar to the above loop, we now do the columns */
	ptrx = &(val[NJ_MAP(0, a, size)]);
	ptry = &(val[NJ_MAP(0, b, size)]);
	r = dmat->r;
	for (i = 0; i<b; i++) {
		long double updated = static_cast<long double>(*r);
		if (i<a) {
			updated = x87_sub(updated, static_cast<long double>(*ptrx));
			ptrx += size - i - 1;
		}

		updated = x87_sub(updated, static_cast<long double>(*ptry));
		*r = static_cast<float>(updated);
		ptry += size - i - 1;
		r++;
	}

	return;
}


static inline
void
NJ_collapse(DMAT *dmat,
	NJ_VERTEX *vertex,
	std::int32_t a,
	std::int32_t b) {


	std::int32_t i;     /* index used for looping */
	std::int32_t size;  /* size of dmat --> reduce pointer dereferencing */
	long double a2clade;  /* distance from a to the new node that joins a and b */
	long double b2clade;  /* distance from b to the new node that joins a and b */
	long double cval;     /* stores distance information during loop */
	float *vptr;    /* pointer to elements in first row of dist matrix */
	float *ptra;    /* pointer to elements in row a of distance matrix */
	float *ptrb;    /* pointer to elements in row b of distance matrix */

	float *val, *r, *r2;  /* simply used to limit pointer dereferencing */


	

	/* some shortcuts to help limit dereferencing */
	val = dmat->val;
	r = dmat->r;
	r2 = dmat->r2;
	size = dmat->size;

	/* compute the distance from the clade components (a, b) to the new node */
	a2clade = x87_div(x87_add(val[NJ_MAP(a, b, size)],
		x87_sub(dmat->r2[a], dmat->r2[b])), 2.0L);
	b2clade = x87_div(x87_add(val[NJ_MAP(a, b, size)],
		x87_sub(dmat->r2[b], dmat->r2[a])), 2.0L);


	r[a] = 0.0;  /* we are removing row a, so clear dist. in r */

				 /*
				 * Fill the horizontal part of the "a" row and finish computing r and r2
				 * we handle the horizontal component first to maximize cache locality
				 */
	ptra = &(val[NJ_MAP(a, a + 1, size)]);   /* start ptra at the horiz. of a  */
	ptrb = &(val[NJ_MAP(a + 1, b, size)]);   /* start ptrb at comparable place */
	for (i = a + 1; i<size; i++) {

		/*
		* Compute distance from new internal node to others in
		* the distance matrix.
		*/
		cval = x87_div(x87_add(x87_sub(*ptra, a2clade),
			x87_sub(*ptrb, b2clade)), 2.0L);

		/* incr.  row b pointer differently depending on where i is in loop */
		if (i<b) {
			ptrb += size - i - 1;  /* traverse vertically  by incrementing by row */
		}
		else {
			ptrb++;            /* traverse horiz. by incrementing by column   */
		}

		/* assign the newly computed distance and increment a ptr by a column */
		*(ptra++) = cval;

		/* accumulate the distance onto the r vector */
		r[a] = x87_add_store(r[a], cval);
		const long double updated_i = x87_add(r[i], cval);
		r[i] = static_cast<float>(updated_i);

		/* scale r2 on the fly here */
		r2[i] = static_cast<float>(x87_div(updated_i, size - 3));
	}

	/* fill the vertical part of the "a" column and finish computing r and r2 */
	ptra = val + a;  /* start at the top of the columb for "a" */
	ptrb = val + b;  /* start at the top of the columb for "b" */
	for (i = 0; i<a; i++) {

		/*
		* Compute distance from new internal node to others in
		* the distance matrix.
		*/
		cval = x87_div(x87_add(x87_sub(*ptra, a2clade),
			x87_sub(*ptrb, b2clade)), 2.0L);

		/* assign the newly computed distance and increment a ptr by a column */
		*ptra = cval;

		/* accumulate the distance onto the r vector */
		r[a] = x87_add_store(r[a], cval);
		const long double updated_i = x87_add(r[i], cval);
		r[i] = static_cast<float>(updated_i);

		/* scale r2 on the fly here */
		r2[i] = static_cast<float>(x87_div(updated_i, size - 3));

		/* here, always increment by an entire row */
		ptra += size - i - 1;
		ptrb += size - i - 1;
	}


	/* scale r2 on the fly here */
	r2[a] = static_cast<float>(x87_div(r[a], size - 3));



	/*
	* Copy row 0 into row b.  Again, the code is structured into two
	* loops to maximize cache locality for writes along the horizontal
	* component of row b.
	*/
	vptr = val;
	ptrb = val + b;
	for (i = 0; i<b; i++) {
		*ptrb = *(vptr++);
		ptrb += size - i - 1;
	}
	vptr++;  /* skip over the diagonal */
	ptrb = &(val[NJ_MAP(b, b + 1, size)]);
	for (i = b + 1; i<size; i++) {
		*(ptrb++) = *(vptr++);
	}

	/*
	* Collapse r here by copying contents of r[0] into r[b] and
	* incrementing pointer to the beginning of r by one row
	*/
	r[b] = r[0];
	dmat->r = r + 1;


	/*
	* Collapse r2 here by copying contents of r2[0] into r2[b] and
	* incrementing pointer to the beginning of r2 by one row
	*/
	r2[b] = r2[0];
	dmat->r2 = r2 + 1;

	/* increment dmat pointer to next row */
	dmat->val += size;

	/* decrement the total size of the distance matrix by one row */
	dmat->size--;

	return;
}




NJ_TREE *NJ_neighbor_joining(DMAT *dmat, int outlyer) {


	NJ_TREE   *tree = NULL;
	NJ_VERTEX *vertex = NULL;

	std::int32_t a, b;
	float min;


	/* initialize the r and r2 vectors */
	NJ_init_r(dmat);

	/* allocate and initialize our vertex vector used for tree construction */
	vertex = NJ_init_vertex(dmat);
	

	/* we iterate until the working distance matrix has only 2 entries */
	while (vertex->nactive > 2) {

		/*
		* Find the global minimum transformed distance from the distance matrix
		*/
		min = NJ_min_transform(dmat, &a, &b);

		/*
		* Build the tree by removing nodes a and b from the vertex array
		* and inserting a new internal node which joins a and b.  Collapse
		* the vertex array similarly to how the distance matrix and r and r2
		* are compacted.
		*/
		NJ_decompose(dmat, vertex, a, b, 0);

		/* decrement the r and r2 vectors by the distances corresponding to a, b */
		NJ_compute_r(dmat, a, b);

		/* compact the distance matrix and the r and r2 vectors */
		NJ_collapse(dmat, vertex, a, b);
	}

	/* Properly join the last two nodes on the vertex list */
	tree = NJ_decompose(dmat, vertex, 0, 1, NJ_LAST);

	/* return the computed tree to the calling function */
	return(tree);
}


static inline
void
NJ_permute(std::int32_t *perm,
	std::int32_t size) {

	std::int32_t i;     /* index used for looping */
	std::int32_t swap;  /* we swap values to generate permutation */
	std::int32_t tmp;   /* used for swapping values */
	double K;

					/* check to see if vector of long ints is valid */
	if (!perm) {
		
		exit(-1);
	}

	/* init permutation as an ordered list of integers */
	for (i = 0; i<size; i++) {
		perm[i] = i;
	}

	/*
	* Iterate across the array from i = 0 to size -1, swapping ith element
	* with a randomly chosen element from a changing range of possible values
	*/
	for (i = 0; i<size; i++) {

		/* choose which element we will swap with */
		K = rand();
		swap = i + (int)((K / RAND_MAX)*(size - 1));
		//swap = i + NJ_genrand_int31_top(size - i);

		/* swap elements here */
		if (i != swap) {
			tmp = perm[swap];
			perm[swap] = perm[i];
			perm[i] = tmp;
		}
	}

	return;
}

static inline
float
NJ_find_hmin(DMAT *dmat,
	std::int32_t a,
	std::int32_t *min,
	std::int32_t *hmincount) {

	std::int32_t i;     /* index variable for looping                    */
	int size;       /* current size of distance matrix               */
	int mindex = 0; /* holds the current index to the chosen minimum */
	float curval;   /* used to hold current transformed values       */
	float hmin;     /* the value of the transformed minimum          */

	float *ptr, *r2, *val;  /* pointers used to reduce dereferencing in inner loop */

							/* values used for stochastic selection among multiple minima */
	float p, x;
	std::int32_t smallcnt;

	/* initialize the min to something large */
	hmin = (float)HUGE_VAL;

	/* setup some pointers to limit dereferencing later */
	r2 = dmat->r2;
	val = dmat->val;
	size = dmat->size;

	/* initialize values associated with minima tie breaking */
	p = 1.0;
	smallcnt = 0;


	ptr = &(val[NJ_MAP(a, a + 1, size)]);   /* at the start of the horiz. part */
	for (i = a + 1; i<size; i++) {

		curval = *(ptr++) - (r2[a] + r2[i]);  /* compute transformed distance */

		if (NJ_FLT_EQ(curval, hmin)) {  /* approx. equal */

			smallcnt++;

			p = 1.0 / (float)smallcnt;
			x = rand()/RAND_MAX;
			//x = genrand_real2();

			/* select this minimum in a way which is proportional to
			the number of minima found along the row so far */
			if (x < p) {
				mindex = i;
			}

		}
		else if (curval < hmin) {

			smallcnt = 1;
			hmin = curval;
			mindex = i;
		}
	}

	/* save off the the minimum index to be returned via reference */
	*min = mindex;

	/* save off the number of minima */
	*hmincount = smallcnt;

	/* return the value of the smallest tranformed distance */
	return(hmin);
}


static inline
float
NJ_find_vmin(DMAT *dmat,
	std::int32_t a,
	std::int32_t *min,
	std::int32_t *vmincount) {

	std::int32_t i;         /* index variable used for looping */
	std::int32_t size;      /* track the size of the matrix    */
	std::int32_t mindex = 0;/* track the index to the minimum  */
	float curval;       /* track value of current transformed distance  */
	float vmin;         /* the index to the smallest "vertical" minimum */

						/* pointers which are used to reduce pointer dereferencing in inner loop */
	float *ptr, *r2, *val;

	/* values used in stochastically breaking ties */
	float p, x;
	std::int32_t smallcnt;

	/* initialize the vertical min to something really big */
	vmin = (float)HUGE_VAL;

	/* save off some values to limit dereferencing later */
	r2 = dmat->r2;
	val = dmat->val;
	size = dmat->size;

	p = 1.0;
	smallcnt = 0;

	/* start on the first row and work down */
	ptr = &(val[NJ_MAP(0, a, size)]);
	for (i = 0; i<a; i++) {

		curval = *ptr - (r2[i] + r2[a]);  /* compute transformed distance */

		if (NJ_FLT_EQ(curval, vmin)) {  /* approx. equal */

			smallcnt++;

			p = 1.0 / (float)smallcnt;
			x = rand()/RAND_MAX;
			//x = genrand_real2();

			/* break ties stochastically to avoid systematic bias */
			if (x < p) {
				mindex = i;
			}

		}
		else if (curval < vmin) {

			smallcnt = 1;
			vmin = curval;
			mindex = i;
		}

		/* increment our working pointer to the next row down */
		ptr += size - i - 1;
	}

	/* pass back the index to the minimum found so far (by reference) */
	*min = mindex;

	/* pass back the number of minima along the vertical */
	*vmincount = smallcnt;

	/* return the value of the smallest transformed distance */
	return(vmin);
}

static inline
int
NJ_check_additivity(DMAT *dmat,
	std::int32_t a,
	std::int32_t b) {

	float a2clade, b2clade;
	float clade_dist;
	std::int32_t target;


	/* determine target taxon here */
	if (b == dmat->size - 1) {
		/* if we can't do a row here, lets do a column */
		if (a == 0) {
			if (b == 1) {
				target = 2;
			}
			else {
				target = 1;
			}
		}
		else {
			target = 0;
		}
	}
	else {
		target = b + 1;
	}


	/* distance between a and the root of clade (a,b) */
	a2clade =
		((dmat->val[NJ_MAP(a, b, dmat->size)]) +
		(dmat->r2[a] - dmat->r2[b])) / 2.0;

	/* distance between b and the root of clade (a,b) */
	b2clade =
		((dmat->val[NJ_MAP(a, b, dmat->size)]) +
		(dmat->r2[b] - dmat->r2[a])) / 2.0;

	/* distance between the clade (a,b) and the target taxon */
	if (b<target) {

		/* compute the distance from the clade root to the target */
		clade_dist =
			((dmat->val[NJ_MAP(a, target, dmat->size)] - a2clade) +
			(dmat->val[NJ_MAP(b, target, dmat->size)] - b2clade)) / 2.0;

		/*
		* Check to see that distance from clade root to target + distance from
		*  b to clade root are equal to the distance from b to the target
		*/
		if (NJ_FLT_EQ(dmat->val[NJ_MAP(b, target, dmat->size)],
			(clade_dist + b2clade))) {
			return(1);  /* join is legitimate   */
		}
		else {
			return(0);  /* join is illigitimate */
		}

	}
	else {

		/* compute the distance from the clade root to the target */
		clade_dist =
			((dmat->val[NJ_MAP(target, a, dmat->size)] - a2clade) +
			(dmat->val[NJ_MAP(target, b, dmat->size)] - b2clade)) / 2.0;

		/*
		* Check to see that distance from clade root to target + distance from
		*  b to clade root are equal to the distance from b to the target
		*/
		if (NJ_FLT_EQ(dmat->val[NJ_MAP(target, b, dmat->size)],
			(clade_dist + b2clade))) {
			return(1);  /* join is legitimate   */
		}
		else {
			return(0);  /* join is illegitimate */
		}
	}
}

static inline
int
NJ_check(int RJ, DMAT *dmat,
	std::int32_t a,
	std::int32_t b,
	float min,
	int additivity) {


	std::int32_t i, size;
	float *ptr, *val, *r2;


	/* some aliases for speed and readability reasons */
	val = dmat->val;
	r2 = dmat->r2;
	size = dmat->size;


	/* now determine if joining a, b will result in broken distances */
	if (additivity) {
		if (!NJ_check_additivity(dmat, a, b)) {
			return(0);
		}
	}

	/* scan the horizontal of row b, punt if anything < min */
	ptr = &(val[NJ_MAP(b, b + 1, size)]);
	for (i = b + 1; i<size; i++) {
		if (NJ_FLT_LT((*ptr - (r2[b] + r2[i])), min)) {
			return(0);
		}
		ptr++;
	}

	/* scan the vertical component of row a, punt if anything < min */
	if (RJ == 1) {  /* if we are doing random joins, we checked this */
		ptr = val + a;
		for (i = 0; i<a; i++) {
			if (NJ_FLT_LT((*ptr - (r2[i] + r2[a])), min)) {
				return(0);
			}
			ptr += size - i - 1;
		}
	}

	/* scan the vertical component of row b, punt if anything < min */
	ptr = val + b;
	for (i = 0; i<b; i++) {
		if (NJ_FLT_LT((*ptr - (r2[i] + r2[b])), min) && i != a) {
			return(0);
		}
		ptr += size - i - 1;
	}

	return(1);
}


void
NJ_free_vertex(NJ_VERTEX *vertex) {

	if (vertex) {
		if (vertex->nodes_handle) {
			free(vertex->nodes_handle);
		}
		free(vertex);
	}

	return;
}

NJ_TREE *
NJ_relaxed_nj(int RJ, DMAT *dmat) {


	NJ_TREE *tree;
	NJ_VERTEX *vertex;
	std::int32_t a, b, t, bh, bv, i;
	float hmin, vmin, hvmin;
	float p, q, x;
	int join_flag;
	int additivity_mode;
	std::int32_t hmincount, vmincount;
	std::int32_t *permutation = NULL;



	/* initialize the r and r2 vectors */
	NJ_init_r(dmat);

	additivity_mode = 1;

	
	/* allocate and initialize our vertex vector used for tree construction */
	vertex = NJ_init_vertex(dmat);

	/* loop until there are only 2 nodes left to join */
	while (vertex->nactive > 2) {

		switch (RJ) {

			/* RANDOMIZED JOINS */
		case 0:

			join_flag = 0;

			NJ_permute(permutation, dmat->size - 1);
			for (i = 0; i<dmat->size - 1 && (vertex->nactive>2); i++) {

				a = permutation[i];

				/* find min trans dist along horiz. of row a */
				hmin = NJ_find_hmin(dmat, a, &bh, &hmincount);
				if (a) {
					/* find min trans dist along vert. of row a */
					vmin = NJ_find_vmin(dmat, a, &bv, &vmincount);
				}
				else {
					vmin = hmin;
					bv = bh;
					vmincount = 0;
				}

				if (NJ_FLT_EQ(hmin, vmin)) {

					/*
					* The minima along the vertical and horizontal are
					* the same.  Compute the proportion of minima along
					* the horizonal (p) and the proportion of minima
					* along the vertical (q).
					*
					* If the same minima exist along the horizonal and
					* vertical, we break the tie in a way which is
					* non-biased.  That is, we break the tie based on the
					* proportion of horiz. minima versus vertical minima.
					*
					*/
					p = (float)hmincount / ((float)hmincount + (float)vmincount);
					q = 1.0 - p;
					x = rand()/RAND_MAX;// genrand_real2();

					if (x < p) {
						hvmin = hmin;
						b = bh;
					}
					else {
						hvmin = vmin;
						b = bv;
					}
				}
				else if (NJ_FLT_LT(hmin, vmin)) {
					hvmin = hmin;
					b = bh;
				}
				else {
					hvmin = vmin;
					b = bv;
				}

				if (NJ_check(RJ, dmat, a, b, hvmin, additivity_mode)) {

					/* swap a and b, if necessary, to make sure a < b */
					if (b < a) {
						t = a;
						a = b;
						b = t;
					}

					join_flag = 1;

					/* join taxa from rows a and b */
					NJ_decompose(dmat, vertex, a, b, 0);

					/* collapse matrix */
					NJ_compute_r(dmat, a, b);
					NJ_collapse(dmat, vertex, a, b);

					NJ_permute(permutation, dmat->size - 1);
				}
			}

			/* turn off additivity if go through an entire cycle without joining */
			if (!join_flag) {
				additivity_mode = 0;
			}

			break;



			/* DETERMINISTIC JOINS */
		case 1:

			join_flag = 0;

			for (a = 0; a<dmat->size - 1 && (vertex->nactive > 2);) {

				/* find the min along the horizontal of row a */
				hmin = NJ_find_hmin(dmat, a, &b, &hmincount);

				if (NJ_check(RJ, dmat, a, b, hmin, additivity_mode)) {

					join_flag = 1;

					/* join taxa from rows a and b */
					NJ_decompose(dmat, vertex, a, b, 0);

					/* collapse matrix */
					NJ_compute_r(dmat, a, b);
					NJ_collapse(dmat, vertex, a, b);

					if (a) {
						a--;
					}

				}
				else {
					a++;
				}
			}

			/* turn off additivity if go through an entire cycle without joining */
			if (!join_flag) {
				additivity_mode = 0;
			}

			break;
		}

	}  /* WHILE */

	   /* Join the last two nodes on the vertex list */
	tree = NJ_decompose(dmat, vertex, 0, 1, NJ_LAST);

	if (vertex) {
		NJ_free_vertex(vertex);
	}

	

	return(tree);
}

void
NJ_output_tree2(int nlen, int *tpos, char *otreex,
	NJ_TREE *tree,
	NJ_TREE *root,
	DMAT *dmat) {


	int   x2, x3, n2, n3, mod;
	float x;

	if (!tree) {
		return;
	}

	if (tree->taxa_index != NJ_INTERNAL_NODE) {

		//print the taxon name (actually its number starting with an S)
		*tpos = *tpos + 1;
		*(otreex + *tpos) = 83;

		mod = nlen;
		n2 = tree->taxa_index;
		n3 = (int)(n2 / mod);
		*tpos = *tpos + 1;
		*(otreex + *tpos) = (char)(48 + n3);
		while (mod > 1) {

			n2 -= n3*mod;
			n3 = (int)(n2 / (mod / 10));
			*tpos = *tpos + 1;
			*(otreex + *tpos) = (char)(48 + n3);
			mod /= 10;
		}
		//Print the branchlength
		*tpos = *tpos + 1;
		*(otreex + *tpos) = 58;//first print the colon
		
		x = tree->dist;
		if (x < 0.0) {
			x *= -1;
			//*tpos = *tpos + 1;
			//*(otreex + *tpos) = 45;//a minus sign
		}
		if (x < 1.0) {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = 48;// a zero
			*tpos = *tpos + 1;
			*(otreex + *tpos) = 46;//a dot
		}
		else {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = (int)(x)+48; //a number>0
			*tpos = *tpos + 1;
			*(otreex + *tpos) = 46;//a dot
		}

		x -= (int)(x);

		mod = 100000;
		x2 = (int)(x*mod);
		x3 = (int)(x2 / (mod / 10));
		while (mod > 10) {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = char(48 + x3);
			mod /= 10;
			x2 -= x3*mod;
			x3 = x2 / (mod / 10);
			//if (x2/10000000 < 1 && (double)(x2/10000000) > 0.1){
			//	tpos++;
			//	*(otree + tpos) = 46;
			//}
		}

//			fprintf(fp, "%s:%f",
//				dmat->taxaname[tree->taxa_index],
//				tree->dist);//print length of terminal branch in normal notation
		

	}
	else {


		if (tree->left && tree->right) {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = char(40);//open bracket
			//fprintf(fp, "(");
		}
		if (tree->left) {
			NJ_output_tree2(nlen, tpos, otreex, tree->left, root, dmat);
		}

		if (tree->left && tree->right) {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = char(44);
			//fprintf(fp, ",");
		}
		if (tree->right) {
			NJ_output_tree2(nlen, tpos, otreex, tree->right, root, dmat);
		}

		if (tree != root->left) {
			if (tree->left && tree->right) {
				if (tree != root) {
					//print close bracket and colon
					*tpos = *tpos + 1;
					*(otreex + *tpos) = char(41);
					*tpos = *tpos + 1;
					*(otreex + *tpos) = char(58);
					//now print the branch length
					x = tree->dist;
					if (x < 0.0) {
						x *= -1;
						//*tpos = *tpos + 1;
						//*(otreex + *tpos) = 45;//a minus sign
					}
					if (x < 1.0) {
						*tpos = *tpos + 1;
						*(otreex + *tpos) = 48;// a zero
						*tpos = *tpos + 1;
						*(otreex + *tpos) = 46;//a dot
					}
					else {
						*tpos = *tpos + 1;
						*(otreex + *tpos) = (int)(x)+48; //a number>0
						*tpos = *tpos + 1;
						*(otreex + *tpos) = 46;//a dot
					}

					x -= (int)(x);

					mod = 100000;
					x2 = (int)(x*mod);
					x3 = (int)(x2 / (mod / 10));
					while (mod > 10) {
						*tpos = *tpos + 1;
						*(otreex + *tpos) = char(48 + x3);
						mod /= 10;
						x2 -= x3*mod;
						x3 = x2 / (mod / 10);
						//if (x2/10000000 < 1 && (double)(x2/10000000) > 0.1){
						//	tpos++;
						//	*(otree + tpos) = 46;
						//}
					}

					//fprintf(fp, "):%f", tree->dist);//close bracket and print branch length in normal notation
					
				}
				else {
					*tpos = *tpos + 1;
					*(otreex + *tpos) = char(41);
					//fprintf(fp, ")");
				}
			}
		}
		else {
			*tpos = *tpos + 1;
			*(otreex + *tpos) = char(41);
			//fprintf(fp, ")");
		}
	}

	return;
}

void
NJ_search_tree(int nlen, int *tpos, char *otreex,
	NJ_TREE *tree,
	NJ_TREE *root,
	DMAT *dmat, NJ_TREE *target, int outlyer) {


	int   x2, x3, n2, n3, mod;
	float x;

	if (!tree) {
		return;
	}

	if (tree->taxa_index == outlyer) {

		target = tree;


	}
	else {


		
		if (tree->left) {
			NJ_search_tree(nlen, tpos, otreex, tree->left, root, dmat,target,outlyer);
		}

		
		if (tree->right) {
			NJ_search_tree(nlen, tpos, otreex, tree->right, root, dmat, target,outlyer);
		}

		
	}

	return;
}

int
NJ_output_tree(int nlen, int *tpos, NJ_TREE *tree,
	DMAT *dmat,
	std::int32_t count, char *outtree,int outlyer) {
	//find the outlyer node
	NJ_TREE *target;

	target = (NJ_TREE *)calloc(1, sizeof(NJ_TREE));
	if (outlyer > 0) {
		NJ_search_tree(nlen, tpos, outtree, tree, tree, dmat, target, outlyer-1);
		NJ_output_tree2(nlen, tpos, outtree, target, tree, dmat);
	}
	else
		NJ_output_tree2(nlen, tpos, outtree, tree, tree, dmat);
	
	*tpos = *tpos + 1;
	*(outtree + *tpos) = char(59);
	// add the final ";"
	//fprintf(fp, ";\n");

	
	return(0);
}

void
NJ_free_tree(NJ_TREE *node) {

	if (!node) {
		return;
	}

	if (node->left) {
		NJ_free_tree(node->left);
	}

	if (node->right) {
		NJ_free_tree(node->right);
	}

	free(node);

	return;
}

void
NJ_free_dmat(DMAT *dmat) {

	std::int32_t i;

	if (dmat) {

		if (dmat->taxaname) {

			for (i = 0; i<dmat->ntaxa; i++) {
				if (dmat->taxaname[i]) {
					free(dmat->taxaname[i]);
				}
			}

			free(dmat->taxaname);
		}

		if (dmat->valhandle) {
			free(dmat->valhandle);
		}

		if (dmat->rhandle) {
			free(dmat->rhandle);
		}

		if (dmat->r2handle) {
			free(dmat->r2handle);
		}

		free(dmat);
	}

	return;
}



//All this clearcut code is written by Luke Sheneman
/*
* clearcut.c
*
* $Id: clearcut.c,v 1.2 2006/08/25 03:58:45 sheneman Exp $
*
*****************************************************************************
*
* Copyright (c) 2004,  Luke Sheneman
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
*  + Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*  + Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in
*    the documentation and/or other materials provided with the
*    distribution.
*  + The names of its contributors may not be used to endorse or promote
*    products derived  from this software without specific prior
*    written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************
*
* An implementation of the Relaxed Neighbor-Joining algorithm
*  of Evans, J., Sheneman, L., and Foster, J.
*
*
* AUTHOR:
*
*   Luke Sheneman
*   sheneman@cs.uidaho.edu
*
*/
float Clearcut(int outlyer, int NextNo, int treetype, int nlen2, int nseed, int RJ, int UBD, float *dists, char *outtree) {
	//treetype = 1 for normal NJ, 0 for rapidNJ
	//RJ = 0 for randomized joins and 1 for deterministic joins during relaxed NJ construction
	int lt, nlen;
	int tpos;
	DMAT *dmat;         /* The working distance matrix */
	//DMAT *dmat_backup = NULL;/* A backup distance matrix    */
	NJ_TREE *tree;      /* The phylogenetic tree       */
	
	std::int32_t i;
	tpos = 0;
	/* some variables for tracking time */
	//struct timeval tv;
	unsigned long long startUs, endUs;


	/* check and parse supplied command-line arguments */
	//nj_args = NJ_handle_args(argc, argv); -these now get passed through the function call

	/* Initialize Mersenne Twister PRNG */
	//init_genrand(seed); already initialised
	srand(nseed);


	//switch (nj_args->input_mode) {

		/* If the input type is a distance matrix */
	///case NJ_INPUT_MODE_DISTANCE: input type always distances

		/* parse the distance matrix */
		dmat = NJ_parse_distance_matrix(dists, NextNo);
		

		/* If the input type is a multiple sequence alignment */
	
	//}

	/*
	* If we are going to generate multiple trees from
	* the same distance matrix, we need to make a backup
	* of the original distance matrix.
	*/
//	if (nj_args->ntrees > 1) {
//		dmat_backup = NJ_dup_dmat(dmat);
//	}

	/* process n trees */
		i = 0;
		

		/* RECORD THE PRECISE TIME OF THE START OF THE NEIGHBOR-JOINING */
		//gettimeofday(&tv, NULL);
		startUs = 10;//((unsigned long long) tv.tv_sec * 1000000ULL)
					 //  + ((unsigned long long) tv.tv_usec);


					 /*
					 * Invoke either the Relaxed Neighbor-Joining (treetype=0)
					 * or the "traditional" Neighbor-Joining algorithm(treetype=1)
					 */
		if (treetype == 1) {
			tree = NJ_neighbor_joining(dmat, outlyer);
		}
		else {
			tree = NJ_relaxed_nj(RJ, dmat);
		}

		

		/* RECORD THE PRECISE TIME OF THE END OF THE NEIGHBOR-JOINING */
		// gettimeofday(&tv, NULL);
		endUs = 20;// ((unsigned long long) tv.tv_sec * 1000000ULL)
				   //  + ((unsigned long long) tv.tv_usec);

		if (NextNo < 100)
			nlen = 10;
		else if (NextNo < 1000)
			nlen = 100;
		else
			nlen = 1000;

		/* Output the neighbor joining tree here */
		NJ_output_tree(nlen, &tpos, tree, dmat, i, outtree,outlyer);

		NJ_free_tree(tree);  /* Free the tree */
		NJ_free_dmat(dmat);  /* Free the working distance matrix */

		lt = tpos;

	

	return(lt);
}

} // namespace MathFuncs::rdp_nj_legacy

namespace MathFuncs {
float rdp_nj_clearcut_compat(int outlyer, int NextNo, int treetype, int nlen2,
                             int nseed, int RJ, int UBD, float *dists,
                             char *outtree) {
    return rdp_nj_legacy::Clearcut(outlyer, NextNo, treetype, nlen2, nseed,
                                   RJ, UBD, dists, outtree);
}
} // namespace MathFuncs
