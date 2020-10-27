#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "util.h"
#include "extract.h"
#include "grad.h"
#include "pi.h"
#include <time.h>
#include <omp.h>


#define POS_RES     0x01   /* 1st bit */
#define NEG_RES     0x02   /* 2nd bit */
#define VISITED     0x04   /* 3rd bit */
#define ACTIVE      0x08   /* 4th bit */
#define BRANCH_CUT  0x10   /* 5th bit */
#define BORDER      0x20   /* 6th bit */
#define UNWRAPPED   0x40   /* 7th bit */
#define POSTPONED   0x80   /* 8th bit */
#define RESIDUE     (POS_RES | NEG_RES)
#define AVOID       (BRANCH_CUT | BORDER)

int NUM_CORES;


double timediff(clock_t t1, clock_t t2) {
    double elapsed;
    elapsed = ((double)t2 - t1) / CLOCKS_PER_SEC * 1000;
    return elapsed;
}

/*
    Exchange two values
 */
void swap(int *a, int *b)
{
    int tmp;

    tmp = *a;
    *a = *b;
    *b = tmp;
}



/*
 *   Computute the x and y gradient.
 */
void Gradxy(float *phase, float *gradx, float *grady, int xsize, int ysize)
{
	int i, j, w;

	#pragma omp parallel for default(none) \
	private(j, i, w) \
	shared(xsize, ysize, phase, gradx, grady)
	for (j=0; j<ysize-1; j++)
	{
		for (i=0; i<xsize-1; i++)
		{
			w = j*xsize + i;
			gradx[w] = Gradient(phase[w], phase[w+1]);
			grady[w] = Gradient(phase[w], phase[w+xsize]);
		}
	}
}




/*
     Returns 0 if no pixels left, 1 otherwise
 */
int GetNextOneToUnwrap(int *a,
                       int *b,
                       int *index_list,
                       int *num_index,
                       int xsize,
                       int ysize)
{
    int index;
    if (*num_index < 1)
        return 0;   /* return if list empty */

    index = index_list[*num_index - 1];
    *a = index%xsize;
    *b = index/xsize;
    --(*num_index);
    return 1;
}






/*
 Insert new pixel into the list.
 Note: qual_map can be NULL
 */
void InsertList(float *soln,
                float val,
                unsigned char *bitflags,
                int a,
                int b,
                int *index_list,
                int *num_index,
                int xsize)
{
    int index;

    index = b*xsize + a;

    soln[index] = val;

    /* add to list */
    index_list[*num_index] = index;

    ++(*num_index);

    bitflags[index] |= UNWRAPPED;

    return;
}




/*
 Insert the four neighboring pixels of the given pixel
 (x,y) into the list.  The quality value of the given
 pixel is "val".
 */
void UpdateList(int x,
                int y,
                float val,
                float *phase,
                float *soln,
                unsigned char *bitflags,
                int xsize,
                int ysize,
                int *index_list,
                int *num_index)
{
    int    i, a, b, k, w;
    float  grad;

    a = x - 1;
    b = y;
    k = b*xsize + a;
    if (a >= 0
        && !(bitflags[k] & (BRANCH_CUT | UNWRAPPED | BORDER))) {
        w = y*xsize + x-1;
        grad = Gradient(phase[w], phase[w+1]);

        InsertList(soln, val + grad, bitflags, a, b,
                   index_list, num_index, xsize);
    }

    a = x + 1;
    b = y;
    k = b*xsize + a;
    if (a < xsize
        && !(bitflags[k] & (BRANCH_CUT | UNWRAPPED | BORDER))) {
        w = y*xsize + x;
        grad = - Gradient(phase[w], phase[w+1]);

        InsertList(soln, val + grad, bitflags, a, b,
                   index_list, num_index, xsize);
    }

    a = x;
    b = y - 1;
    k = b*xsize + a;
    if (b >= 0
        && !(bitflags[k] & (BRANCH_CUT | UNWRAPPED | BORDER))) {
        w = (y-1)*xsize + x;
        grad = Gradient(phase[w], phase[w+xsize]);

        InsertList(soln, val + grad, bitflags, a, b,
                   index_list, num_index, xsize);
    }

    a = x;
    b = y + 1;
    k = b*xsize + a;
    if (b < ysize
        && !(bitflags[k] & (BRANCH_CUT | UNWRAPPED | BORDER))) {
        w = y*xsize + x;
        grad = - Gradient(phase[w], phase[w+xsize]);

        InsertList(soln, val + grad, bitflags, a, b,
                   index_list, num_index, xsize);
    }
}



/* Unwrap the phase data (by Itoh's method) without crossing
 * any branch cuts.  Return number of disconnected pieces.
 */
int UnwrapAroundCutsGoldstein(float *phase,
                     unsigned char *bitflags,
                     float *soln,
                     int xsize,
                     int ysize,
                     int *path_order)
{
    int    i, j, k, a, b, c, n, num_pieces=0;
    float  value;
    int    num_index, max_list_size;
    int    *index_list;
    char   filename[300];

    max_list_size = xsize*ysize;
    AllocateInt(&index_list, max_list_size + 1, "bookkeeping list (index)");


    /* find starting point */
    n = 0;
    num_index = 0;

    for (j=0; j<ysize; j++)
    {
        for (i=0; i<xsize; i++)
        {
            k = j*xsize + i;
            if (!(bitflags[k] & (BRANCH_CUT | UNWRAPPED | BORDER)))
            {
                bitflags[k] |= UNWRAPPED;
                if (bitflags[k] & POSTPONED) /* soln[k] already stores the unwrapped value */
                    value = soln[k];
                else
                {
                    ++num_pieces;
                    value = soln[k] = phase[k];
                }

                UpdateList(i, j, value, phase, soln, bitflags, xsize,
                           ysize, index_list, &num_index);

                while (num_index > 0)
                {
                    ++n;

                    if (!GetNextOneToUnwrap(&a, &b, index_list,
                                            &num_index, xsize, ysize))
                        break;   /* no more to unwrap */

                    c = b*xsize + a;

                    // * * * * * * * * * * *
                    //   Save path order
                    //
                    path_order[c] = n;

                    bitflags[c] |= UNWRAPPED;
                    value = soln[c];
                    UpdateList(a, b, value, phase, soln, bitflags,
                               xsize, ysize, index_list, &num_index);
                }
            }
        }
    }

    free(index_list);

    /* unwrap branch cut pixels */
    for (j=1; j<ysize; j++)
    {
        for (i=1; i<xsize; i++)
        {
            k = j*xsize + i;

            if (bitflags[k] & AVOID)
            {
                if (!(bitflags[k-1] & AVOID))
                {
                    soln[k] = soln[k-1] + Gradient(phase[k], phase[k-1]);
                    path_order[k] = ++n;
                }
                else if (!(bitflags[k-xsize] & AVOID))
                {
                    soln[k] = soln[k-xsize] + Gradient(phase[k], phase[k-xsize]);
                    path_order[k] = ++n;
                }

            }
        }
    }



    return num_pieces;
}








/* Unwrap the phase data (by Itoh's method) without crossing
 * any branch cuts.  Return number of disconnected pieces.
 */
int UnwrapAroundCutsFrontier(float *phase,
                             unsigned char *bitflags,
                             float *soln,
                             int xsize,
                             int ysize,
                             int *path_order,
                             float *grady,
                             float *gradx,
                             int *list,
                             int length)
{
    int    i, j, k, kk, x, y, l, index, n=0, num_pieces=0;
    int    flag, base_in, base_out, top_in, top_out;
    float  value;


    /* find starting point */

    for (k=0; k<length; k++)
    {
        if (!(*(bitflags + k) & (BRANCH_CUT | UNWRAPPED | BORDER)))
        {
        	++num_pieces;

        	/* starting pixel */
            soln[k] = phase[k];
            flag = 1;

            /* base and top indexes */
            base_in = 0;
            base_out = xsize + ysize;
            top_in = base_in;
            top_out = base_out;

            *(list + top_in++) = k;

            /* Solve each level */
            while (flag)
            {
            	/* start integration level */
                for (l=base_in; l<top_in; l++)
                {
                    kk = *(list + l);
                    x = kk%xsize;
                    y = kk/xsize;

                    *(bitflags + kk) |= UNWRAPPED;
                    value = *(soln + kk);

                    /* save path order */
                    //*(path_order + kk) = n++;

                    /* neighbor pixels */

                    index = kk - 1;

                    if (x-1 >= 0
                        && !(*(bitflags + index) & (BRANCH_CUT | UNWRAPPED | BORDER)))
                    {
                        /* solution */
                        *(bitflags + index) |= UNWRAPPED;
                        *(soln + index) = value + *(gradx + index);
                        *(list + top_out++) = index;
                    }


                    index = kk + 1;

                    if (x+1 < xsize
                        && !(*(bitflags + index) & (BRANCH_CUT | UNWRAPPED | BORDER)))
                    {
                        /* solution */
                    	*(bitflags + index) |= UNWRAPPED;
                        *(soln + index) = value - *(gradx + kk);
                        *(list + top_out++) = index;
                    }


                    index = kk - xsize;

                    if (y-1 >= 0
                        && !(*(bitflags + index) & (BRANCH_CUT | UNWRAPPED | BORDER)))
                    {
                        /* solution */
                    	*(bitflags + index) |= UNWRAPPED;
                        *(soln + index) = value + *(grady + index);
                        *(list + top_out++) = index;
                    }


                    index = kk + xsize;

                    if (y+1 < ysize
                        && !(*(bitflags + index) & (BRANCH_CUT | UNWRAPPED | BORDER)))
                    {
                        /* solution */
                    	*(bitflags + index) |= UNWRAPPED;
                        *(soln + index) = value - *(grady + kk);
                        *(list + top_out++) = index;
                    }

                }
                /* end of level loop for */


                /* Exchange limits of the list */
                if (base_out==top_out)
                    flag = 0;
                else
                {
                    swap(&base_in, &base_out);
                    swap(&top_in, &top_out);
                    top_out = base_out;
                }
            }
            /* end of in loop while */
        }
        /* end of starting pixel search if */

    }
    /* end of starting for */



    /* unwrap branch cut pixels */

	#pragma omp parallel for default(none) \
	private(i, j, k) \
	shared(ysize, xsize, bitflags, soln, path_order, phase, n)
    for (j=1; j<ysize; j++)
    {
        for (i=1; i<xsize; i++)
        {
            k = j*xsize + i;

            if (bitflags[k] & AVOID)
            {
                if (!(bitflags[k-1] & AVOID))
                {
                    *(soln + k) = *(soln + k - 1) + Gradient(phase[k], phase[k-1]);
					//#pragma omp critical
                    //*(path_order + k) = n++;
                }
                else if (!(bitflags[k-xsize] & AVOID))
                {
                	*(soln + k) = *(soln + k - xsize) + Gradient(phase[k], phase[k-xsize]);
					//#pragma omp critical
                	//*(path_order + k) = n++;
                }
            }
        }
    }



    return num_pieces;
}








/* Place a branch cut in the bitflags array from pixel (a,b) */
/* to pixel (c,d).  The bit for the branch cut pixels is     */
/* given by the value of "code".                             */
void PlaceCut(unsigned char *array,
              int a,
              int b,
              int c,
              int d,
              int xsize,
              int ysize,
              int code)
{
    int  i, j, k, ii, jj, m, n, istep, jstep;
    double  r;

    /* residue location is upper-left corner of 4-square */
    if (c > a && a > 0) a++;
    else if (c < a && c > 0) c++;
    if (d > b && b > 0) b++;
    else if (d < b && d > 0) d++;

    if (a==c && b==d) {
        array[b*xsize + a] |= code;
        return;
    }
    m = (a < c) ? c - a : a - c;
    n = (b < d) ? d - b : b - d;
    if (m > n) {
        istep = (a < c) ? +1 : -1;
        r = ((double)(d - b))/((double)(c - a));
        for (i=a; i!=c+istep; i+=istep) {
            j = b + (i - a)*r + 0.5;
            array[j*xsize + i] |= code;
        }
    }
    else {   /* n < m */
        jstep = (b < d) ? +1 : -1;
        r = ((double)(c - a))/((double)(d - b));
        for (j=b; j!=d+jstep; j+=jstep) {
            i = a + (j - b)*r + 0.5;
            array[j*xsize + i] |= code;
        }
    }
    return;
}



/* Return the squared distance between the pixel (a,b) and the */
/* nearest border pixel.  The border pixels are encoded in the */
/* bitflags array by the value of "border_code".               */
int DistToBorder(unsigned char *bitflags,
                 int border_code,
                 int a,
                 int b,
                 int *ra,
                 int *rb, int xsize,
                 int ysize)
{
    int  besta, bestb, found, dist2, best_dist2;
    int  i, j, k, bs;
    *ra = *rb = 0;
    for (bs=0; bs<xsize + ysize; bs++) {
        found = 0;
        best_dist2 = 1000000;  /* initialize to large value */
        /* search boxes of increasing size until border pixel found */
        for (j=b - bs; j<=b + bs; j++) {
            for (i=a - bs; i<=a + bs; i++) {
                k = j*xsize + i;
                if (i<=0 || i>=xsize - 1 || j<=0 || j>=ysize - 1
                    || (bitflags[k] & border_code)) {
                    found = 1;
                    dist2 = (j - b)*(j - b) + (i - a)*(i - a);
                    if (dist2 < best_dist2) {
                        best_dist2 = dist2;
                        besta = i;
                        bestb = j;
                    }
                }
            }
        }
        if (found) {
            *ra = besta;
            *rb = bestb;
            break;
        }
    }
    return best_dist2;
}




/* Goldstein's phase-unwrapping algorithm.  The bitflags store */
/* the masked pixels (to be ignored) and the residues and      */
/* accumulates other info such as the branch cut pixels.       */
void BranchCuts_parallel(unsigned char *bitflags,
                         int MaxCutLen,
                         int NumRes,
                         int xsize,
                         int ysize,
                         int iniy,
                         int endy)
{
    int            i, j, k, ii, jj, kk, m, n, ri, rj;
    int            charge, boxctr_i, boxctr_j, boxsize, bs2;
    int            dist, min_dist, rim_i, rim_j, near_i, near_j;
    int            ka, num_active, max_active, *active_list;
    int            bench;
    int            draw_cut_line;
    double         r;

    if (MaxCutLen < 2) MaxCutLen = 2;
    max_active = NumRes + 10;
    AllocateInt(&active_list, max_active + 1, "book keeping data");

    /* branch cuts */

    for (j=iniy; j<endy; j++)
    {
        for (i=0; i<xsize; i++)
        {
        	k = j*xsize + i;

            if ((bitflags[k] & (POS_RES | NEG_RES))
                && !(bitflags[k] & VISITED))
            {
                bitflags[k] |= VISITED;  /* turn on visited flag */
                bitflags[k] |= ACTIVE;   /* turn on active flag */
                charge = (bitflags[k] & POS_RES) ? 1 : -1;
                num_active = 0;
                active_list[num_active++] = k;

                if (num_active > max_active)
                    num_active = max_active;

                for (boxsize = 3; boxsize<=2*MaxCutLen; boxsize += 2)
                {
                    bs2 = boxsize/2;
                    for (ka=0; ka<num_active; ka++)
                    {
                        boxctr_i = active_list[ka]%xsize;
                        boxctr_j = active_list[ka]/xsize;
                        for (jj=boxctr_j - bs2; jj<=boxctr_j + bs2; jj++)
                        {
                            for (ii=boxctr_i - bs2; ii<=boxctr_i + bs2; ii++)
                            {
                                kk = jj*xsize + ii;
                                if (ii<0 || ii>=xsize || jj<0 || jj>=ysize)
                                {
                                    continue;
                                }
                                else
                                {
                                    if (ii==0 || ii==xsize-1 || jj==0 || jj==ysize-1
                                        || (bitflags[kk] & BORDER))
                                    {
                                        charge = 0;
                                        DistToBorder(bitflags, BORDER, boxctr_i,
                                                     boxctr_j, &ri, &rj, xsize, ysize);
                                        PlaceCut(bitflags, ri, rj, boxctr_i, boxctr_j,
                                                 xsize, ysize, BRANCH_CUT);
                                    }
                                    else if ((bitflags[kk] & (POS_RES | NEG_RES))
                                             && !(bitflags[kk] & ACTIVE))
                                    {
                                        if (!(bitflags[kk] & VISITED))
                                        {
                                            charge += (bitflags[kk] & POS_RES) ? 1 : -1;
                                            bitflags[kk] |= VISITED;   /* set flag */
                                        }
                                        active_list[num_active++] = kk;
                                        if (num_active > max_active)
                                            num_active = max_active;
                                        bitflags[kk] |= ACTIVE;  /* set active flag */
                                        PlaceCut(bitflags, ii, jj, boxctr_i, boxctr_j,
                                                 xsize, ysize, BRANCH_CUT);
                                    }
                                    if (charge==0)
                                        goto continue_scan;
                                }  /* else */
                            }   /* for (ii ... */
                        }   /* for (jj ... */
                    }  /* for (ka ... */
                }   /* for (boxsize ... */

                if (charge != 0)
                {   /* connect branch cuts to rim */
                    min_dist = xsize + ysize;  /* large value */
                    for (ka=0; ka<num_active; ka++)
                    {
                        ii = active_list[ka]%xsize;
                        jj = active_list[ka]/xsize;
                        if ((dist = DistToBorder(bitflags, BORDER,
                                                 ii, jj, &ri, &rj, xsize, ysize))<min_dist)
                        {
                            min_dist = dist;
                            near_i = ii;
                            near_j = jj;
                            rim_i = ri;
                            rim_j = rj;
                        }
                    }

                    PlaceCut(bitflags, near_i, near_j, rim_i, rim_j,
                             xsize, ysize, BRANCH_CUT);
                }
                continue_scan :
                /* mark all active pixels inactive */
                for (ka=0; ka<num_active; ka++)
                    bitflags[active_list[ka]] &= ~ACTIVE;  /* turn flag off */
            }  /* if (bitflags ... */

        }
    }



    free(active_list);

    return;
}




/* Goldstein's phase-unwrapping algorithm.  The bitflags store */
/* the masked pixels (to be ignored) and the residues and      */
/* accumulates other info such as the branch cut pixels.       */
void GoldsteinBranchCuts_parallel(unsigned char *bitflags,
                         int MaxCutLen,
                         int NumRes,
                         int xsize,
                         int ysize)
{
    int band, b, iniy, endy;
    int MaxCutLen2;


    /* length of a band */
    band = ceil((double)ysize/(double)NUM_CORES);

    MaxCutLen2 = (xsize + band)/2;

    /*
     * Place branch cuts per band
     */
	#pragma omp parallel for default(none) \
	private(b, iniy, endy) \
	shared(NUM_CORES, band, bitflags, MaxCutLen2, NumRes, xsize, ysize)
    for (b=0; b<NUM_CORES; b++)
    {
    	iniy = b*band;

    	if (b<NUM_CORES-1)
    	    endy = iniy + band;
    	else
    	    endy = ysize;

    	BranchCuts_parallel(bitflags, MaxCutLen2, NumRes, xsize, ysize, iniy, endy);
    }
}







/* Goldstein's phase-unwrapping algorithm.  The bitflags store */
/* the masked pixels (to be ignored) and the residues and      */
/* accumulates other info such as the branch cut pixels.       */
void GoldsteinBranchCuts_serial(unsigned char *bitflags,
                         int MaxCutLen,
                         int NumRes,
                         int xsize,
                         int ysize)
{
    int            i, j, k, ii, jj, kk, m, n, ri, rj;
    int            charge, boxctr_i, boxctr_j, boxsize, bs2;
    int            dist, min_dist, rim_i, rim_j, near_i, near_j;
    int            ka, num_active, max_active, *active_list;
    int            bench;
    int            draw_cut_line;
    double         r;

    if (MaxCutLen < 2) MaxCutLen = 2;
    max_active = NumRes + 10;
    AllocateInt(&active_list, max_active + 1, "book keeping data");

    /* branch cuts */
    printf("Computing branch cuts\n");

    for (j=0; j<ysize; j++)
    {
        for (i=0; i<xsize; i++)
        {
            k = j*xsize + i;
            if ((bitflags[k] & (POS_RES | NEG_RES))
                && !(bitflags[k] & VISITED))
            {
                bitflags[k] |= VISITED;  /* turn on visited flag */
                bitflags[k] |= ACTIVE;   /* turn on active flag */
                charge = (bitflags[k] & POS_RES) ? 1 : -1;
                num_active = 0;
                active_list[num_active++] = k;

                if (num_active > max_active)
                    num_active = max_active;

                for (boxsize = 3; boxsize<=2*MaxCutLen; boxsize += 2)
                {
                    bs2 = boxsize/2;
                    for (ka=0; ka<num_active; ka++)
                    {
                        boxctr_i = active_list[ka]%xsize;
                        boxctr_j = active_list[ka]/xsize;
                        for (jj=boxctr_j - bs2; jj<=boxctr_j + bs2; jj++)
                        {
                            for (ii=boxctr_i - bs2; ii<=boxctr_i + bs2; ii++)
                            {
                                kk = jj*xsize + ii;
                                if (ii<0 || ii>=xsize || jj<0 || jj>=ysize)
                                {
                                    continue;
                                }
                                else
                                {
                                    if (ii==0 || ii==xsize-1 || jj==0 || jj==ysize-1
                                        || (bitflags[kk] & BORDER))
                                    {
                                        charge = 0;
                                        DistToBorder(bitflags, BORDER, boxctr_i,
                                                     boxctr_j, &ri, &rj, xsize, ysize);
                                        PlaceCut(bitflags, ri, rj, boxctr_i, boxctr_j,
                                                 xsize, ysize, BRANCH_CUT);
                                    }
                                    else if ((bitflags[kk] & (POS_RES | NEG_RES))
                                             && !(bitflags[kk] & ACTIVE))
                                    {
                                        if (!(bitflags[kk] & VISITED))
                                        {
                                            charge += (bitflags[kk] & POS_RES) ? 1 : -1;
                                            bitflags[kk] |= VISITED;   /* set flag */
                                        }
                                        active_list[num_active++] = kk;
                                        if (num_active > max_active)
                                            num_active = max_active;
                                        bitflags[kk] |= ACTIVE;  /* set active flag */
                                        PlaceCut(bitflags, ii, jj, boxctr_i, boxctr_j,
                                                 xsize, ysize, BRANCH_CUT);
                                    }
                                    if (charge==0)
                                        goto continue_scan;
                                }  /* else */
                            }   /* for (ii ... */
                        }   /* for (jj ... */
                    }  /* for (ka ... */
                }   /* for (boxsize ... */

                if (charge != 0)
                {   /* connect branch cuts to rim */
                    min_dist = xsize + ysize;  /* large value */
                    for (ka=0; ka<num_active; ka++)
                    {
                        ii = active_list[ka]%xsize;
                        jj = active_list[ka]/xsize;
                        if ((dist = DistToBorder(bitflags, BORDER,
                                                 ii, jj, &ri, &rj, xsize, ysize))<min_dist)
                        {
                            min_dist = dist;
                            near_i = ii;
                            near_j = jj;
                            rim_i = ri;
                            rim_j = rj;
                        }
                    }

                    PlaceCut(bitflags, near_i, near_j, rim_i, rim_j,
                             xsize, ysize, BRANCH_CUT);
                }
                continue_scan :
                /* mark all active pixels inactive */
                for (ka=0; ka<num_active; ka++)
                    bitflags[active_list[ka]] &= ~ACTIVE;  /* turn flag off */
            }  /* if (bitflags ... */
        }  /* for (i ... */
    }  /* for (j ... */


    free(active_list);

    return;
}




/* Detect residues in phase data and mark them as positive or  */
/* negative residues in the bitflags array.  Ignore the pixels */
/* marked with avoid_code in the bitflags araay.               */

int Residues_parallel(float *phase,
             unsigned char *bitflags,
             int xsize,
             int ysize)
{
    int  i, j, k, NumRes=0;
    double  r;

    #pragma omp parallel for default(none) \
    private(j, i, k, r) \
    shared(ysize, xsize, bitflags, phase) \
    reduction(+ : NumRes) \
    collapse(2)
    for (j=0; j<ysize - 1; j++)
    {
        for (i=0; i<xsize - 1; i++)
        {
            k = j*xsize + i;

            if (bitflags && ((bitflags[k] & AVOID)
                             || (bitflags[k+1] & AVOID)
                             || (bitflags[k+1+xsize] & AVOID)
                             || (bitflags[k+xsize] & AVOID))) {
                continue; /* masked region: don't unwrap */
            }
            r = Gradient(phase[k+1], phase[k])
            + Gradient(phase[k+1+xsize], phase[k+1])
            + Gradient(phase[k+xsize], phase[k+1+xsize])
            + Gradient(phase[k], phase[k+xsize]);
            if (bitflags) {
                if (r > 0.01) bitflags[k] |= POS_RES;
                else if (r < -0.01) bitflags[k] |= NEG_RES;
            }
            if (r*r > 0.01)
                ++NumRes;
        }
    }
    return NumRes;
}



/* Detect residues in phase data and mark them as positive or  */
/* negative residues in the bitflags array.  Ignore the pixels */
/* marked with avoid_code in the bitflags araay.               */

int Residues_serial(float *phase,
             unsigned char *bitflags,
             int xsize,
             int ysize)
{
    int  i, j, k, NumRes=0;
    double  r;

    for (j=0; j<ysize - 1; j++)
    {
        for (i=0; i<xsize - 1; i++)
        {
            k = j*xsize + i;

            if (bitflags && ((bitflags[k] & AVOID)
                             || (bitflags[k+1] & AVOID)
                             || (bitflags[k+1+xsize] & AVOID)
                             || (bitflags[k+xsize] & AVOID))) {
                continue; /* masked region: don't unwrap */
            }
            r = Gradient(phase[k+1], phase[k])
            + Gradient(phase[k+1+xsize], phase[k+1])
            + Gradient(phase[k+xsize], phase[k+1+xsize])
            + Gradient(phase[k], phase[k+xsize]);
            if (bitflags) {
                if (r > 0.01) bitflags[k] |= POS_RES;
                else if (r < -0.01) bitflags[k] |= NEG_RES;
            }
            if (r*r > 0.01)
                ++NumRes;
        }
    }
    return NumRes;
}






double goldstein_phase_unwrapping(const char *pname,
								  const char *data_path,
                                  int type,
                                  int xsize,
                                  int ysize,
                                  int mask_flag)
{
    /* data variables */

    int           *path_order;
    float         *phase;
    float         *soln;
    float         *grady, *gradx;
    float         *mask;
    unsigned char *unwrap, *bitflags;
    clock_t       t1, t2;
    double        elapsed_time;


    /* other variables */


    FILE          *ifp, *ofp, *ifq, *ifm;
    char          prefix[120], fname[120];
    float         grad, minval, maxval, high_qual, valf, aux;
    int           k, length, index, p, q, iter=0, bin, l;
    int           a, b, w, x, y, i, j, num_pieces;
    int           imin, imax, binarg, NumRes, MaxCutLen;
    int           *list;

	

    /*
     *     ALLOCATE MEMORY
     */

    length = xsize*ysize;

    AllocateFloat(&phase, length, "phase data");
    AllocateFloat(&soln, length, "solution array");
    AllocateFloat(&grady, length, "vertical gradient");
    AllocateFloat(&gradx, length, "horizontal gradient");
    AllocateByte(&unwrap, length, "flag array");
    AllocateByte(&bitflags, length, "flag array");
    AllocateInt(&path_order, length, "integration path");
    AllocateFloat(&mask, length, "mask array");
    AllocateInt(&list, 2*(xsize+ysize), "in-out list");


    /*
     *    READ MASK
     */

    if (mask_flag) {
		strcpy(prefix, data_path);
		strcat(prefix,"\\data\\");
		strcat(prefix,pname);
		strcat(prefix, ".mask");
		
        OpenFile(&ifm, prefix, "rb");
        GetPhase(type, ifm, fname, mask, xsize, ysize);

        for (k=0; k<length; k++) {
            if (mask[k]>0)
                mask[k] = 1;
            else
                mask[k] = 0;
        }
    }
    else
    {
        for (k=0; k<length; k++)
            mask[k] = 1;
    }




    /*
     *    READ PHASE
     */
    prefix[0]='\0';
	strcpy(prefix, data_path);
	strcat(prefix,"\\data\\");
	strcat(prefix,pname);
	strcat(prefix, ".phase");
	
    OpenFile(&ifp, prefix, "rb");
    GetPhase(type, ifp, fname, phase, xsize, ysize);


    //for (k=0; k<length; k++) phase[k] = phase[k]*TWOPI;


    /*
     *    INIT BITFLAGS
     */
	#pragma omp parallel for default(none) \
	private(k) \
	shared(length, mask, bitflags)
    for (k=0; k<length; k++) {
        if (mask[k]==0)
            bitflags[k] |= BORDER;
        else
            bitflags[k] = 0;
    }


    /* Pre-compute the vertical and horizontal gradients */
    Gradxy(phase, gradx, grady, xsize, ysize);


    // ** starting time

    t1 = clock();

    /*
     *    LOCATE AND PROCESS RESIDUES
     */

    //NumRes = Residues_serial(phase, bitflags, xsize, ysize);
    NumRes = Residues_parallel(phase, bitflags, xsize, ysize);


    printf("Number of residues: %d\n", NumRes);
    
	// Save residues image
	prefix[0]='\0';
	strcpy(prefix, data_path);
	strcat(prefix,"\\data\\");
	strcat(prefix,pname);
	strcat(prefix, ".res");
    SaveByteToImage(bitflags, "residues", prefix, xsize, ysize, 1, 1, 0);


    /*
     *    GENERATE BRANCH CUTS
     */

    MaxCutLen = (xsize + ysize)/2;

    //GoldsteinBranchCuts_serial(bitflags, MaxCutLen, NumRes, xsize, ysize);
    GoldsteinBranchCuts_parallel(bitflags, MaxCutLen, NumRes, xsize, ysize);


	// Save branch cuts image
	
	prefix[0]='\0';
	strcpy(prefix, data_path);
	strcat(prefix,"\\data\\");
	strcat(prefix,pname);
	strcat(prefix, ".brc");
    SaveByteToImage(bitflags, "branch cuts", prefix, xsize, ysize, 1, 1, BRANCH_CUT | BORDER);



    /*
     *    UNWRAP AROUND CUTS
     */

    //num_pieces = UnwrapAroundCutsGoldstein(phase,bitflags, soln, xsize, ysize, path_order);
    num_pieces = UnwrapAroundCutsFrontier(phase,bitflags, soln, xsize, ysize, path_order, grady, gradx, list, length);


    // ** end time

    t2 = clock();

    printf("Number of pieces: %d\n", num_pieces);

    // compute and print the elapsed time in millisec
    elapsed_time = timediff(t1, t2);

    printf("Elapsed time: %f ms\n", elapsed_time);


    /*
     *    SAVE RESULTS
     */

    for (k=0; k<length; k++)
        soln[k] *= TWOPI;


    /* save solution */
	
	strcpy(prefix, data_path);
    strcat(prefix, "\\data\\");
    strcat(prefix, pname);
	strcat(prefix, ".out");
    OpenFile(&ofp, prefix, "w");
    WriteFloat(ofp, soln, length, fname);
	
	
	/* save path inegration */
	strcpy(prefix, data_path);
    strcat(prefix, "\\data\\");
    strcat(prefix, pname);
    strcat(prefix, ".path");
    SaveIntToImage(path_order, "path integration", prefix, xsize, ysize);	
	printf("\n");
	
	
	/* DEALLOCATE MEMORY */

    free(phase);
    free(soln);
    free(unwrap);
    free(path_order);
    free(grady);
    free(gradx);
    free(list);
    free(bitflags);
    free(mask);


    return elapsed_time;
}



int main()
{
    int    mask_flag       = 0;
    int    type            = 3;
    int    MAX_ITERATIONS  = 2;
    double elapsed_time    = 0;
	char   data_path[1024];
    
	
    chdir("..");
    getcwd(data_path,1024);
	

    NUM_CORES = omp_get_num_procs()/2;

    /* set number of threads */
    omp_set_num_threads(NUM_CORES);


    printf("Number of threads: %d\n",NUM_CORES);

    for (int i=0; i<MAX_ITERATIONS; i++)
        elapsed_time += goldstein_phase_unwrapping("peaks.1024x1024", data_path, type, 1024, 1024, mask_flag);

    elapsed_time /= (double)MAX_ITERATIONS;
    printf("\nAverage elapsed time: %f ms\n", elapsed_time);

    return 0;
}
