/*
 *      ANSI C implementation of vector operations.
 *
 * Copyright (c) 2007-2010 Naoaki Okazaki
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 *This file has been modified from the original arithmetic_ansi.h file so
 *that it uses mxMalloc and mxFree if MATLAB_MEX_FILE is defined, making
 *the integration of the library into Matlab simpler.
 *David F. Crouse 6 August 2015
 */

/* $Id$ */

#include <stdlib.h>
#include <memory.h>

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#endif

#if     LBFGS_FLOAT == 32 && LBFGS_IEEE_FLOAT
#define fsigndiff(x, y) (((*(uint32_t*)(x)) ^ (*(uint32_t*)(y))) & 0x80000000U)
#else
#define fsigndiff(x, y) (*(x) * (*(y) / fabs(*(y))) < 0.)
#endif/*LBFGS_IEEE_FLOAT*/

inline static void* vecalloc(size_t size)
{
    #ifdef MATLAB_MEX_FILE
    void *memblock = mxMalloc(size);
    #else
    void *memblock = malloc(size);
    #endif
    
    if (memblock) {
        memset(memblock, 0, size);
    }
    return memblock;
}

inline static void vecfree(void *memblock)
{
    #ifdef MATLAB_MEX_FILE
    mxFree(memblock);
    #else
    free(memblock);
    #endif
}

inline static void vecset(lbfgsfloatval_t *x, const lbfgsfloatval_t c, const int n)
{
    int i;
    
    for (i = 0;i < n;++i) {
        x[i] = c;
    }
}

inline static void veccpy(lbfgsfloatval_t *y, const lbfgsfloatval_t *x, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] = x[i];
    }
}

inline static void vecncpy(lbfgsfloatval_t *y, const lbfgsfloatval_t *x, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] = -x[i];
    }
}

inline static void vecadd(lbfgsfloatval_t *y, const lbfgsfloatval_t *x, const lbfgsfloatval_t c, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] += c * x[i];
    }
}

inline static void vecdiff(lbfgsfloatval_t *z, const lbfgsfloatval_t *x, const lbfgsfloatval_t *y, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        z[i] = x[i] - y[i];
    }
}

inline static void vecscale(lbfgsfloatval_t *y, const lbfgsfloatval_t c, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] *= c;
    }
}

inline static void vecmul(lbfgsfloatval_t *y, const lbfgsfloatval_t *x, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] *= x[i];
    }
}

inline static void vecdot(lbfgsfloatval_t* s, const lbfgsfloatval_t *x, const lbfgsfloatval_t *y, const int n)
{
    int i;
    *s = 0.;
    for (i = 0;i < n;++i) {
        *s += x[i] * y[i];
    }
}

inline static void vec2norm(lbfgsfloatval_t* s, const lbfgsfloatval_t *x, const int n)
{
    vecdot(s, x, x, n);
    *s = (lbfgsfloatval_t)sqrt(*s);
}

inline static void vec2norminv(lbfgsfloatval_t* s, const lbfgsfloatval_t *x, const int n)
{
    vec2norm(s, x, n);
    *s = (lbfgsfloatval_t)(1.0 / *s);
}
