/*
Copyright (c) 2015, Cisco Systems
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <memory.h>
#include <assert.h>

#include "global.h"
#include "common_block.h"

int zigzag16[16] = {
    0, 1, 5, 6, 
    2, 4, 7, 12, 
    3, 8, 11, 13, 
    9, 10, 14, 15
};

int zigzag64[64] = {
     0,  1,  5,  6, 14, 15, 27, 28,
     2,  4,  7, 13, 16, 26, 29, 42,
     3,  8, 12, 17, 25, 30, 41, 43,
     9, 11, 18, 24, 31, 40, 44, 53,
    10, 19, 23, 32, 39, 45, 52, 54,
    20, 22, 33, 38, 46, 51, 55, 60,
    21, 34, 37, 47, 50, 56, 59, 61,
    35, 36, 48, 49, 57, 58, 62, 63
};

int zigzag256[256] = {
    0,  1,  5,  6, 14, 15, 27, 28, 44, 45, 65, 66, 90, 91,119,120,
    2,  4,  7, 13, 16, 26, 29, 43, 46, 64, 67, 89, 92,118,121,150,
    3,  8, 12, 17, 25, 30, 42, 47, 63, 68, 88, 93,117,122,149,151,
    9, 11, 18, 24, 31, 41, 48, 62, 69, 87, 94,116,123,148,152,177,
   10, 19, 23, 32, 40, 49, 61, 70, 86, 95,115,124,147,153,176,178,
   20, 22, 33, 39, 50, 60, 71, 85, 96,114,125,146,154,175,179,200,
   21, 34, 38, 51, 59, 72, 84, 97,113,126,145,155,174,180,199,201,
   35, 37, 52, 58, 73, 83, 98,112,127,144,156,173,181,198,202,219,
   36, 53, 57, 74, 82, 99,111,128,143,157,172,182,197,203,218,220,
   54, 56, 75, 81,100,110,129,142,158,171,183,196,204,217,221,234,
   55, 76, 80,101,109,130,141,159,170,184,195,205,216,222,233,235,
   77, 79,102,108,131,140,160,169,185,194,206,215,223,232,236,245,
   78,103,107,132,139,161,168,186,193,207,214,224,231,237,244,246,
  104,106,133,138,162,167,187,192,208,213,225,230,238,243,247,252,
  105,134,137,163,166,188,191,209,212,226,229,239,242,248,251,253,
  135,136,164,165,189,190,210,211,227,228,240,241,249,250,254,255
};




int chroma_qp[52] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
        17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 29,
        30, 31, 32, 33, 33, 34, 34, 35, 35, 36, 36, 37, 37, 38,
        39, 40, 41, 42, 43, 44, 45
};

int super_table[8][20] = {
  {-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,  1, 0, 5, 2, 6, 3, 7, 4, 8,-1},
  {-1, 0, -1,-1,-1,-1,-1,-1,-1,-1,  2, 1, 6, 3, 7, 5, 8, 4, 9,-1},
  {-1, 0, -1,-1,-1,-1,-1,-1,-1,-1,  2, 1, 6, 3, 7, 5, 8, 4, 9,-1},
  {-1, 0, -1,-1,-1,-1,-1,-1,-1,-1,  2, 1, 6, 3, 7, 5, 8, 4, 9,-1},

  { 0,-1,  2, 1,12, 7,13, 5,16,11,  3, 4,14, 8, 9, 6,15,10,17,18},
  { 0, 1,  3, 2,10, 7,11, 6,16, 9,  5, 4,15,13,14, 8,17,12,18,19},
  { 0, 1,  3, 2,10, 4,12, 5,14, 6,  8, 7,15,13,16,11,17, 9,18,19},
  { 0, 1,  3, 2, 7, 4, 8, 5, 9, 6, 11,10,15,14,16,13,17,12,18,19}
};

const uint16_t gquant_table[6] = {26214,23302,20560,18396,16384,14564};
const uint16_t gdequant_table[6] = {40,45,51,57,64,72};

int get_left_available(int ypos, int xpos, int size, int width){
  int left_available = xpos > 0;
  return left_available;
}

int get_up_available(int ypos, int xpos, int size, int width){
  int up_available = ypos > 0;
  return up_available;
}

int get_upright_available(int ypos, int xpos, int size, int width){

  int upright_available = (ypos > 0) && (xpos + size < width);
  if (size==32 && (ypos%64)==32) upright_available = 0;
  if (size==16 && ((ypos%32)==16 || ((ypos%64)==32 && (xpos%32)==16))) upright_available = 0;
  if (size== 8 && ((ypos%16)==8 || ((ypos%32)==16 && (xpos%16)==8) || ((ypos%64)==32 && (xpos%32)==24))) upright_available = 0;

  return upright_available;
}

int get_downleft_available(int ypos, int xpos, int size, int height){

  int downleft_available = (xpos > 0) && (ypos + size < height);
  if (size==64) downleft_available = 0;
  if (size==32 && (ypos%64)==32) downleft_available = 0;
  if (size==16 && ((ypos%64)==48 || ((ypos%64)==16 && (xpos%32)==16))) downleft_available = 0;
  if (size== 8 && ((ypos%64)==56 || ((ypos%16)==8 && (xpos%16)==8) || ((ypos%64)==24 && (xpos%32)==16))) downleft_available = 0;

  return downleft_available;
}


void dequantize (int16_t *coeff, int16_t *rcoeff, int qp, int size)
{
  int tr_log2size = log2i(size);
  const int lshift = qp / 6;
  const int rshift = tr_log2size - 1;
  const int scale = gdequant_table[qp % 6];
  const int add = 1<<(rshift-1);

  for (int i = 0; i < size ; i++){
    for (int j = 0; j < size; j++){
      int c = coeff[i*size+j];
      rcoeff[i*size+j] = ((c * scale << lshift) + add) >> rshift;
    }
  }
}

void reconstruct_block(int16_t *block, uint8_t *pblock, uint8_t *rec, int size, int stride)
{ 
  int i,j;
  for(i=0;i<size;i++){    
    for (j=0;j<size;j++){
      rec[i*stride+j] = (uint8_t)clip255(block[i*size+j] + (int16_t)pblock[i*size+j]);      
    }
  }
}

void find_block_contexts(int ypos, int xpos, int height, int width, int size, deblock_data_t *deblock_data, block_context_t *block_context, int enable){

  if (ypos >= MIN_BLOCK_SIZE && xpos >= MIN_BLOCK_SIZE && ypos + size < height && xpos + size < width && enable && size <= MAX_TR_SIZE) {
    int by = ypos/MIN_PB_SIZE;
    int bx = xpos/MIN_PB_SIZE;
    int bs = width/MIN_PB_SIZE;
    int bindex = by*bs+bx;
    block_context->split = (deblock_data[bindex-bs].size < size) + (deblock_data[bindex-1].size < size);
    int cbp1;
    cbp1 = (deblock_data[bindex-bs].cbp.y > 0) + (deblock_data[bindex-1].cbp.y > 0);
    block_context->cbp = cbp1;
    int cbp2 = (deblock_data[bindex-bs].cbp.y > 0 || deblock_data[bindex-bs].cbp.u > 0 || deblock_data[bindex-bs].cbp.v > 0) +
           (deblock_data[bindex-1].cbp.y > 0 || deblock_data[bindex-1].cbp.u > 0 || deblock_data[bindex-1].cbp.v > 0);
    block_context->index = 3*block_context->split + cbp2;
  }
  else{
    block_context->split = -1;
    block_context->cbp = -1;
    block_context->index = -1;
  }
}

void clpf_block(const uint8_t *src, uint8_t *dst, int sstride, int dstride, int x0, int y0, int size, int width, int height) {
  int left = x0 & ~(dstride-1);
  int top = y0 & ~(dstride-1);
  int right = min(width-1, left + dstride-1);
  int bottom = min(height-1, top + dstride-1);

  for (int y=y0;y<y0+size;y++){
    for (int x=x0;x<x0+size;x++) {
      int X = src[(y+0)*sstride + x+0];
      int A = y == top ? X : src[(y-1)*sstride + x+0];
      int B = x == left ? X : src[(y+0)*sstride + x-1];
      int C = x == right ? X : src[(y+0)*sstride + x+1];
      int D = y == bottom ? X : src[(y+1)*sstride + x+0];
      int delta = ((A>X)+(B>X)+(C>X)+(D>X) > 2) - ((A<X)+(B<X)+(C<X)+(D<X) > 2);
      dst[(y-top)*dstride + x-left] = X + delta;
    }
  }
}

const od_filter_dering_direction_func
 OD_DERING_DIRECTION_C[OD_DERINGSIZES] = {
  od_filter_dering_direction_4x4_c,
  od_filter_dering_direction_8x8_c
};

const od_filter_dering_orthogonal_func
 OD_DERING_ORTHOGONAL_C[OD_DERINGSIZES] = {
  od_filter_dering_orthogonal_4x4_c,
  od_filter_dering_orthogonal_8x8_c
};

int direction_offsets_table[16][3] = {
  { -69,-138,-207 },
  {   1, -68, -67 },
  {   1,   2,   3 },
  {   1,  72,  73 },
  {  71, 142, 213 },
  {  70, 141, 211 },
  {  70, 140, 210 },
  {  70, 139, 209 },
  {  69, 138, 207 },
  {  69, 137, 206 },
  {  68, 136, 204 },
  {  68, 135, 203 },
  {  67, 134, 201 },
  {  67, 133, 200 },
  {  66, 132, 198 },
  {  66, 131, 197 }
};

/* Detect direction. 0 means 45-degree up-right, 2 is horizontal, and so on.
   The search minimizes the weighted variance along all the lines in a
   particular direction, i.e. the squared error between the input and a
   "predicted" block where each pixel is replaced by the average along a line
   in a particular direction. Since each direction have the same sum(x^2) term,
   that term is never computed. See Section 2, step 2, of:
   http://jmvalin.ca/notes/intra_paint.pdf */
static int od_dir_find8(const uint8_t *img, int stride, int32_t *var) {
  int i;
  int cost[8] = {0};
  int partial[8][15] = {{0}};
  int best_cost = 0;
  int best_dir = 0;
  for (i = 0; i < 8; i++) {
    int j;
    for (j = 0; j < 8; j++) {
      int x;
      /*x = img[i*stride + j] >> OD_COEFF_SHIFT;*/
      x = img[i*stride + j];
      partial[0][i + j] += x;
      partial[1][i + j/2] += x;
      partial[2][i] += x;
      partial[3][3 + i - j/2] += x;
      partial[4][7 + i - j] += x;
      partial[5][3 - i/2 + j] += x;
      partial[6][j] += x;
      partial[7][i/2 + j] += x;
    }
  }
  for (i = 0; i < 8; i++) {
    cost[2] += partial[2][i]*partial[2][i] >> 3;
    cost[6] += partial[6][i]*partial[6][i] >> 3;
  }
  for (i = 0; i < 7; i++) {
    cost[0] += OD_DIVU_SMALL(partial[0][i]*partial[0][i], i + 1)
     + OD_DIVU_SMALL(partial[0][14 - i]*partial[0][14 - i], i + 1);
    cost[4] += OD_DIVU_SMALL(partial[4][i]*partial[4][i], i + 1)
     + OD_DIVU_SMALL(partial[4][14 - i]*partial[4][14 - i], i + 1);
  }
  cost[0] += partial[0][7]*partial[0][8 - 1] >> 3;
  cost[4] += partial[4][7]*partial[4][8 - 1] >> 3;
  for (i = 1; i < 8; i += 2) {
    int j;
    for (j = 0; j < 4 + 1; j++) {
      cost[i] += partial[i][3 + j]*partial[i][3 + j] >> 3;
    }
    for (j = 0; j < 4 - 1; j++) {
      cost[i] += OD_DIVU_SMALL(partial[i][j]*partial[i][j], 2*j + 2)
       + OD_DIVU_SMALL(partial[i][10 - j]*partial[i][10 - j], 2*j + 2);
    }
  }
  for (i = 0; i < 8; i++) {
    if (cost[i] > best_cost) {
      best_cost = cost[i];
      best_dir = i;
    }
  }
  /* Difference between the optimal variance and the variance along the
     orthogonal direction. Again, the sum(x^2) terms cancel out. */
  *var = best_cost - cost[(best_dir + 4) & 7];
  return best_dir;
}

#define OD_DERING_VERY_LARGE (30000)
#define OD_DERING_INBUF_SIZE ((OD_BSIZE_MAX + 2*OD_FILT_BORDER)*\
 (OD_BSIZE_MAX + 2*OD_FILT_BORDER))

/* Smooth in the direction detected. */
void od_filter_dering_direction_c(int16_t *y, int ystride, int16_t *in,
 int ln, int threshold, int dir) {
  int i;
  int j;
  int k;
  static const int taps[3] = {3, 2, 2};
  for (i = 0; i < 1 << ln; i++) {
    for (j = 0; j < 1 << ln; j++) {
      od_coeff sum;
      od_coeff xx;
      od_coeff yy;
      xx = in[i*OD_FILT_BSTRIDE + j];
      sum= 0;
      for (k = 0; k < 3; k++) {
        od_coeff p0;
        od_coeff p1;
        p0 = in[i*OD_FILT_BSTRIDE + j + direction_offsets_table[dir][k]] - xx;
        p1 = in[i*OD_FILT_BSTRIDE + j - direction_offsets_table[dir][k]] - xx;
        if (abs(p0) < threshold) sum += taps[k]*p0;
        if (abs(p1) < threshold) sum += taps[k]*p1;
      }
      yy = xx + ((sum + 8) >> 4);
      y[i*ystride + j] = yy;
    }
  }
}

void od_filter_dering_direction_4x4_c(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir) {
  od_filter_dering_direction_c(y, ystride, in, 2, threshold, dir);
}

void od_filter_dering_direction_8x8_c(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir) {
  od_filter_dering_direction_c(y, ystride, in, 3, threshold, dir);
}

/* Smooth in the direction orthogonal to what was detected. */
void od_filter_dering_orthogonal_c(uint8_t *y, int ystride, int16_t *in,
 uint8_t *x, int xstride, int ln, int threshold, int dir) {
  int i;
  int j;
  int offset;
  if (dir <= 4) offset = OD_FILT_BSTRIDE;
  else offset = 1;
  for (i = 0; i < 1 << ln; i++) {
    for (j = 0; j < 1 << ln; j++) {
      od_coeff athresh;
      od_coeff yy;
      od_coeff sum;
      od_coeff p;
      /* Deringing orthogonal to the direction uses a tighter threshold
         because we want to be conservative. We've presumably already
         achieved some deringing, so the amount of change is expected
         to be low. Also, since we might be filtering across an edge, we
         want to make sure not to blur it. That being said, we might want
         to be a little bit more aggressive on pure horizontal/vertical
         since the ringing there tends to be directional, so it doesn't
         get removed by the directional filtering. */
      athresh = OD_MINI(threshold, threshold/3
       + abs(in[i*OD_FILT_BSTRIDE + j] - x[i*xstride + j]));
      yy = in[i*OD_FILT_BSTRIDE + j];
      sum = 0;
      p = in[i*OD_FILT_BSTRIDE + j + offset] - yy;
      if (abs(p) < athresh) sum += p;
      p = in[i*OD_FILT_BSTRIDE + j - offset] - yy;
      if (abs(p) < athresh) sum += p;
      p = in[i*OD_FILT_BSTRIDE + j + 2*offset] - yy;
      if (abs(p) < athresh) sum += p;
      p = in[i*OD_FILT_BSTRIDE + j - 2*offset] - yy;
      if (abs(p) < athresh) sum += p;
      y[i*ystride + j] = (yy + ((3*sum + 8) >> 4) + (1 << OD_COEFF_SHIFT >> 1)) >> OD_COEFF_SHIFT;
    }
  }
}

void od_filter_dering_orthogonal_4x4_c(uint8_t *y, int ystride, int16_t *in,
 uint8_t *x, int xstride, int threshold, int dir) {
  od_filter_dering_orthogonal_c(y, ystride, in, x, xstride, 2, threshold, dir);
}

void od_filter_dering_orthogonal_8x8_c(uint8_t *y, int ystride, int16_t *in,
 uint8_t *x, int xstride, int threshold, int dir) {
  od_filter_dering_orthogonal_c(y, ystride, in, x, xstride, 3, threshold, dir);
}

/* This table approximates x^0.16 with the index being log2(x). It is clamped
   to [-.5, 3]. The table is computed as:
   round(256*min(3, max(.5, 1.08*(sqrt(2)*2.^([0:17]+8)/256/256).^.16))) */
static int16_t od_thresh_table_q8[18] = {
  128, 134, 150, 168, 188, 210, 234, 262,
  292, 327, 365, 408, 455, 509, 569, 635,
  710, 768,
};

/* Compute deringing filter threshold for each 8x8 block based on the
   directional variance difference. A high variance difference means that we
   have a highly directional pattern (e.g. a high contrast edge), so we can
   apply more deringing. A low variance means that we either have a low
   contrast edge, or a non-directional texture, so we want to be careful not
   to blur. */
static void od_compute_thresh(int thresh[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
 int threshold, int32_t var[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
 int32_t sb_var, int nhb, int nvb) {
  int bx;
  int by;
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      int v1;
      int v2;
      /* We use both the variance of 8x8 blocks and the variance of the
         entire superblock to determine the threshold. */
      v1 = OD_MINI(32767, var[by][bx] >> 6);
      v2 = OD_MINI(32767, sb_var/(OD_BSIZE_MAX*OD_BSIZE_MAX));
      thresh[by][bx] = threshold*od_thresh_table_q8[OD_CLAMPI(0,
       OD_ILOG(v1*v2) - 9, 17)] >> 8;
    }
  }
}

void od_dering(uint8_t *y, int ystride, uint8_t *x, int
 xstride, int ln, int sbx, int sby, int nhsb, int nvsb, int q, int xdec,
 int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
               int pli, unsigned char *bskip, int skip_stride) {
  int i;
  int j;
  int n;
  int threshold;
  int bx;
  int by;
  int16_t inbuf[OD_DERING_INBUF_SIZE];
  int16_t *in;
  int16_t tmp[64*64];
  int tmp_stride = 64;
  int nhb;
  int nvb;
  int bsize;
  int varsum = 0;
  int32_t var[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS];
  int thresh[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS];

  int dering_thresholds[56] = {
    17, 17, 17, 17, 18, 19, 20, 23,
    25, 28, 31, 34, 37, 40, 44, 48,
    53, 58, 64, 70, 77, 85, 93, 102,
    112, 123, 135, 148, 162, 178, 195, 214,
    235, 257, 282, 309, 340, 372, 408, 448,
    492, 539, 591, 649, 712, 781, 856, 939,
    1031, 1130, 1240, 1360, 1492, 1637, 1795, 1969
  };

  n = 1 << ln;
  bsize = 3 - xdec;
  nhb = nvb = n >> bsize;
  in = inbuf + OD_FILT_BORDER*OD_FILT_BSTRIDE + OD_FILT_BORDER;
  /* We avoid filtering the pixels for which some of the pixels to average
     are outside the frame. We could change the filter instead, but it would
     add special cases for any future vectorization. */
  for (i = 0; i < OD_DERING_INBUF_SIZE; i++) inbuf[i] = OD_DERING_VERY_LARGE;
  for (i = -OD_FILT_BORDER*(sby != 0); i < n
   + OD_FILT_BORDER*(sby != nvsb - 1); i++) {
    for (j = -OD_FILT_BORDER*(sbx != 0); j < n
     + OD_FILT_BORDER*(sbx != nhsb - 1); j++) {
      in[i*OD_FILT_BSTRIDE + j] = x[i*xstride + j] << OD_COEFF_SHIFT;
    }
  }
  threshold = dering_thresholds[q];
  if (pli == 0) {
    for (by = 0; by < nvb; by++) {
      for (bx = 0; bx < nhb; bx++) {
        dir[by][bx] = od_dir_find8(&x[8*by*xstride + 8*bx], xstride,
         &var[by][bx]);
        varsum += var[by][bx];
      }
    }
    od_compute_thresh(thresh, threshold, var, varsum, nhb, nvb);
  }
  else {
    for (by = 0; by < nvb; by++) {
      for (bx = 0; bx < nhb; bx++) {
        thresh[by][bx] = threshold;
      }
    }
  }
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      if (0 && bskip[by*skip_stride + bx])
        thresh[by][bx] = 0;
    }
  }
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      OD_DERING_DIRECTION_C[bsize - 2](
       &tmp[(by*tmp_stride << bsize) + (bx << bsize)], tmp_stride,
       &in[(by*OD_FILT_BSTRIDE << bsize) + (bx << bsize)],
       thresh[by][bx], dir[by][bx]);
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      in[i*OD_FILT_BSTRIDE + j] = tmp[i*tmp_stride + j];
    }
  }
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      OD_DERING_ORTHOGONAL_C[bsize - 2](
       &y[(by*ystride << bsize) + (bx << bsize)], ystride,
       &in[(by*OD_FILT_BSTRIDE << bsize) + (bx << bsize)],
       &x[(by*xstride << bsize) + (bx << bsize)], xstride,
       thresh[by][bx], dir[by][bx]);
    }
  }
}
