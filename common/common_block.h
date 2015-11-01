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

#if !defined(_COMMON_BLOCK_H_)
#define _COMMON_BLOCK_H_

#include "types.h"
#include "simd.h"

#include "daala_compat.h"
#define OD_DERING_NBLOCKS (OD_BSIZE_MAX/8)

int get_left_available(int ypos, int xpos, int size, int width);
int get_up_available(int ypos, int xpos, int size, int width);
int get_upright_available(int ypos, int xpos, int size, int width);
int get_downleft_available(int ypos, int xpos, int size, int height);

void dequantize (int16_t *coeff,int16_t *rcoeff,int quant,int size);
void reconstruct_block(int16_t *block, uint8_t *pblock, uint8_t *rec, int size, int stride);

void find_block_contexts(int ypos, int xpos, int height, int width, int size, deblock_data_t *deblock_data, block_context_t *block_context, int enable);

void clpf_block(const uint8_t *src, uint8_t *dst, int sstride, int dstride, int x0, int y0, int size, int width, int height);

typedef void (*od_filter_dering_direction_func)(int16_t *y, int ystride,
 int16_t *in, int threshold, int dir);
typedef void (*od_filter_dering_orthogonal_func)(uint8_t *y, int ystride,
 int16_t *in, uint8_t *x, int xstride, int threshold, int dir);

void od_filter_dering_direction_4x4_c(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir);
void od_filter_dering_direction_8x8_c(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir);
void od_filter_dering_orthogonal_4x4_c(uint8_t *y, int ystride, int16_t *in,
 uint8_t *x, int xstride, int threshold, int dir);
void od_filter_dering_orthogonal_8x8_c(uint8_t *y, int ystride, int16_t *in,
 uint8_t *x, int xstride, int threshold, int dir);

#define OD_DERING_NBLOCKS (OD_BSIZE_MAX/8)

#define OD_FILT_BORDER (3)
#define OD_FILT_BSTRIDE (OD_BSIZE_MAX + 2*OD_FILT_BORDER)

extern const od_filter_dering_direction_func
 OD_DERING_DIRECTION_C[OD_DERINGSIZES];
extern const od_filter_dering_orthogonal_func
 OD_DERING_ORTHOGONAL_C[OD_DERINGSIZES];

void od_dering(uint8_t *y, int ystride, uint8_t *x, int
 xstride, int ln, int sbx, int sby, int nhsb, int nvsb, int q, int xdec,
 int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
               int pli, unsigned char *bskip, int skip_stride);

#endif
