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

#include "global.h"
#include <string.h>

#include "mainenc.h"
#include "encode_block.h"
#include "common_block.h"
#include "common_frame.h"
#include "enc_kernels.h"

#include "daala_compat.h"

extern int chroma_qp[52];
const double squared_lambda_QP [52] = {
    0.0382, 0.0485, 0.0615, 0.0781, 0.0990, 0.1257, 0.1595, 0.2023, 0.2567,
    0.3257, 0.4132, 0.5243, 0.6652, 0.8440, 1.0709, 1.3588, 1.7240, 2.1874,
    2.7754, 3.5214, 4.4679, 5.6688, 7.1926, 9.1259, 11.5789, 14.6912, 18.6402,
    23.6505, 30.0076, 38.0735, 48.3075, 61.2922, 77.7672, 98.6706, 125.1926, 158.8437,
    201.5399, 255.7126, 324.4467, 411.6560, 522.3067, 662.6996, 840.8294, 1066.8393, 1353.5994,
    1717.4389, 2179.0763, 2764.7991, 3507.9607, 4450.8797, 5647.2498, 7165.1970
};

static int dering_decision(int k, int l, yuv_frame_t *rec, yuv_frame_t *org, const deblock_data_t *deblock_data, int block_size, uint8_t qp, double lambda, void *stream) {
  int skip_stride = 8;
  unsigned char bskip_y[8*8];
  for (int m=0;m<MAX_BLOCK_SIZE/block_size;m++){
    for (int n=0;n<MAX_BLOCK_SIZE/block_size;n++){
      int xpos = l*MAX_BLOCK_SIZE + n*block_size;
      int ypos = k*MAX_BLOCK_SIZE + m*block_size;
      int index = (ypos/MIN_PB_SIZE)*(rec->width/MIN_PB_SIZE) + (xpos/MIN_PB_SIZE);
      bskip_y[m*skip_stride + n] =
          deblock_data[index].mode == MODE_BIPRED || !deblock_data[index].cbp.y;
    }
  }
  int num_sb_hor = rec->width/MAX_BLOCK_SIZE;
  int num_sb_ver = rec->height/MAX_BLOCK_SIZE;
  int xpos = l*MAX_BLOCK_SIZE;
  int ypos = k*MAX_BLOCK_SIZE;
  uint8_t buf[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int ln = 6;
  int xdec = 0;
  int pli = 0;
  int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS];

  od_dering(buf, OD_BSIZE_MAX, &rec->y[ypos*rec->stride_y + xpos], rec->stride_y,
            ln, l, k, num_sb_hor, num_sb_ver, qp, xdec,
            dir, pli, bskip_y, skip_stride);


  int n = 1 << ln;
  double unfiltered_error = 0;
  double filtered_error = 0;
  double filtered_rate = 1;
  double unfiltered_rate = 1;

  /* Optimize deringing for PSNR. */
  for (int y = 0; y < n; y++) {
    for (int x = 0; x < n; x++) {
      int r;
      od_coeff p;
      od_coeff o;
      r = org->y[(ypos + y)*org->stride_y + (xpos + x)];
      p = buf[y*OD_BSIZE_MAX + x];
      filtered_error += (r - p)*(double)(r - p);
      o = rec->y[(ypos + y)*rec->stride_y + (xpos + x)];
      unfiltered_error += (r - o)*(double)(r - o);
    }
  }
  int filtered = (filtered_error + (int32_t)(lambda*filtered_rate + 0.5)) <
                 (unfiltered_error + (int32_t)(lambda*unfiltered_rate + 0.5));

  putbits(1, filtered, (stream_t*)stream);
  return filtered;
}

void encode_frame(encoder_info_t *encoder_info)
{
  int k,l;
  int width = encoder_info->width;
  int height = encoder_info->height;  
  int num_sb_hor = (width + MAX_BLOCK_SIZE - 1)/MAX_BLOCK_SIZE;
  int num_sb_ver = (height + MAX_BLOCK_SIZE - 1)/MAX_BLOCK_SIZE;
  stream_t *stream = encoder_info->stream;

  memset(encoder_info->deblock_data, 0, ((height/MIN_PB_SIZE) * (width/MIN_PB_SIZE) * sizeof(deblock_data_t)) );

  frame_info_t *frame_info = &(encoder_info->frame_info);
  double lambda_coeff;
  if (frame_info->frame_type == I_FRAME)
    lambda_coeff = encoder_info->params->lambda_coeffI;
  else if(frame_info->frame_type == P_FRAME)
    lambda_coeff = encoder_info->params->lambda_coeffP;
  else{
    if (frame_info->b_level==0)
      lambda_coeff = encoder_info->params->lambda_coeffB0;
    else if(frame_info->b_level == 1)
      lambda_coeff = encoder_info->params->lambda_coeffB1;
    else if (frame_info->b_level == 2)
      lambda_coeff = encoder_info->params->lambda_coeffB2;
    else if (frame_info->b_level == 3)
      lambda_coeff = encoder_info->params->lambda_coeffB3;
    else
      lambda_coeff = encoder_info->params->lambda_coeffB;
  }
  frame_info->lambda = lambda_coeff*squared_lambda_QP[frame_info->qp];

  putbits(1,encoder_info->frame_info.frame_type!=I_FRAME,stream);
  uint8_t qp = encoder_info->frame_info.qp;
  putbits(8,(int)qp,stream);
  putbits(4,(int)encoder_info->frame_info.num_intra_modes,stream);

  // Signal actual number of reference frames
  if (frame_info->frame_type!=I_FRAME)
    putbits(2,encoder_info->frame_info.num_ref-1,stream);

  int r;
  for (r=0;r<encoder_info->frame_info.num_ref;r++){
    putbits(6,encoder_info->frame_info.ref_array[r]+1,stream);
  }

  for (k=0;k<num_sb_ver;k++){
    for (l=0;l<num_sb_hor;l++){
      int xposY = l*MAX_BLOCK_SIZE;
      int yposY = k*MAX_BLOCK_SIZE;

      for (int ref_idx = 0; ref_idx <= frame_info->num_ref - 1; ref_idx++){
        frame_info->mvcand_num[ref_idx] = 0;
        frame_info->mvcand_mask[ref_idx] = 0;
      }
      frame_info->best_ref = -1;

      int max_delta_qp = encoder_info->params->max_delta_qp;
      if (max_delta_qp){
        /* RDO-based search for best QP value */
        int cost,min_cost,best_qp,qp0,max_delta_qp,min_qp,max_qp;
        max_delta_qp = encoder_info->params->max_delta_qp;
        min_cost = 1<<30;
        stream_pos_t stream_pos_ref;
        read_stream_pos(&stream_pos_ref,stream);
        best_qp = qp;
        min_qp = qp-max_delta_qp;
        max_qp = qp;
        for (qp0=min_qp;qp0<=max_qp;qp0++){
          cost = process_block(encoder_info,MAX_BLOCK_SIZE,yposY,xposY,qp0);
          if (cost < min_cost){
            min_cost = cost;
            best_qp = qp0;
          }
        }
        write_stream_pos(stream,&stream_pos_ref);
        process_block(encoder_info,MAX_BLOCK_SIZE,yposY,xposY,best_qp);
      }
      else{
        process_block(encoder_info,MAX_BLOCK_SIZE,yposY,xposY,qp);
      }
    }
  }

  int qpc = chroma_qp[qp];

  if (encoder_info->params->deblocking){
    deblock_frame_y(encoder_info->rec, encoder_info->deblock_data, width, height, qp);
    deblock_frame_uv(encoder_info->rec, encoder_info->deblock_data, width, height, qpc);
  }

  int sb_signal = 1;
  double lambda = frame_info->lambda;

  if (encoder_info->params->clpf){
    putbits(1, 1, stream);
    putbits(1, !sb_signal, stream);
    clpf_frame(encoder_info->rec, encoder_info->orig, encoder_info->deblock_data, stream,
               dering_decision, qp, qpc, lambda);
  }

  /* Sliding window operation for reference frame buffer by circular buffer */

  /* Store pointer to reference frame that is shifted out of reference buffer */
  yuv_frame_t *tmp = encoder_info->ref[MAX_REF_FRAMES-1];

  /* Update remaining pointers to implement sliding window reference buffer operation */
  memmove(encoder_info->ref+1, encoder_info->ref, sizeof(yuv_frame_t*)*(MAX_REF_FRAMES-1));

  /* Set ref[0] to the memory slot where the new current reconstructed frame wil replace reference frame being shifted out */
  encoder_info->ref[0] = tmp;

  /* Pad the reconstructed frame and write into ref[0] */
  create_reference_frame(encoder_info->ref[0],encoder_info->rec);

#if 0
  /* To test sliding window operation */
  int offsetx = 500;
  int offsety = 200;
  int offset_rec = offsety * encoder_info->rec->stride_y +  offsetx;
  int offset_ref = offsety * encoder_info->ref[0]->stride_y +  offsetx;
  printf("rec: %3d ",encoder_info->rec->y[offset_rec]);
  printf("ref: ");
  for (r=0;r<MAX_REF_FRAMES;r++){
    printf("%3d ",encoder_info->ref[r]->y[offset_ref]);
  }
#endif

}


