#ifndef _IM_func_h
#define _IM_func_h


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <set>
#include <vector>


#include "Imgx.h"
#include "IM_struct.h"

/* copyright David SInclair 2020.

   released under the Apache2.0 license:
   https://www.apache.org/licenses/LICENSE-2.0



 */



int read_pgm( Matrix<short> &im, char *fname);

void C_to_BW( CxImage &cim, Matrix<short> &out );

int  canny(Matrix<short> &im, float low, float high, Matrix<short> &edge);

// colour segmentation bits.
  
void colour_segmentation( CxImage &cim, Matrix<short> &grey, Matrix<short> &bx,  Matrix<int> &scratch );

void Redge_model( Matrix<short> &im, Matrix<short> &bx,  Matrix<int> &scratch );

void Redge_model_B(Matrix<short> &im, Matrix<short> &bx, Matrix<short> &by, Matrix<int> &scratch );

void Redge_cleanup( Matrix<short> &lim,  Matrix<short> &fog, Matrix<int> &scratch, int hyst);

void edge_morph( Matrix<short> &im, Matrix<short> &lim,  Matrix<short> &edg);

void  edtrans(Matrix<short> &im);

int dtrans_peaks(Matrix<short> &dim, std::vector<PT>  &peaks);

int  peak_test(Matrix<short> &dim, PT &pk,int label, int minimum, Matrix<short> &taken);

int vsregion_grow3(Matrix<short> &im, Matrix<short> &lim,  Matrix<short> &edge, Matrix<short> &vog,
		   std::vector<Dreg> &vregs, std::vector<PT> &seeds);

int vsregion_grow3_colour(CxImage &cim, Matrix<short> &lim, Matrix<short> &edge, Matrix<short> &vog,
			  std::vector<Dreg> &vregs, std::vector<PT> &seeds);

int vsregion_grow3_xyz(CxImage &cim, Matrix<short> &lim, Matrix<short> &edge, Matrix<short> &vog,
		       std::vector<Dreg> &vregs, std::vector<PT> &seeds);

int br_colour(CxImage &cim, Matrix<short>  &lim, CxImage &cim2 );

void RGB2XYZ( CxImage &cim, CxImage &xyz);

void colour_seg_XYZ( CxImage &cim, Matrix<short> &grey, Matrix<short> &bx,  Matrix<int> &scratch );

void colour_segmentation_Canny_XYZ( CxImage &XYZ, Matrix<short> &rmap,  Matrix<int> &scratch );

void find_black_regions(CxImage &im, Matrix<short> &lim, Matrix<short> &edge, std::vector<Dreg> &vregs, Matrix<int> scratch);



#endif


