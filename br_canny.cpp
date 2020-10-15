
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <vector>
#include <set>
#include <string>
#include <cfloat>
#include "Imgx.h"
#include "IM_struct.h"
#include "IM_func.h"



// copyright David Sinclair 2020


/* released under the Apache2 license.

   https://www.apache.org/licenses/LICENSE-2.0


   dr.sinclair@@gmail.com
*/



#define debug_flag 1


  
int  canny(Matrix<short> &im, float low, float high, Matrix<short> &edge)
{
  int nr = im.rows(); int br = nr-2;
  int nc = im.cols(); int bc = nc-2;
  int np1 = nc+1, nm1 = nc-1;


  // scratch space for intermediate results.
  
  Matrix<int> magx(nr,nc);
  Matrix<short> drx(nr,nc);
  Matrix<short> dcx(nr,nc);

  // traditionally canny has a gaussian smoothing step
  // but modern cameras have low enough noise for this to be ignored.
  
  
  // compute derivatives
  for(int r=2; r<br; r++){
    short *pr = im[r], *pend;
    pend = pr + bc;
    pr+=2;
		
    for(int c=1; pr<pend; pr++, c++ ){
      int dc = (im[r][c+1] -  im[r][c-1])*2
	+ im[r-1][c+1] - im[r-1][c-1]
	+ im[r+1][c+1] -  im[r+1][c-1];
      dcx[r][c] = dc;

      int dr =  (im[r+1][c] -  im[r-1][c])*2 +
	im[r+1][c-1] - im[r-1][c-1]+ 
	im[r+1][c+1] -  im[r-1][c+1];
	
      drx[r][c] = dr;
	
      int mag = dr*dr + dc*dc;
      if( mag > low ){
	magx[r][c] = mag;
      }

    }
  }


  if( debug_flag > 0 ){
    // magx.write_rlm_ascii("magx.rlm");
  }
	
  edge.clear();

	
  // non-maximum supression
  int count;
  int val;
  int nms_strength=0, vr, vc;
  std::vector<int> vi, vj, vals;
  vi.reserve((nr*nc)/8);
  vj.reserve((nr*nc)/8);
  vals.reserve((nr*nc)/8);
	
  count = 0;
  for(int i=1; i<br; i++){
    int *ptf = magx[i];
    int *pend  = ptf + bc;
    ptf  += 1;
		
    for(int j=1 ; ptf<pend; ptf++, j++){
      if( *ptf > low ){
	vr =  drx[i][j];
	vc =  dcx[i][j];
	val = *ptf;
				
	if( vc>0 & vr >= 0){   // 1st quadrant
	  if( vr/vc < 0.414 ) {
	    if( magx[i][j-1]<val & val >=magx[i][j+1] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	  else if( vr/vc > 2.415 ){
	    if( magx[i-1][j]<val & val >=magx[i+1][j] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	  else{
	    if( magx[i-1][j-1]<val & val >=magx[i+1][j+1] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	}
	else if( vr<0 & vc>=0){ // 4th quadrant
	  if( -vc/vr > 2.415 ){
	    if( magx[i][j-1]<val & val >=magx[i][j+1] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	  else if( -vc/vr < 0.414 ){
	    if( magx[i-1][j]<val & val >=magx[i+1][j] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	  else{
	    if( magx[i-1][j+1]<val & val >=magx[i+1][j-1] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	}
	else if( vr>0 & vc<=0){ // 2nd quadrant
	  if( -vc/vr > 2.415 ){
	    if( magx[i][j-1]<val & val >=magx[i][j+1] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	  else if( -vc/vr < 0.414 ){
	    if( magx[i-1][j]<val & val >=magx[i+1][j] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	  else{
	    if( magx[i-1][j+1]<val & val >=magx[i+1][j-1] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	}
	else if(vc != 0){ // 3rd quadrant
	  if( vr/vc < 0.414 ){
	    if( magx[i][j-1]<val & val >=magx[i][j+1] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	  else if( vr/vc > 2.415 ){
	    if( magx[i-1][j]<val & val >=magx[i+1][j] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	  else{
	    if( magx[i-1][j-1]<val & val >=magx[i+1][j+1] ){
	      vi.push_back( i);
	      vj.push_back( j);
	      vals.push_back( val);
	      count++;
	    }
	  }
	}
      }
    }
  }

  magx.clear();
  Matrix<int> &smot = magx;

  for(int i=0; i<count; i++){
    smot[vi[i]][vj[i]] = vals[i];
  }
			
  int num, tr,tc;
 
  int *ii, *jj;
  ii = new int [count];
  jj = new int [count];
		
  for(int i=0; i<count; i++){
		
    if( smot[vi[i] ][vj[i] ] > high){ // begin edge follow 8-connected!
      num = 1;
      ii[0] = vi[i]; jj[0] = vj[i];
      for( int k=0; k<num; k++){  // four connected ordering first!
	tr = ii[k]-1;
	tc = jj[k];
	if(smot[tr][tc]>0){
	  ii[num] = tr; jj[num] = tc;
	  smot[ tr ][ tc ] = 0;
	  num ++;
	}
	tr = ii[k]+1;
	if(smot[tr][tc]>0){
	  ii[num] = tr; jj[num] = tc;
	  smot[ tr ][ tc ] = 0;
	  num ++;
	}
	tr = ii[k]; tc = jj[k]-1;
	if(smot[tr][tc]>0){
	  ii[num] = tr; jj[num] = tc;
	  smot[ tr ][ tc ] = 0;
	  num ++;
	}
	tc = jj[k]+1;
	if(smot[tr][tc]>0){
	  ii[num] = tr; jj[num] = tc;
	  smot[ tr ][ tc ] = 0;
	  num ++;
	}
				
	// 8-connected bits
				
	tr = ii[k]-1;
	if(smot[tr][tc]>0){
	  ii[num] = tr; jj[num] = tc;
	  smot[ tr ][ tc ] = 0;
	  num ++;
	}
	tr = ii[k]+1;
	if(smot[tr][tc]>0){
	  ii[num] = tr; jj[num] = tc;
	  smot[ tr ][ tc ] = 0;
	  num ++;
	}
	tc = jj[k]-1;
	if(smot[tr][tc]>0){
	  ii[num] = tr; jj[num] = tc;
	  smot[ tr ][ tc ] = 0;
	  num ++;
	}
	tr = ii[k]-1;
	if(smot[tr][tc]>0){
	  ii[num] = tr; jj[num] = tc;
	  smot[ tr ][ tc ] = 0;
	  num ++;
	}
      }
      if( num > 8 ){ // drop short edges.
				
	for(int q=0; q<num; q++){
	  edge[ii[q]][jj[q]] = 1;
	}

	  
	
      }
    }
  }

  if( debug_flag > 0 ) {
    // edge.write_rlm_ascii("edge.rlm");
  }	

  delete [] ii;
  delete [] jj;
	
  return(count);
}


