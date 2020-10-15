
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include <vector>
#include <set>
#include <string>
#include <cfloat>
#include <stack>
#include <queue>


#include "Imgx.h"
#include "IM_struct.h"
#include "IM_func.h"



// copyright David Sinclair 2020


/* released under the Apache2 license.

   https://www.apache.org/licenses/LICENSE-2.0


   dr.sinclair@@gmail.com
*/



#define debug_flag 1


/* segment colour images into coloured regions.

*/

// use XYZ colour space for colour metirc.
// allow more brightness than Y,Z  variation
// 

void colour_segmentation_Canny_XYZ( CxImage &xyz, Matrix<short> &rmap,  Matrix<int> &scratch )
{

  int nr = xyz.rows(), nc = xyz.cols();
  Matrix<short> edg(nr,nc);

  
  canny(xyz.red(), 800, 2000, edg);


  if( debug_flag > 0 ){
    // edg.write_rlm_ascii("edges.rlm");
    
  }

  Matrix<short> mag(edg);
    
    
  edtrans( mag );

  if( debug_flag > 0 ){
    mag.write_rlm_ascii("dist.rlm");
  }
  
  vector<PT> seeds; 
  int num_seeds = dtrans_peaks(mag, seeds);
  cerr << "seed points " << num_seeds << std::endl;


  // GROW solid colour REGIONS
  vector<Dreg> dregs;
  
  unsigned int nv = vsregion_grow3_xyz(xyz, rmap, edg, mag, dregs, seeds);
  std::cerr << "regions-0 " << nv << std::endl;



  find_black_regions(xyz, rmap, edg, dregs, scratch);

  
  if( debug_flag > 0 ){
    rmap.write_rlm_ascii("region_map.rlm");  
  }

  // you can write any region merging or edge incorporation bits you like yourself.


  
  return;
}


/* cim is assumed to be xyz colour spaced

 */


int vsregion_grow3_xyz(CxImage &cim, Matrix<short> &lim,
			  Matrix<short> &edge, Matrix<short> &vog,
			  std::vector<Dreg> &vregs, std::vector<PT> &seeds)
{
  int i,r,c, tr, tc, rr,cc;
  int cent, v0, v1, v2;
  int label, assim = 0;
  float dr, dg, db, drs, dgs, dbs, R, G, B, tmp, tmp2;
  int ii[40], jj[40], num;
	
  queue<UDL> bndry, failed, fx, blank;
  
  UDL bnd;
  
  int nr = lim.rows(); int br = nr-1;
  int nc = lim.cols(); int bc = nc-1;
	

  Matrix<short> &red = cim.red();     // X
  Matrix<short> &green = cim.green(); // Y 
  Matrix<short> &blue = cim.blue();   // Z

  
  
  cent = 1;
	
  Dreg vrug;
  vregs.push_back(vrug);
  int num_vor = seeds.size();
	
  for(i=0; i<num_vor; i++){
    if( lim[seeds[i].r ][ seeds[i].c ] == 0){
      r = seeds[i].r; c = seeds[i].c;
      v0 = vog[seeds[i].r ][ seeds[i].c ];
			
      //cerr << v0 << " " ;
			
      if( (r>0) && (r<br) && (c>0) && (c<bc) ){
				
	lim[ r ][ c ] = cent;  // seed label
	num = 0;
	ii[num] = r; jj[num] = c;
				
	vrug.area = 1;
	vrug.r = red[r][c];
	vrug.g = green[r][c];
	vrug.b = blue[r][c];
	vrug.label = cent;
				
	//				vregs.push_back(vrug);
				
	bnd.lab = cent;
	bnd.ui = r;
	bnd.uj = c;
				
	tr = r-1; tc = c;
	if( edge[tr][tc] == 0){  // good boundary point
	  bnd.di = tr;
	  bnd.dj = tc;
					
	  bndry.push(bnd);
	}
				
	tr = r+1;
	if( edge[tr][tc] == 0){  // good boundary point
	  bnd.di = tr;
					
	  bndry.push(bnd);
	}
	tr = r; tc = c-1;
	if( edge[tr][tc] == 0){  // good boundary point
	  bnd.di = tr;
	  bnd.dj = tc;
					
	  bndry.push(bnd);
	}
	tc = c+1;
	if( edge[tr][tc] == 0){  // good boundary point
	  bnd.dj = tc;
					
	  bndry.push(bnd);
	}
				
	label = cent;
	Dreg &vreg = vrug;
	// grow regions using the boundary que
	// store pixels up to a certain size?
	// and then delete small regions?
	failed = blank;
				
	while( bndry.size() != 0){
	  UDL &crnt = bndry.front();
	  rr = crnt.di;   cc = crnt.dj;
					
	  // check to see if down pixel is taken
	  if( lim[rr][cc] == 0){
	    r = crnt.ui;    c = crnt.uj;
	    //v1 = vog[r][c];
	    v2 = vog[rr][cc];
	    if( (v0 < -3) || ( v2 >= v0 )){ // do not match away from edges.
	      //if( 1){
							
	      R = red[rr][cc];
	      G = green[rr][cc];
	      B = blue[rr][cc];

	      // test local colour difference
	      dr = red[r][c] - R; dg = green[r][c] - G; db = blue[r][c] - B;
	      drs = vreg.r - R;	 dgs = vreg.g - G;     dbs = vreg.b - B;
	      
	      //G = im[rr][cc];
	      //dg = im[r][c] - G;
	      //dgs = vreg.g - G;
	      
	      
	      if( (fabs(dr) < 4) && (fabs(dg) < 3) &&(fabs(db) < 3) &&
		  fabs(drs) < 30  &&
		  fabs(dgs) < 10  &&
		  fabs(dbs) < 10 
		  ){ 
		//if( fabs(dg) < 3 && abs(dgs < 30 ) ){ 
		assim ++;
		//add pixel into region
		tmp = 1.0/( 1.0 + vreg.area);
		tmp2 = vreg.area*tmp;
		
		vreg.r = vreg.r * tmp2 + R*tmp;
		vreg.g = vreg.g * tmp2 + G*tmp;
		vreg.b = vreg.b * tmp2 + B*tmp;
		vreg.area ++;
								
		lim[rr][cc] = label;
		//if( num < 8 ){
		//num++;
		//ii[num] = rr;	jj[num] = cc;
		//}
		// look for new bndry pixels
								
		bnd.lab = label;
		bnd.ui = rr;	bnd.uj = cc;
								
		tr = rr-1; tc = cc;
		if( (tr>=0) && (tr<nr) && (tc>=0) && (tc<nc) ){
		  if( (edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){ // good boundary point
		    bnd.di = tr;
		    bnd.dj = tc;
										
		    bndry.push(bnd);
		  }
		}
		tr = rr+1;
		if( (tr>=0) && (tr<nr) && (tc>=0) && (tc<nc) ){
		  if(( edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){
		    bnd.di = tr;
		    bnd.dj = tc;
										
		    bndry.push(bnd);
		  }
		}
		tr = rr; tc = cc-1;
		if( (tr>0) && (tr<nr) && (tc>0) && (tc<nc) ){
		  if(( edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){
		    bnd.di = tr;
		    bnd.dj = tc;
										
		    bndry.push(bnd);
		  }
		}
		tc = cc+1;
		if( (tr>0) && (tr<nr) && (tc>0) && (tc<nc) ){
		  if(( edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){
		    bnd.di = tr;
		    bnd.dj = tc;
										
		    bndry.push(bnd);
		  }
		}
	      }
	      else{
		failed.push(crnt);
	      }
	    }
	    else{
	      failed.push(crnt);
	    }
	  }
	  bndry.pop();
	}
	// now grow out region from failed boundary.
				
	while( failed.size() != 0){
	  UDL &crnt = failed.front();
	  rr = crnt.di;   cc = crnt.dj;
					
	  // check to see if down pixel is taken
	  if( lim[rr][cc] == 0){
	    r = crnt.ui;    c = crnt.uj;
	    v2 = vog[rr][cc];
	    if( (v0 < -3) || ( v2 >= v0 )){ // do not match away from edges.
	      
	      R = red[rr][cc]; G = green[rr][cc]; B = blue[rr][cc];
	      // test local colour difference
	      dr = red[r][c] - R; dg = green[r][c] - G; db = blue[r][c] - B;
	      drs = vreg.r - R;	 dgs = vreg.g - G;     dbs = vreg.b - B;

	      //G = im[rr][cc];
	      //dg = im[r][c] - G;	
	      //dgs = vreg.g - G;
	      
	      if( (fabs(dr) < 10) && (fabs(dg) < 4) &&(fabs(db) < 4) &&
		  fabs(drs) < 30  &&
		  fabs(dgs) < 10  &&
		  fabs(dbs) < 10  
		  ){
		//	        ( fabs(dr-dg)< 12 )&&(fabs(db -dg)< 12 )&&(fabs(dr-db)< 12 )){
		// look at central colour difference too
		
		//if( fabs(dg) < 8){
		// look at central colour difference too
		
		if( fabs(drs)<30 ){
		  // add in new pixel
		  tmp = 1.0/( 1.0 + vreg.area);
		  tmp2 = vreg.area*tmp;
		  
		  vreg.r = vreg.r * tmp2 + R*tmp;
		  vreg.g = vreg.g * tmp2 + G*tmp;
		  vreg.b = vreg.b * tmp2 + B*tmp;
		  vreg.area ++;
		  lim[rr][cc] = label;
									
		  // add in new possible boundary pixels
		  bnd.ui = rr;	bnd.uj = cc;
		  bnd.lab = label;
									
		  tr = rr-1; tc = cc;
		  if( (tr>=0) && (tr<nr) && (tc>=0) && (tc<nc) ){
		    if( (edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){  // good boundary point
		      bnd.di = tr;
		      bnd.dj = tc;
											
		      failed.push(bnd);
		    }
		  }
		  tr = rr+1;
		  if( (tr>=0) && (tr<nr) && (tc>=0) && (tc<nc) ){
		    if(( edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){
		      bnd.di = tr;
		      bnd.dj = tc;
											
		      failed.push(bnd);
		    }
		  }
		  tr = rr; tc = cc-1;
		  if( (tr>0) && (tr<nr) && (tc>0) && (tc<nc) ){
		    if(( edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){
		      bnd.di = tr;
		      bnd.dj = tc;
											
		      failed.push(bnd);
		    }
		  }
		  tc = cc+1;
		  if( (tr>0) && (tr<nr) && (tc>0) && (tc<nc) ){
		    if(( edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){
		      bnd.di = tr;
		      bnd.dj = tc;
											
		      failed.push(bnd);
		    }
		  }
									
		}
		else{
		  fx.push(crnt);
		}
	      }
	      else{
		fx.push(crnt);
	      }
	    }
	  }
	  failed.pop();
	}
			
	//if( num < 8 ){ // after region growing: delete labeled poxels
	//for(; num >=0; num--){
	//lim[ii[num]][jj[num]] = 0;
	//}
	//}
	//else{
	vregs.push_back(vrug);
	cent ++;
      }
    }
  }
	

  //	cerr << " init " << vregs.size() << " " << endl;
	
  //	lim.write_rlm_ascii("rmap1.rlm");
  //write_vregs(vregs, "rmap1.txl");

	
  for( r=1; r<br; r++){
    for( c=1; c<bc; c++){
      if( (lim[ r ][ c ] == 0 )&&( edge[ r ][ c ] == 0)){
	bnd.di = r;
	bnd.dj = c;
	if( lim[ r+1 ][ c ] > 0 ){
	  bnd.lab = lim[ r+1 ][ c ];
	  bnd.ui = r+1;
	  bnd.uj = c;
	  failed.push(bnd);
	}
	if( lim[ r-1 ][ c ] > 0 ){
	  bnd.lab = lim[ r-1 ][ c ];
	  bnd.ui = r-1;
	  bnd.uj = c;
	  failed.push(bnd);
	}
	if( lim[ r ][ c+1 ] > 0 ){
	  bnd.lab = lim[ r ][ c+1 ];
	  bnd.ui = r;
	  bnd.uj = c+1;
	  failed.push(bnd);
	}
	if( lim[ r][ c-1 ] > 0 ){
	  bnd.lab = lim[ r ][ c-1 ];
	  bnd.ui = r;
	  bnd.uj = c-1;
	  failed.push(bnd);
	}
      }
    }
  }
	
	
	
  // morph boundaries out to edges (with colour gate!)

  while( failed.size() != 0){
    UDL &crnt = failed.front();
    rr = crnt.di;   cc = crnt.dj;
		
    // check to see if down pixel is taken
    if( lim[rr][cc] == 0){
      r = crnt.ui;    c = crnt.uj;
      label = crnt.lab;
      Dreg &vreg = vregs[label];
						
      R = red[rr][cc]; G = green[rr][cc]; B = blue[rr][cc];
      drs = vreg.r - R;	 dgs = vreg.g - G;     dbs = vreg.b - B;
      dr = red[r][c] - R; dg = green[r][c] - G; db = blue[r][c] - B;

      //G = im[rr][cc];
      //dgs = vreg.g - G;
      //dg = im[r][c] - G;
       
							
      if( ( (fabs(drs)<50) && (fabs(dgs)<10) &&(fabs(dbs)<10) ) || ( (fabs(dr) < 14) && (fabs(dg) < 8) &&(fabs(db) < 8) )){
	//if( fabs(dgs)<30  && fabs(dg) < 14 ){
	tmp = 1.0/( 1.0 + vreg.area);
	tmp2 = vreg.area*tmp;

	vreg.r = vreg.r * tmp2 + red[rr][cc]*tmp;
	vreg.g = vreg.g * tmp2 + green[rr][cc]*tmp;
	vreg.b = vreg.b * tmp2 + blue[rr][cc]*tmp;
	vreg.area ++;
	lim[rr][cc] = label;
			
	// add in new possible boundary pixels
	lim[rr][cc] = label;
	bnd.ui = rr;	bnd.uj = cc;
	bnd.lab = label;
			
	tr = rr-1; tc = cc;
	if( (tr>=0) && (tr<nr) && (tc>=0) && (tc<nc) ){
	  if( (edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){  // good boundary point
	    bnd.di = tr;
	    bnd.dj = tc;

	    failed.push(bnd);
	  }
	}
	tr = rr+1;
	if( (tr>=0) && (tr<nr) && (tc>=0) && (tc<nc) ){
	  if(( edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){
	    bnd.di = tr;
	    bnd.dj = tc;

	    failed.push(bnd);
	  }
	}
	tr = rr; tc = cc-1;
	if( (tr>0) && (tr<nr) && (tc>0) && (tc<nc) ){
	  if(( edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){
	    bnd.di = tr;
	    bnd.dj = tc;

	    failed.push(bnd);
	  }
	}
	tc = cc+1;
	if( (tr>0) && (tr<nr) && (tc>0) && (tc<nc) ){
	  if(( edge[tr][tc] == 0) && (lim[tr][tc] == 0) ){
	    bnd.di = tr;
	    bnd.dj = tc;

	    failed.push(bnd);
	  }
	}
      }
    }
    failed.pop();
  }
	
	
  return(cent-1);
}

__attribute__((visibility("hidden")))  
int br_colour(CxImage &cim, Matrix<short>  &lim, CxImage &cim2 )
{
  int nr = cim.rows(), nc = cim.cols(); 
  int br = nr-1, bc = nc-1;
  int ncx = nc-1, nc1 = nc+1;



  
  struct timeval t0, t1;
  
  std::vector<PT> ptx, ptm;
  int cnt = 0;
  
  for( int r=1; r<br; r++){
    short *p = lim[r];
    short *pend = p+nc;
    for( ; p<pend; p++){
      if( *p >  cnt){
	cnt = *p;
      }
    }
  }

  cnt +=4;
  long *R, *G, *B, *nums;
  R = new long [cnt];
  G = new long [cnt];
  B = new long [cnt];
  nums = new long [cnt];

  memset(R,0,cnt*sizeof(long));
  memset(G,0,cnt*sizeof(long));
  memset(B,0,cnt*sizeof(long));
  memset(nums,0,cnt*sizeof(long));

  for( int r=0; r<nr; r++){
    short *p = lim[r];
    short *pend = p+nc;
    for( int c =0; p<pend; p++, c++){
      //if( *p > 0){
      if( abs(*p) > 0){
	int id = abs(*p);
	R[id] += cim.r[r][c];
	G[id] += cim.g[r][c];
	B[id] += cim.b[r][c];
	nums[id] ++;
      }
    }
  }

  
  for( int id=1; id<cnt; id++){
    if(	nums[id] > 0 ){
      R[id] /= nums[id];
      G[id] /= nums[id];
      B[id] /= nums[id];
    }
  }

  for( int r=0; r<nr; r++){
    short *p = lim[r];
    short *pend = p+nc;
    for( int c =0; p<pend; p++, c++){
      //if( *p > 0){
      if( abs(*p) > 0){
	int id = abs(*p);
	cim2.r[r][c] = R[id];
	cim2.g[r][c] = G[id];
	cim2.b[r][c] = B[id];
      }
    }
  }


  delete [] R;
  delete [] G;
  delete [] B;
  delete [] nums;

  
  return(cnt);
  
}
  
/*
RGB to XYZ & XYZ to RGB

RGB values in a particular set of primaries can be transformed to and from CIE XYZ via a 3x3 matrix transform. These transforms involve tristimulus values, that is a set of three linear-light components that conform to the CIE color-matching functions. CIE XYZ is a special set of tristimulus values. In XYZ, any color is represented as a set of positive values.

To transform from XYZ to RGB (with D65 white point), the matrix transform used is [3]:

   [ R ]   [  3.240479 -1.537150 -0.498535 ]   [ X ]
   [ G ] = [ -0.969256  1.875992  0.041556 ] * [ Y ]
   [ B ]   [  0.055648 -0.204043  1.057311 ]   [ Z ].

The range for valid R, G, B values is [0,1]. Note, this matrix has negative coefficients. Some XYZ color may be transformed to RGB values that are negative or greater than one. This means that not all visible colors can be produced using the RGB system.

The inverse transformation matrix is as follows:

   [ X ]   [  0.412453  0.357580  0.180423 ]   [ R ] **
   [ Y ] = [  0.212671  0.715160  0.072169 ] * [ G ]
   [ Z ]   [  0.019334  0.119193  0.950227 ]   [ B ].

*/


// approximate conversion of RGB to XYZ

void RGB2XYZ( CxImage &cim, CxImage &xyz)
{

  int a11, a12, a13,
    a21, a22, a23,
    a31, a32, a33;
  a11 = (int) (1024*0.412453);   a12 = (int) (1024*0.357580);  a13 = (int) (1024*0.180423);
  a21 = (int) (1024*0.212671);   a22 = (int) (1024*0.715160);  a23 = (int) (1024*0.072169);
  a31 = (int) (1024*0.019334);   a32 = (int) (1024*0.119193);  a33 = (int) (1024*0.950227);

  int nr = cim.red().rows(), nc = cim.red().cols();

  Matrix<short> &R = cim.red(); 
  Matrix<short> &G = cim.green(); 
  Matrix<short> &B = cim.blue(); 

  Matrix<short> &X = xyz.red(); 
  Matrix<short> &Y = xyz.green(); 
  Matrix<short> &Z = xyz.blue(); 


  short *pR = R[0], *pG = G[0], *pB = B[0];
  short *pX = X[0], *pY = Y[0], *pZ = Z[0];

  short *pond = pR + nr*nc;
  for( ; pR<pond; pR++, pG++, pB++, pX++, pY++, pZ++){
    int v = *pR * a11 +  *pG * a12 + *pB * a13;  
    *pX = v>>10;

    v = *pR * a21 +  *pG * a22 + *pB * a23;  
    *pY = v>>10;

    v = *pR * a31 +  *pG * a32 + *pB * a33;  
    *pZ = v>>10;


  }
  
}
