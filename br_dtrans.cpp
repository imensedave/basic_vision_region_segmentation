
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

 
/* distance transform based segmentation routines.
 *
 * copyright 2006 cambridge ontology
 */

using namespace std;


/*
  
int van_grey_seg( Matrix<short> &grey, Matrix<short> &bin ){

  Matrix<short> & edg = g_bin;
  Matrix<short> & dimg = g_bin2;

  int nr=grey.rows(), nc = grey.cols();
  
  canny_edge(grey, 800, 2000, edg);
  edg.write_rlm_ascii("edg.rlm");

  dimg = edg;
  edtrans( dimg );


  dimg.write_rlm_ascii("dist.rlm");

    
  vector<PT> seeds; 
  int num_seeds = dtrans_peaks(dimg, seeds);
  //cerr << "seed points " << num_seeds << std::endl;


  // GROW REGIONS
  vector<Dreg> dregs;
  //std::cerr << "regions ";

  
  Matrix<short> rmap(nr,nc);
  
  unsigned int nv = vsregion_grow3(grey, rmap, edg, dimg, dregs, seeds);
  std::cerr << "regions-0 " << nv << std::endl;


  rmap.write_rlm_ascii("rmap0.rlm");
  bin = rmap;
  
  return(0);
  
  small_gobit_merge(rmap, dregs, 16);
  // rmap.write_rlm_ascii("rmap1.rlm");
 
  
  int numletx = 0, last_poxel=0, maxid = 0;
  dark_murk(grey, rmap, dregs, g_bin2, bin, g_lets, numletx, last_poxel, maxid);


  bin = rmap; 
  
  //  small_merge(rmap, dregs, 150 );
  // small_merge(rmap, dregs, 50 );
  
  
  //fast_grey_edge_incorporate(grey, rmap, dregs);
  rmap.write_rlm_ascii("rmap1.rlm");
  
  
}
*/



// two pass distance transform
// edges >0, distance as -ve integers.
void  edtrans(Matrix<short> &im)
{
  int nr = im.rows(), nr2 = nr-2, nc = im.cols(), nc2 = nc-2;
  int v, mx, m=-nr-nc, nc1 = nc+1, nx1 = nc-1;
  
  for (int r = 1; r < nr - 1; r++) {
    short *p = im[r] + 1;
    
    for (short *pend = p + nc2; p < pend; p++)
      if (*p == 0){
         mx = m;
  
         v = *(p - nc1);
         if (v > 0){ *p= -1;  goto jump1;}
         else if ((v < 0) && (v > mx))  mx = v;
                  
         v = *(p - nc);
         if (v > 0){    *p= -1;  goto jump1;}
         else if ((v < 0) && (v > mx))    mx = v;
  
         v = *(p - nx1);
         if (v > 0){     *p= -1;  goto jump1;}
         else if ((v < 0) && (v > mx))    mx = v;

         v = *(p - 1);
         if (v > 0) {   *p= -1;  goto jump1;}
         v = *(p+1);
         if (v > 0) {   *p= -1;  goto jump1;}
         v = *(p + nx1);
         if (v > 0){    *p= -1;  goto jump1;}
         v = *(p + nc);
         if (v > 0){    *p= -1;  goto jump1;}
         v = *(p + nc1);
         if (v > 0){     *p= -1;  goto jump1;}

         *p = mx-1;
         jump1:;
      }
  }
  for ( int j = nc2; j >= 1; --j){
    short *p = im[nr2] + j;
    
    for (int i = nr2; i >= 1; --i, p -= nc)
      if (*p <= 0){
         mx = m;
         
         v = *(p + nc1);
         if (v > 0){     *p= -1;  goto jump2;}
         else if ((v < 0) && (v > mx))    mx = v;
         
         v = *(p + nc);
         if (v > 0){    *p= -1;  goto jump2;}
         else if ((v < 0) && (v > mx))    mx = v;
                           
         v = *(p + nx1);
         if (v > 0){    *p= -1;  goto jump2;}
         else if ((v < 0) && (v > mx))    mx = v;
         
         v = *(p - 1);
         if (v > 0) {   *p= -1;  goto jump2;}
         else if (( v < 0) && (v > mx))    mx = v;
         
         v = *(p - nc);
         if (v > 0){    *p= -1;  goto jump2;}
         else if ((v < 0) && (v > mx))    mx = v;
  
         v = *(p - nx1);
         if (v > 0){     *p= -1;  goto jump2;}
         else if ((v < 0) && (v > mx))    mx = v;

         v = *(p+1);
         if (v > 0) {   *p= -1;  goto jump2;}
         else if ((v < 0) && (v > mx))    mx = v;

         v = *(p - nc1);
         if (v > 0){ *p= -1;  goto jump2;}
         else if ((v < 0) && (v > mx))  mx = v;

  
         *p = mx-1;
         jump2:;
      }
  }
  
  return;
}



/* pick out a set of peaks */


int dtrans_peaks(Matrix<short> &dim, std::vector<PT>  &peaks)
{
  int nr = dim.rows(), nc = dim.cols(), br=nr-1, bc=nc-1;
  int label = 1;

  Matrix<short> taken(nr,nc);  // profligate...

  for (PT i(1, 1); i.r < br; ++i.r) {
    short *p = dim[i.r] + 1;
    short *s = taken[i.r] + 1;
    
    for (i.c = 1; i.c < bc; ++i.c, ++p, ++s)
      if ((*p <= -2) && (*s == 0))
	if( (*(p - 1) >= *p) && (*(p + 1) >= *p) &&
	    (*(p - nc) >= *p) && (*(p + nc) >= *p) )
	  if( ( (*(p - 1-nc) >= *p) && (*(p + 1+nc) >= *p) &&
		(*(p - nc+1) >= *p) && (*(p + nc-1) >= *p) )){
	    if (peak_test(dim, i, label, *p, taken))  {
	      i.v = *p;
	      peaks.push_back(i);
	    }
	    ++label;
	  }
  }
  sort(peaks.begin(), peaks.end());
  
  return peaks.size();
}


int  peak_test(Matrix<short> &dim, PT &pk,int label,
	       int minimum, Matrix<short> &taken)
{
   int r,c;
   int result = 1;
   std::queue<PT> points;
   int width = dim.cols();
   int height = dim.rows();
   
   for (points.push(pk); points.size() != 0; points.pop())    {
      PT p = points.front();
      
      r = p.r;
      c = p.c+1;
      if( (c>=0) && (r<width) ){
         if( dim[r][c] < minimum ){
            result = 0;
            return( result );
         }
         else if( dim[r][c] == minimum ){
            if( taken[r][c] == 0 ){
               taken[r][c] = label;
               PT tmp(r,c);
               points.push( tmp );
            }
         }
      }
      
      c = p.c-1;
      if( (c>=0) && (r<width) ){
         if( dim[r][c] < minimum ){
            result = 0;
            return( result );
         }
         else if( dim[r][c] == minimum ){
            if( taken[r][c] == 0 ){
               taken[r][c] = label;
               PT tmp(r,c);
               points.push( tmp );
            }
         }
      }
      
      
      r = p.r-1;
      if( (r>=0) && (r<height)){
         c = p.c;
         if( (c>=0) && (r<width) ){
            if( dim[r][c] < minimum ){
               result = 0;
               return( result );
            }
            else if( dim[r][c] == minimum ){
               if( taken[r][c] == 0 ){
                  taken[r][c] = label;
                  PT tmp(r,c);
                  points.push( tmp );
               }
            }
         }
         
         c = p.c+1;
         if( (c>=0) && (r<width) ){
            if( dim[r][c] < minimum ){
               result = 0;
               return( result );
            }
            else if( dim[r][c] == minimum ){
               if( taken[r][c] == 0 ){
                  taken[r][c] = label;
                  PT tmp(r,c);
                  points.push( tmp );
               }
            }
         }
         
         c = p.c-1;
         if( (c>=0) && (r<width) ){
            if( dim[r][c] < minimum ){
               result = 0;
               return( result );
            }
            else if( dim[r][c] == minimum ){
               if( taken[r][c] == 0 ){
                  taken[r][c] = label;
                  PT tmp(r,c);
                  points.push( tmp );
               }
            }
         }
      }
      
      
      r = p.r+1;
      if( (r>=0) && (r<height)){
         c = p.c;
         if( (c>=0) && (r<width) ){
            if( dim[r][c] < minimum ){
               result = 0;
               return( result );
            }
            else if( dim[r][c] == minimum ){
               if( taken[r][c] == 0 ){
                  taken[r][c] = label;
                  PT tmp(r,c);
                  points.push( tmp );
               }
            }
         }
         
         c = p.c+1;
         if( (c>=0) && (r<width) ){
            if( dim[r][c] < minimum ){
               result = 0;
               return( result );
            }
            else if( dim[r][c] == minimum ){
               if( taken[r][c] == 0 ){
                  taken[r][c] = label;
                  PT tmp(r,c);
                  points.push( tmp );
               }
            }
         }
         
         c = p.c-1;
         if( (c>=0) && (r<width) ){
            if( dim[r][c] < minimum ){
               result = 0;
               return( result );
            }
            else if( dim[r][c] == minimum ){
               if( taken[r][c] == 0 ){
                  taken[r][c] = label;
                  PT tmp(r,c);
                  points.push( tmp );
               }
            }
         }
      }
   }
   return( result );
}

void fast_grey_edge_incorporate( Matrix<short> &im, Matrix<short> &lab,
				 vector<Dreg> &vregs)
{

  int I[4] = {-1, 0, 0, 1};
  int J[4] = { 0, 1,-1, 0};
  float d1, d2, d3, rr, gg, bb, tmp2, tmp, min;
  int idx, ik, jk, mink;

  
  int nc = im.cols(); int bc = nc-1;
  int nr = im.rows(); int br = nr-1;

  //g_bin3.clear();
  //Matrix<short> &scr = g_bin3;

  Matrix<short> scr(nr,nc);

  short *pl;

  for( int i=1; i<br; i++){
    pl = lab[i]+1;
    for( int j=1; j<bc; j++, pl++){
      
      if( (*pl) == 0 ){

        min = 100000.0;
        //rr = (float) im.red()[i][j];
        gg = (float) im[i][j];
	//        bb = (float) im.blue()[i][j];

        for( int k=0; k<4; k++){
          ik = i+I[k];
          jk = j+J[k];
          if( ( lab[ik ][jk ] != 0 ) && (scr[ik][jk] == 0)){ 
            
	    //  d1 = im.red()[ik ][jk] - rr;
            d2 = im[ik ][jk] - gg;
            //d3 = im.blue()[ik ][jk] - bb;
	    tmp = fabs(d2);// + fabs(d3);

            if( tmp < min ) {
              min = tmp;
              idx = lab[ik ][jk];
              mink = k;
            }
          }
        }
        if( min < 10000.0 ){
          lab[i][j] = idx;
          scr[i][j] = 1;

          // update region colour
          Dreg &vreg = vregs[idx];   // just a reference to the region in question
          ik = i+I[mink];
          jk = j+J[mink];
          tmp = 1.0/(vreg.area +1);
          tmp2 = tmp*vreg.area;
          vreg.area++;
          //vreg.r = vreg.r*tmp2 + im.red()[ik ][jk]*tmp;
          vreg.g = vreg.g*tmp2 + im[ik ][jk]*tmp;
          //vreg.b = vreg.b*tmp2 + im.blue()[ik ][jk]*tmp;

        }
      }
    }
  }
}

/* use distance transform as a constraining factor in regions growing
 * do not allow wide open regions to propagate up hill 
 * in a distance sense
 *
 * grow smooth regions first
 *
 */

int vsregion_grow3(Matrix<short> &im, Matrix<short> &lim,
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
	
  int nr = im.rows(); int br = nr-1;
  int nc = im.cols(); int bc = nc-1;
	
  	
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
	//vrug.r = red[r][c];
	vrug.g = im[r][c];
	//vrug.b = blue[r][c];
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
							
	      //R = red[rr][cc]; G = green[rr][cc]; B = blue[rr][cc];
	      // test local colour difference
	      //dr = red[r][c] - R; dg = green[r][c] - G; db = blue[r][c] - B;
	      //drs = vreg.r - R;	 dgs = vreg.g - G;     dbs = vreg.b - B;
	      
	      G = im[rr][cc];
	      dg = im[r][c] - G;
	      dgs = vreg.g - G;
	      
	      
	      //if( (fabs(dr) < 3) && (fabs(dg) < 3) &&(fabs(db) < 3) ){ 
	      if( fabs(dg) < 3 && abs(dgs < 30 ) ){ 
		assim ++;
		//add pixel into region
		tmp = 1.0/( 1.0 + vreg.area);
		tmp2 = vreg.area*tmp;
		
		//vreg.r = vreg.r * tmp2 + R*tmp;
		vreg.g = vreg.g * tmp2 + G*tmp;
		//vreg.b = vreg.b * tmp2 + B*tmp;
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
							
	      //R = red[rr][cc]; G = green[rr][cc]; B = blue[rr][cc];
	      // test local colour difference
	      //dr = red[r][c] - R; dg = green[r][c] - G; db = blue[r][c] - B;
	      //drs = vreg.r - R;	 dgs = vreg.g - G;     dbs = vreg.b - B;

	      G = im[rr][cc];
	      dg = im[r][c] - G;	
	      dgs = vreg.g - G;
	      
	      //if( (fabs(dr) < 8) && (fabs(dg) < 8) &&(fabs(db) < 10)	&&
	      //  ( fabs(dr-dg)< 12 )&&(fabs(db -dg)< 12 )&&(fabs(dr-db)< 12 )){
		// look at central colour difference too

		if( fabs(dg) < 8){
		// look at central colour difference too
								
		  if( fabs(dgs)<30 ){
		  // add in new pixel
		  tmp = 1.0/( 1.0 + vreg.area);
		  tmp2 = vreg.area*tmp;
		  
		  //vreg.r = vreg.r * tmp2 + R*tmp;
		  vreg.g = vreg.g * tmp2 + G*tmp;
		  //vreg.b = vreg.b * tmp2 + B*tmp;
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
						
      //  R = red[rr][cc]; G = green[rr][cc]; B = blue[rr][cc];
      //drs = vreg.r - R;	 dgs = vreg.g - G;     dbs = vreg.b - B;
      //dr = red[r][c] - R; dg = green[r][c] - G; db = blue[r][c] - B;

      G = im[rr][cc];
      dgs = vreg.g - G;
       dg = im[r][c] - G;
       
							
       //if( ( (fabs(drs)<50) && (fabs(dgs)<50) &&(fabs(dbs)<50) ) || ( (fabs(dr) < 14) && (fabs(dg) < 14) &&(fabs(db) < 14) )){
      if( fabs(dgs)<30  && fabs(dg) < 14 ){
	tmp = 1.0/( 1.0 + vreg.area);
	tmp2 = vreg.area*tmp;

	//vreg.r = vreg.r * tmp2 + red[rr][cc]*tmp;
	vreg.g = vreg.g * tmp2 + im[rr][cc]*tmp;
	//vreg.b = vreg.b * tmp2 + blue[rr][cc]*tmp;
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




void find_black_regions(CxImage &im, Matrix<short> &lim, Matrix<short> &edge, std::vector<Dreg> &vregs, Matrix<int> scratch)
{
  int num, ii, jj;


  int nc = im.cols(); int bc = nc-1;
  int nr = im.rows(); int br = nr-1;

  int idx = vregs.size(); 

  int *e1=scratch[0], *e2= scratch[nr/2];


  int I[4] = {-1, 0, 0, 1};
  int J[4] = { 0, 1,-1, 0};

  short *pl, *pe;
  int ti, tj;

  for( int r=0; r<nr; r++){
    pl = lim[r];
    pe = edge[r];
    for( int c=0; c<nc; c++, pl++, pe++){

      // start to grow a region from seed point (i,j).
      if(( *pl == 0 ) && ( *pe == 0 ) ){
	num = 1;
	e1[0] = r;
	e2[0] = c;

	*pl = idx;
	
	int i=0;
	while( i<num ){
	  for( int j=0; j<4; j++){ // look in 4-nhbd of edge for next pixels
	    ii = e1[i] + I[j];
	    if(( ii >= 0 ) && ( ii < nr)){
	      jj = e2[i] + J[j];
	      
	      if(( jj >= 0 ) && ( jj < nc)){
		
		if(  lim[ii][jj] == 0  &&  edge[ii][jj] == 0 ){      // not trough pixel
		  lim[ii][jj] = idx;
		  edge[ii][jj] == 8;
		  
		  e1[num] = ii;
		  e2[num] = jj;
		  num++;
		}
	      }
	    }
	  }
	  i++;
	}
	if( num >= 10 ){  // make new region if area > 20 pixels
	  Dreg vreg;
	  vreg.label = idx;
	  vreg.area = num;
	  for( int i=0; i<num; i++){
	    ti = e1[i]; tj = e2[i];
	    lim[ ti ][ tj ] = idx;
	    vreg.r += im.red()[ti][tj];
	    vreg.g += im.green()[ti][tj];
	    vreg.b += im.blue()[ti][tj];
	    
	  }
	  vreg.r /= num;
	  vreg.g /= num;
	  vreg.b /= num;
	  vregs.push_back(vreg);
	  idx++;
	}
	else{
	  for( int i=0; i<num; i++){
	    lim[ e1[i] ][ e2[i] ] = idx;
	  }
	}
      }
    }
  }

  //  std::cerr << "now have " << vregs.size() << " regions" << std::endl;
  return;
  
} 
