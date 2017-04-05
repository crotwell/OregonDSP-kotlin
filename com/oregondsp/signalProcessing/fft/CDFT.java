// Copyright (c) 2011  Deschutes Signal Processing LLC
// Author:  David B. Harris

//  This file is part of OregonDSP.
//
//    OregonDSP is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Lesser General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    OregonDSP is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with OregonDSP.  If not, see <http://www.gnu.org/licenses/>.


package com.oregondsp.signalProcessing.fft;


/**
 *  Class to calculate the complex discrete Fourier transform of a complex sequence and its inverse using the split-radix algorithm.
 *  
 *  <p>This class is designed for efficient calculation of many discrete Fourier transforms of the
 *  same length.  It is limited to transform lengths that are powers of two and greater than or 
 *  equal to 32.  The class recursively constructs and links smaller DFTs with hard-wired array indices 
 *  to minimize index calculations during the overall DFT evaluation.  This approach may produce large run-time 
 *  images for very large DFTs (> 32768).  Special hand-coded implementations of length 8 and 16 DFTs eliminate 
 *  many unnecessary calculations.  The code uses precomputed sine and cosine tables and does not implement 
 *  in-place calculations in order to eliminate the bit reversal step.  Consequently, this implementation 
 *  trades memory for speed.</p>
 *  
 *  <p> Example of use:</p>
 *  <p>
 *  <font face="courier">
 *   int N &nbsp&nbsp&nbsp&nbsp&nbsp= 16384;<BR>
 *   int log2N &nbsp= 14;<BR>
 *   float[] xr = new float[N];<BR>
 *   float[] xi = new float[N];<BR>
 *   float[] Xr = new float[N];<BR>
 *   float[] Xi = new float[N];<BR>
 *   CDFT Xfm = new CDFT( log2N );<BR>
 *   <BR>
 *   // load data<BR>
 *   for ( int i = 0;  i < N;  i++ ) {<BR>
 *     &nbsp xr[i] = ...<BR>
 *     &nbsp xi[i] = ...<BR>
 *   }<BR>
 *   <BR>
 *   // evaluate transform of data<BR>
 *   Xfm.evaluate( xr, xi, Xr, Xi );<BR>
 *  </font>
 *  </p>
 *  
 *  <p>The real and imaginary parts of the transform are stored in Xr and Xi in natural order, with the zeroth
 *  discrete frequency value in Xr(0) and Xi(0), and the N-1st value ( 2*pi*(N-1)/N ) in Xr(N-1) and Xi(N-1).
 *  </p>
 *  
 *  <p>As long as the transform size does not change, the CDFT object does not need to be reinstantiated.  
 *  Consequently, the data arrays can be reloaded and the evaluate method invoked to compute additional 
 *  DFTs without incurring the cost of CDFT object instantiation.</p>
 *  
 *  <p>It may happen in some applications that the array arguments in the evaluate() and evaluateInverse() 
 *  methods never change, i.e. the same arrays are used repeatedly.  Since this implementation is recursive,
 *  the input and output arrays are recursively linked down the chain of smaller DFTs that implement the full
 *  DFT.  This linking operation can be avoided when the arguments to evaluate() and evaluateInverse() never vary.
 *  For this circumstance an alternative constructor is provided, that links the input and output arrays at 
 *  construction time (for a slight performance improvement).  To avoid relinking arrays, this constructor should 
 *  be paired with the evaluate() and evaluateInverse() methods that have NO arguments.  Example:
 *  
 *  <p>
 *  <font face="courier">
 *   CDFT Xfm = new CDFT( xr, xi, Xr, Xi, log2N );<BR>
 *   <BR>
 *   // load data<BR>
 *   for ( int i = 0;  i < N;  i++ ) {<BR>
 *    &nbsp xr[i] = ...<BR>
 *    &nbsp xi[i] = ...<BR>
 *   }<BR>
 *   <BR>
 *   // evaluate transform of data<BR>
 *   Xfm.evaluate();<BR>
 *  </font>
 *  </p>  
 *   
 *  <p>For the inverse transform in this usage, the roles of (xr,xi) and (Xr,Xi) are reversed.  The pair 
 *  (xr,xi) contains the transform real and imaginary parts in natural order, and upon execution of 
 *  evaluateInverse(), the pair (Xr,Xi) contains the real and imaginary parts of the corresponding sequence 
 *  (inverse transform).
 *  </p>
 *  
 *  <p>See "On Computing the Split-Radix FFT", Sorensen, H. V., Heideman, M. T. and Burrus, C. S.
 *  IEEE TRANSACTIONS ON ACOUSTICS, SPEECH, AND SIGNAL PROCESSING, VOL. ASSP-34, NO. 1, 
 *  FEBRUARY, 1986, pp. 152-156.</p>
 *  
 *  @author David B. Harris,   Deschutes Signal Processing LLC
 */
public class CDFT {
  
  private float[] yr;
  private float[] yi;
  private boolean arraysUnlinked;

  private float[] c;
  private float[] c3;
  private float[] s;
  private float[] s3;

  int N;
  int log2N;

  private CDFTsr dft;



  /**
   * Default constructor.
   */
  public CDFT() {
  }
  
  
 
  /** 
   * Constructs a CDFT instance without references to sequence and transform arrays
   * @param log2N       base-2 logarithm of the length of the transform
   */
  public CDFT( int log2N ) {
    
    if ( log2N < 3 ) throw new IllegalArgumentException( "DFT size must be >= 8" );
    arraysUnlinked = true;
    
    this.log2N = log2N;
    N = 1 << log2N;

    createTable();

    if (      log2N == 3 )
      dft = new CDFTsr8(  0, 1, 0 );
    else if ( log2N == 4 )
      dft = new CDFTsr16( 0, 1, 0 );
    else if ( log2N >= 5 ) {
      dft = new CDFTsr( log2N, c, c3, s, s3 );
    } 
    
  }
  
  
  
  /** 
   * evaluates the DFT with specified sequence and transform arrays
   * @param xr          float array containing sequence real part
   * @param xi          float array containing sequence imaginary part
   * @param Xr          float array containing transform real part
   * @param Xi          float array containing transform imaginary part
   */
  public void evaluate( float[] xr, float[] xi, float[] Xr, float[] Xi ) {
    this.yr = Xr;
    this.yi = Xi;
    dft.link( xr, xi, Xr, Xi );
    arraysUnlinked = false;
    dft.evaluate();
  }
  
  
  
  /**
   * evaluates the inverse DFT with specified transform and sequence arrays
   * @param Xr          float array containing transform real part
   * @param Xi          float array containing transform imaginary part
   * @param xr          float array containing sequence real part 
   * @param xi          float array containing sequence imaginary part
   */
  public void evaluateInverse( float[] Xr, float[] Xi, float[] xr, float[] xi ) {
    this.yr = xr;
    this.yi = xi; 
    dft.link( Xr, Xi, xr, xi );
    arraysUnlinked = false;
    evaluateInverse(); 
  }
  

  
  /** 
   * constructs a CDFT instance with references to sequence and transform arrays
   * @param xr          float array containing sequence real part on forward evaluation,
   *                    transform real part on inverse evaluation
   * @param xi          float array containing sequence imaginary part on forward evaluation,
   *                    transform imaginary part on inverse evaluation
   * @param yr          float array containing transform real part on forward evaluation,
   *                    sequence real part on inverse evaluation
   * @param yi          float array containing transform imaginary part on forward evaluation,
   *                    sequence imaginary part on inverse evaluation
   * @param log2N       base-2 logarithm of the length of the transform
   */
  public CDFT( float[] xr, float[] xi, float[] yr, float[] yi, int log2N ) {
    
    if ( log2N < 3 ) throw new IllegalArgumentException( "DFT size must be >= 8" );
    
    this.yr = yr;
    this.yi = yi;

    this.log2N = log2N;
    N = 1 << log2N;

    createTable();

    if (      log2N == 3 )
      dft = new CDFTsr8(  0, 1, 0 );
    else if ( log2N == 4 )
      dft = new CDFTsr16( 0, 1, 0 );
    else if ( log2N >= 5 ) 
      dft = new CDFTsr( log2N, c, c3, s, s3 );
    
    dft.link(  xr, xi, yr, yi );
    arraysUnlinked = false;

  }
  


  /**
   *  evaluates the DFT assuming sequence and transformed arrays have been linked at construction time
   * 
   */
  public void evaluate() {
    if ( arraysUnlinked ) 
      throw new IllegalStateException( "Sequence and transform arrays are not linked" );
    dft.evaluate();
  }

  
  
  
  /**
   * evaluates the inverse DFT assuming the sequence and transform arrays have been linked at construction time
   */
  public void evaluateInverse() {
    
    if ( arraysUnlinked ) 
      throw new IllegalStateException( "Sequence and transform arrays are not linked" );
    
    dft.evaluate();
    
    float scale = 1.0f / (float) N;
    int N2 = N/2;
    
    yr[0]  *= scale;
    yi[0]  *= scale;
    yr[N2] *= scale;
    yi[N2] *= scale;
    
    int i = 1;  
    int j = N-1;
    
    float tmp;
    
    while ( i < j ) {
      tmp   = yr[i];
      yr[i] = yr[j]*scale;
      yr[j] = tmp*scale;
      tmp   = yi[i];
      yi[i] = yi[j]*scale;
      yi[j] = tmp*scale;

      i++;
      j--;
    }
    
  }

  
  
  private void createTable() {
    
    int N8 = N/8;
    
    c  = new float[N8];
    c3 = new float[N8];
    s  = new float[N8];
    s3 = new float[N8];
    
    for ( int i = 0; i < N8; i++ ) {
      c[  i ] =  (float) Math.cos( 2 * Math.PI * i / N );
      c3[ i ] =  (float) Math.cos( 2 * Math.PI * 3 * i / N );
      s[  i ] = -(float) Math.sin( 2 * Math.PI * i / N );
      s3[ i ] = -(float) Math.sin( 2 * Math.PI * 3 * i / N );
    }

  }



  /**
   * Convenience method to multiply two complex transforms of the same size.
   * @param Xr     float array containing the real part of the first transform 
   * @param Xi     float array containing the imaginary part of the first transform 
   * @param Yr     float array containing the real part of the second transform before call, real part of the product after call
   * @param Yi     float array containing the imaginary part of the second transform before call, imaginary part of the product after call
   * @param sign   +1 for convolution type product, -1 for correlation type product
   */
  public static void dftProduct( float[] Xr, float[] Xi, float[] Yr, float[] Yi, float sign ) {
    
    if ( Xr.length != Yr.length  ||  Xi.length != Yi.length  ||  Xr.length != Xi.length )
      throw new IllegalArgumentException( "Transform array lengths are not equal" );
    
    float tmp;
    for ( int i = 0;  i < Xr.length;  i++ ) {
      tmp   = Xr[i]*Yr[i]  -  sign*Xi[i]*Yi[i];
      Yi[i] = Xr[i]*Yi[i]  +  sign*Xi[i]*Yr[i];
      Yr[i] = tmp;
    }
    
  }

}
