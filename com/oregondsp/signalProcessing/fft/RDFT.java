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
 * 
 *  Class to calculate the discrete Fourier transform of a real sequence using symmetries and a complex DFT of half the length.
 *  
 *  <p>This class is designed for efficient calculation of many discrete Fourier transforms of the
 *  same length.  It is limited to transform lengths that are powers of two and are greater than or 
 *  equal to 32.  It uses the identity:</p>  
 *  
 *  <p>
 *  x[2n] + j*x[2n+1]  <->  X[k] + X[k+N/2] + j*W<sup>k</sup>*( X[k] - X[k+N/2] )
 *  </p>
 *  <p>
 *  and the fact that, for real sequences x[n], the transform exhibits conjugate symmetry:
 *  </p>
 *  <p>X[k] = conjg( X[N-k] )</p>
 *  <p>
 *  to calculate a real DFT with a complex DFT of half the length.
 *  </p>
 *  
 *  <p> Example of use: </p>
 *  <p>
 *  <font face="courier">
 *   int N &nbsp&nbsp&nbsp&nbsp= 16384;<BR>
 *   int log2N = 14;<BR>
 *   float[] x = new float[N];<BR>
 *   float[] X = new float[N];<BR>
 *   RDFT Xfm = new RDFT( log2N );<BR>
 *   <BR>
 *   // load data<BR>
 *   for ( int i = 0;  i < N;  i++ ) {<BR>
 *     &nbsp x[i] = ...<BR>
 *   }<BR>
 *   <BR>
 *   // evaluate transform of data<BR>
 *   Xfm.evaluate( x, X );<BR>
 *  </font>
 *  </p>
 *  <p>
 *  The input sequence, x[n], is in natural order, and the output transform is in the packed form used
 *  by Sorensen et al. (1987) for conjugate symmetric transforms: <p>
 * 
 *  <p>
 *  <font face="courier">
 *   __0_____1_____2_____3_____..._______N/2-1______N/2_______N/2+1______N/2+2____...____N-1 <BR>
 *   Xr(0)_Xr(1)_Xr(2)_Xr(3)___..._____Xr(N/2-1)__Xr(N/2)___Xi(N/2-1)__Xi(N/2-2)__...___Xi(1)
 *   </font>
 *   </p>
 *   <p> where Xr is the real part of the transform and Xi is the imaginary part.</p>
 *  
 *  <p>As long as the transform size does not change, the RDFT object does not need to be reinstantiated.  
 *  Consequently, the data arrays can be reloaded and the evaluate method invoked to compute additional 
 *  DFTs without incurring the cost of RDFT object instantiation.</p>
 *  
 *  <p>
 *  The inverse DFT is calculated with a call to evaluateInverse():
 *  </p>
 *  <p>
 *  <font face="courier">
 *  Xfm.evaluateInverse( X, x );
 *   </font>
 *   </p> 
 *
 * <p>
 * The input transform, X, is in conjugate symmetric packed form and the output sequence, x, is in natural order.
 * </p>
 * <p>
 * This class also contains a convenience method to support convolution and correlation via the FFT.
 * </p>
 *  
 *  <p> See "Real-valued Fast Fourier Transform Algorithms", Sorensen, H. V., et al., IEEE TRANSACTIONS ON 
 *  ACOUSTICS, SPEECH, AND SIGNAL PROCESSING, VOL. ASSP-35, NO. 6, JUNE 1987, pp. 849-863.</p>
 *  
 * @author David B. Harris,  Deschutes Signal Processing LLC
 *
 */
public class RDFT {
  
  private int     N;
  private int     N2;
  private int     N4;
  
  private float[] xr;
  private float[] xi;
  private float[] Xr;
  private float[] Xi;
  
  private CDFT    dft;
  
  private float[] c;
  private float[] s;
  
  
  public RDFT( int log2N ) {
    
    if ( log2N < 4 ) throw new IllegalArgumentException( "DFT size must be >= 16" );
    
    N  = 1 << log2N;
    N2 = N/2;
    N4 = N/4;
    
    xr = new float[ N2 ];
    xi = new float[ N2 ];
    Xr = new float[ N2 ];
    Xi = new float[ N2 ];
    
    s  = new float[ N4 ];
    c  = new float[ N4 ];
    
    for ( int i = 0;  i < N4;  i++ ) {
      s[i] = (float)  Math.sin( 2.0*Math.PI/N * i );
      c[i] = (float)  Math.cos( 2.0*Math.PI/N * i );
    }
    
    
    dft = new CDFT( log2N-1 );
    
  }
  
  
  
  /**
   * Evaluates the DFT of a real sequence x.
   * @param x     float[] containing the real sequence in natural order.
   * @param X     float[] containing the transform of the sequence in conjugate symmetric packed form.
   */
  public void evaluate( float[] x, float[] X ) {
     
    // Uses symmetries to perform the real length-N DFT with a special length-N set of butterflies
    // and one length-N/2 complex DFT.
    //
    // uses, specifically, the identity:  x[2n] + j*x[2n+1]  <->  X[k] + X[k+N/2] + j*W^k*( X[k] - X[k+N/2] )
    // and the fact that, for real sequences, X[k] = conjg( X[N-k] )
    
    for ( int i = 0;  i < N2;  i++ ) {
      int j = i << 1;
      xr[i] = x[j++];
      xi[i] = x[j];
    }
    
    dft.evaluate( xr, xi, Xr, Xi );
    
    // special case at k = 0
    
    X[0]  = Xr[0] + Xi[0];
    X[N2] = Xr[0] - Xi[0];
    
    // 1 <= k < N/4
    
    int N2pk = N2 + 1;
    int N2mk = N2 - 1;
    int Nmk  = N - 1;
    for ( int k = 1;  k < N4;  k++ ) {
      
      float Xrk = Xr[k];
      float Xik = Xi[k];
      float XrN2mk = Xr[N2mk];
      float XiN2mk = Xi[N2mk];
      
      float Sr = ( Xrk + XrN2mk )/2;
      float Si = ( Xik - XiN2mk )/2;
      
      float Dr = ( Xik + XiN2mk )/2;
      float Di = ( XrN2mk - Xrk )/2;
      
      float tmp = c[k]*Dr + s[k]*Di;
      Di        = c[k]*Di - s[k]*Dr;
      Dr        = tmp;
      
      X[k]      = Sr + Dr;
      X[Nmk]    = Si + Di;
      
      X[N2mk]   = Sr - Dr;
      X[N2pk]   = Di - Si;
      
      N2pk++;
      N2mk--;
      Nmk--;
    }
    
    // special case at k = N/4
    
    //  cos( 2*pi/N * k ) = cos( pi/2 ) = 0
    //  sin( 2*pi/N * k ) = sin( pi/2 ) = 1
    
    X[N4]    =  Xr[N4];
    X[N2+N4] = -Xi[N4];
     
  }   
  
  

  /**
   * Evaluates the inverse DFT of a conjugate symmetric transform.
   * @param X     float[] containing the input transform of the sequence in conjugate symmetric packed form.
   * @param x     float[] containing the output real sequence in natural order.
   */
  public void evaluateInverse( float[] X, float[] x ) {
    
    // Assumed input storage:
    //   0     1     2     3     ...       N/2-1      N/2       N/2+1      N/2+2    ...    N-1
    // Xr(0) Xr(1) Xr(2) Xr(3)   ...     Xr(N/2-1)  Xr(N/2)   Xi(N/2-1)  Xi(N/2-2)        Xi(1)
    
    
    // Uses symmetries to perform the real length-N inverse DFT with a special length-N set of butterflies
    // and one length-N/2 complex DFT.
    //
    // uses, specifically, the identity:  x[2n] + j*x[2n+1]  <->  X[k] + X[k+N/2] + j*W^k*( X[k] - X[k+N/2] )
    // and the fact that, for real sequences, X[k] = conjg( X[N-k] )
    
    
    // special case at k = 0
    
    Xr[0] = X[0] + X[N2];
    Xi[0] = X[0] - X[N2];
    
    // 1 <= k < N/4
    
    int N2pk = N2 + 1;
    int N2mk = N2 - 1;
    int Nmk  = N - 1;
    for ( int k = 1;  k < N4;  k++ ) {
      
      float Xrk    =  X[k];
      float Xik    =  X[Nmk];
      float XrkpN2 =  X[N2mk];
      float XikpN2 = -X[N2pk];
      
      float Dr = Xrk - XrkpN2;
      float Di = Xik - XikpN2;
      
      Xr[k] = Xrk + XrkpN2  -  s[k]*Dr  -  c[k]*Di;
      Xi[k] = Xik + XikpN2  +  c[k]*Dr  -  s[k]*Di;
      
      N2pk++;
      N2mk--;
      Nmk--;
    }
    
    // special case at k = N/4
    
    //  cos( 2*pi/N * k ) = cos( pi/2 ) = 0
    //  sin( 2*pi/N * k ) = sin( pi/2 ) = 1
    
    Xr[N4] =  2.0f*X[N4];
    Xi[N4] = -2.0f*X[N2+N4];
    
    // N/4 + 1  <=  k  <  N/2 - 1;
    
    //  cos( 2*pi/N * (N/4+m) ) = cos( 2*pi/N*m + pi/2 ) = -cos( 2*pi/N*(N/4-m) )
    //  sin( 2*pi/N * (N/4+m) ) = sin( 2*pi/N * (N/4-m) )
    
    N2pk = N2 + N4 + 1;
    N2mk = N4 - 1;
    Nmk  = N - N4 - 1;
    int reflect = N4 - 1;
    for ( int k = N4 + 1;  k < N2;  k++ ) {
      
      float Xrk    =  X[k];
      float Xik    =  X[Nmk];
      float XrkpN2 =  X[N2mk];
      float XikpN2 = -X[N2pk];
      
      float Dr = Xrk - XrkpN2;
      float Di = Xik - XikpN2;

      Xr[k] = Xrk + XrkpN2  -  s[reflect]*Dr  +  c[reflect]*Di;
      Xi[k] = Xik + XikpN2  -  c[reflect]*Dr  -  s[reflect]*Di;
      
      N2pk++;
      N2mk--;
      Nmk--;
      reflect--;
    }
    
    dft.evaluate( Xr, Xi, xr, xi );
    
    x[0] = xr[0]/N;
    x[1] = xi[0]/N;
    
    int j = N2-1;
    for ( int k = 1;  k < N2;  k++ ) {
      int i = k << 1;
      x[i++] = xr[j]/N;
      x[i]   = xi[j]/N;
      j--;
    }
    
  } 
  
  
  
  /**
   * Calculates the product of two conjugate symmetric dfts of the same length and stores the result in the second dft.
   * 
   * <p>This is a convenience method to support convolution and correlation operations using the Fast Fourier Transform.
   *   Example of use to compute the convolution of two sequences: </p>
 *  <p>
 *  <font face="courier">
 *   int N &nbsp&nbsp&nbsp&nbsp= 1024;<BR>
 *   int log2N = 10;<BR>
 *   float[] x = new float[N];<BR>
 *   float[] y = new float[N];<BR> 
 *   float[] X = new float[N];<BR>
 *   float[] Y = new float[N];<BR>
 *   RDFT Xfm = new RDFT( log2N );<BR>
 *   <BR>
 *   // load sequences<BR>
 *   for ( int i = 0;  i < N;  i++ ) {<BR>
 *     &nbsp x[i] = ...<BR>
 *     &nbsp y[i] = ...<BR>
 *   }<BR>
 *   <BR>
 *   // evaluate transforms of sequences<BR>
 *   Xfm.evaluate( x, X );<BR>
 *   Xfm.evaluate( y, Y );<BR>
 *   <BR>
 *   // product of transforms <BR>
 *   RDFT.dftProduct( X, Y, 1 );<BR>
 *   <BR>
 *   // inverse transform to obtain convolution<BR>
 *   float[] xy = new float[1024];<BR>
 *   RDFT.evaluateInverse( Y, xy );<BR>
 *  </font>
 *  </p>
   *
   * @param kernel       first DFT
   * @param transform    second DFT before call, contains product after the call
   * @param sign         +1 if a convolution type product, -1 if a correlation type product
   */
  public static void dftProduct( float[] kernel, float[] transform, float sign ) {
    
    if ( kernel.length != transform.length )
      throw new IllegalArgumentException( "kernel and transform arrays must have the same size" );
    
    int n    = kernel.length;
    int half = n/2;
    transform[0]    *= kernel[0];
    transform[half] *= kernel[half];
    
    float tmp;
    for ( int i = 1;  i < half;  i++ ) {
      int im = n-i;
      tmp           = kernel[i]*transform[i]  - sign*kernel[im]*transform[im];
      transform[im] = kernel[i]*transform[im] + sign*kernel[im]*transform[i];
      transform[i]  = tmp;
    }
    
  }
  
}
