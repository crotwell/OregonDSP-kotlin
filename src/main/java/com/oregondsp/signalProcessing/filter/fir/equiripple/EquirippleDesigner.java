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

package com.oregondsp.signalProcessing.filter.fir.equiripple;

import java.util.ArrayList;



import com.oregondsp.signalProcessing.fft.RDFT;
import com.oregondsp.signalProcessing.filter.LagrangePolynomial;


/**
 * Implements the Parks-McClellan algorithm for designing equiripple FIR digital filters.
 * 
 * <p>See</p>
 * 
 * <p>Chebyshev Approximation for Nonrecursive Digital Filters with Linear Phase,
 * Thomas W. Parks and James H. McClellan (1972), IEEE Transactions on Circuit Theory, Vol. CT-19, no. 2,
 * pp. 184-194.</p>
 * 
 * <p>A Unified Approach to the Design of Optimum Linear Phase FIR Digital Filters,
 * James H. McClellan and Thomas W. Parks (1973), IEEE Transactions on Circuit Theory, Vol. CT-20, 
 * No. 6, pp. 697-701.</p>
 * 
 * <p>A Program for the Design of Linear Phase Finite Impulse Response Digital Filters,
 * Thomas W. Parks and James H. McClellan (1972), IEEE Transactions on Audio and Electroacoustics,
 * Vol. AU-20 no. 3, pp. 195-199.</p>
 * 
 *  @author David B. Harris,   Deschutes Signal Processing LLC
 */
class EquirippleDesigner {

  /** Constant specifying the maximum number of iterations of the Remez exchange algorithm */
  private static final int    MAXITER   = 25;
  
  
  /**
   * Remez exchange algorithm.
   *
   * @param G     DesignGrid object containing the finite frequency-sampling grid used by the Remez
   *                exchange algorithm.
   */
  static void remez( DesignGrid G ) {
    
    int nextrema = G.extremaIndices.length;
    
    ArrayList< Integer > newExtrema = new ArrayList< Integer >();
    double[] E  = new double[ G.gridSize ];
    double[] GA = new double[ G.gridSize ];
    
   
    int    niter = 0;
    
    do {
      
      double delta = computeDelta( G );
      
      System.out.println( "delta: " + delta );
       
      LagrangePolynomial LP = constructInterpolatingPolynomial( G, delta );
    
      //  Compute current approximant (GA) and error function (E) on grid
    
      for ( int i = 0;  i < G.gridSize;  i++ ) {
        GA[i] = LP.evaluate( G.X[i] );
        E[i]  = GA[i] - G.H[i];
      }
      
      // Search for new extrema starting from old extrema
      
      newExtrema.clear();
    
      int change = 0;
      
      for ( int currentExtremum = 0;  currentExtremum < nextrema;  currentExtremum++ ) {
        
        int currentGridPt = G.extremaIndices[ currentExtremum ];
        int s             = sgn( E[currentGridPt ] );
        
        // search forward
        
        int ptr = currentGridPt + 1;
        if ( ptr < G.gridSize ) {
          while( sgn( E[ptr] - E[ptr-1] )  ==  s ) {
            ptr++;
            if ( ptr == G.gridSize ) break;
          }
        }
        ptr--;
        
        if ( ptr == currentGridPt ) {
          
          // forward search failed, try backward search
          
          ptr = currentGridPt - 1;
          if ( ptr >= 0 ) {
            while( sgn( E[ptr] - E[ptr+1] )  ==  s ) {
              ptr--;
              if ( ptr < 0 ) break;
            }
          }
          ptr++;
          
        }
        
        newExtrema.add( ptr );
        if ( ptr != currentGridPt ) change++;
      }
      
      // test for exchanges at 0 and pi
      
      
      if ( G.containsZero  &&  G.containsPi ) {
        
        int gridPi = G.gridSize-1;
      
        if ( newExtrema.contains( 0 ) ) {
        
          if ( !newExtrema.contains( gridPi ) ) {
            if ( sgn( E[ gridPi ] )  !=  sgn( E[ G.extremaIndices[nextrema-1] ] ) ) {
              if ( Math.abs( E[ gridPi ] )  >  Math.abs( E[ 0 ] ) ) {
                newExtrema.remove( 0 );
                newExtrema.add( gridPi );
                change++;
              }
            }
          }
        }
      
        else {
        
          if ( newExtrema.contains( gridPi ) ) {
          
            if ( sgn( E[0] ) != sgn( E[ G.extremaIndices[0] ] ) ) {
              if ( Math.abs( E[ 0 ] )  >  Math.abs( E[ gridPi ] ) ) {
                newExtrema.remove( newExtrema.size()-1 );
                newExtrema.add( 0, 0 );
                change++;
              }
            }
          }
        
        }
        
      }
      
      if ( change == 0 ) break;

      // exchange extrema
    
      for ( int i = 0;  i < nextrema;  i++ ) {
        G.extremaIndices[ i ] = newExtrema.get(i);
      }
      
      niter++;
    } while ( niter < MAXITER );
    
  }
  
  
  
  /**
   * Method to compute the Linfinity norm best approximation error on the current set of extrema.
   *
   * @param G      DesignGrid instance, which contains a list of current extrema.
   * @return       double containing the error on the current set of extrema.
   */
  static double computeDelta( DesignGrid G ) {
    
    int nextrema = G.extremaIndices.length;
    double[] extrema = new double[ nextrema ];
    for ( int i = 0;  i < nextrema;  i++ ) {
      extrema[i] = G.X[ G.extremaIndices[i] ];
    }
    double[] gamma = LagrangePolynomial.BarycentricWeights( extrema );
    
    double num   = 0.0;
    double denom = 0.0;
    double s     = 1.0;
    for ( int i = 0;  i < nextrema;  i++ ) {
      int j = G.extremaIndices[i];
      num   += gamma[i]*G.H[j];
      denom += s*gamma[i]/G.W[j];
      s = -s;
    }
    
    return num/denom;
  }
  
  
  
  /**
   * Constructs the Lagrange polynomial on the set of extrema.
   *
   * @param  G        DesignGrid instance which contains the current set of extrema.
   * @param  delta    Current deviation on that set of extrema.
   * @return          Lagrange polynomial instance that interpolates the extrema.
   */
  static LagrangePolynomial constructInterpolatingPolynomial( DesignGrid G, double delta ) {
    
    double[] extremaSubset = new double[ G.extremaIndices.length - 1 ];
    int n = extremaSubset.length;
    double[] x = new double[ n ];
    double[] f = new double[ n ];
    double s = 1.0;
    for ( int i = 0;  i < n;  i++ ) {
      int j = G.extremaIndices[i];
      x[i] = G.X[ j ];
      f[i] = G.H[j] - s*delta/G.W[j];
      s = -s;
    }
    
    return new LagrangePolynomial( x, f );
  }
  
  
  
  /**
   * Calculates coefficients of the best Chebyshev approximation out of a cosine basis.
   *
   * @param   G    DesignGrid instance that contains a list of the current extrema. 
   * @param   Nc   The number of coefficients of the corresponding FIR filter.  
   * @return       float[] containing the coefficients of the FIR filter.
   */
  static float[] calculateCoefficients( DesignGrid G, int Nc ) {
    
    LagrangePolynomial LP = constructInterpolatingPolynomial( G, computeDelta( G ) );
   
    int     log2nfft = 6;
    int     nfft     = 64;
    while ( nfft < Nc ) {
      nfft *= 2;
      log2nfft++;
    }
    float[] X        = new float[ nfft ];
    float[] x        = new float[ nfft ];
    for ( int i = 0;  i <= nfft/2;  i++ ) {
      X[i] = (float) LP.evaluate( Math.cos( 2.0*Math.PI*i/nfft ) );
    }
    
    RDFT dft = new RDFT( log2nfft );
    dft.evaluateInverse( X, x );
    
    return x;
  }
  
  
  
  /**
   * Method to compute the sign of a double.
   *
   * @param x    the double in question
   * @return     int ( = 1 if x > 0, = 0 if x = 0, = -1 if x < 0 )
   */
  static int sgn( double x ) {
    if ( x > 0.0 ) 
      return 1;
    else if ( x < 0.0 ) 
      return -1;
    else 
      return 0;
  }

  
}
