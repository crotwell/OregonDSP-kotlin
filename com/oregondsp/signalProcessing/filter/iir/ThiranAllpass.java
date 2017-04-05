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

package com.oregondsp.signalProcessing.filter.iir;


import com.oregondsp.signalProcessing.filter.Polynomial;


/**
 * Designs and implements Thiran allpass filters.
 * 
 * Thiran allpass filters are used to interpolate signals by a fractional sample.  
 * They have unit amplitude response, thus have no amplitude distortion, and approximate
 * a flat group delay specified by D.  The group delay function is maximally flat at 0 Hz.
 * 
 * @author David B. Harris   Deschutes Signal Processing LLC
 *
 */
public class ThiranAllpass extends Allpass {

  
  /**
   * constructs a Thiran allpass filter.
   *
   * @param N     the order of the allpass filter, typically 3 or 4
   * @param D     the delay, in samples, best between N-1 and N
   */
  public ThiranAllpass( int N, double D ) {
    
    super( N );
    
    double[] a  = new double[N+1];
    
    a[0]        = 1.0;
    for ( int i = 1;  i <= N;  i++ ) {
      double prod = 1.0;
      for ( int n = 0;  n <= N;  n++ ) {
        prod *= ( (double)(D - N + n ) ) / ( (double) (D - N + i + n ) );
      }
      a[i] = Math.pow( -1, i ) * ( factorial(N) / ( factorial(N-i) * factorial(i) ) ) * prod;
    }
    
    Polynomial P = new Polynomial( a );
    k            = P.reflectionCoefficients();
    constructRationalRepresentation();
    
  }
    

  /**
   * Factorial function required to evaluate the coefficients of the Thiran allpass filter.
   *
   * @param n       int argument of the factorial function.
   * @return        int n!
   */
  private int factorial( int n ) {

    int retval = 1;
    if ( n > 1 )
      for ( int i = 2;  i <= n;  i++ ) retval *= i;
    
    return retval;
  }

}
