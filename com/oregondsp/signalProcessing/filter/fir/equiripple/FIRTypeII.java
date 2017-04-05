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

import com.oregondsp.signalProcessing.Sequence;


/**
 * Class for designing FIR type II digital filters.  Even length filters with even symmetry.
 * 
 * <p>See </p>
 * 
 * <p>A Unified Approach to the Design of Optimum Linear Phase FIR Digital Filters,
 * James H. McClellan and Thomas W. Parks (1973), IEEE Transactions on Circuit Theory, Vol. CT-20, 
 * No. 6, pp. 697-701.</p>
 * 
 * <p>and</p>
 * 
 * <p>FIR Digital Filter Design Techniques Using Weighted Chebyshev Approximation, 
 * Lawrence R. Rabiner, James H. McClellan and Thomas W. Parks (1975) PROCEEDINGS OF THE IEEE,
 * VOL. 63, NO. 4, pp. 595-610.</p>
 * 
 *  @author David B. Harris,   Deschutes Signal Processing LLC
 */
abstract class FIRTypeII extends EquirippleFIRFilter {
  
  // even length filters with even symmetry

  /**
   * Instantiates a new fIR type ii.
   *
   * @param numBands     int specifying the number of pass and stop bands.
   * @param nHalf        int specifying the half size of the filter - equal to the number of
   *                       approximating basis functions in this case.
   */
  FIRTypeII( int numBands, int nHalf ) {
    
    super( numBands, nHalf, 2*nHalf );
    
  }

  
  
  /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#populateGrid(com.oregondsp.signalProcessing.filter.fir.equiripple.DesignGrid)
   */
  void populateGrid( DesignGrid G ) {
    
    for ( int i = 0;  i < G.gridSize;  i++ ) {
      G.H[i]    = desiredResponse( G.grid[i] ) / Math.cos( G.grid[i]*Math.PI/2.0 );
      G.W[i]    = weight( G.grid[i] ) * Math.cos( G.grid[i]*Math.PI/2.0 );
    }
    
    if ( Math.abs( G.grid[0] ) < 1.0E-6 )
      G.containsZero = true;
    else
      G.containsZero = false;
    
    G.containsPi   = false;
   }
   
  
  
  /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#interpretCoefficients(float[])
   */
  float[] interpretCoefficients( float[] coefficients ) {
    float[] retval = new float[ Nc ];
    Sequence.circularShift( coefficients, N-1 );
    retval[0] = 0.5f*coefficients[0];
    for ( int i = 1;  i < Nc-1;  i++ ) {
      retval[i] = 0.5f*( coefficients[i] + coefficients[i-1] );
    }
    retval[ Nc-1 ] = 0.5f*coefficients[ Nc-2 ];
    return retval;
  }

}
