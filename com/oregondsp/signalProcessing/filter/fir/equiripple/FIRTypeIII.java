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
 * Class for designing FIR type III digital filters.  Odd length filters with odd symmetry.
 * 
 * <p>See</p>
 * 
 * A Unified Approach to the Design of Optimum Linear Phase FIR Digital Filters,
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
abstract class FIRTypeIII extends EquirippleFIRFilter {

  
  /**
   * Instantiates a new FIR type III filter.
   *
   * @param numBands      int specifying the number of pass and stop bands.
   * @param nHalf         int specifying the half-length of the filter - equal to one less
   *                        than the number of approximating basis functions in the Remez algorithm.
   */
  FIRTypeIII( int numBands, int nHalf ) {
    
    super( numBands, nHalf, 2*nHalf+1 );
    
  }
  
  
  
  /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#populateGrid(com.oregondsp.signalProcessing.filter.fir.equiripple.DesignGrid)
   */
  void populateGrid( DesignGrid G ) {
    
    for ( int i = 0;  i < G.gridSize;  i++ ) {
      G.H[i]    = desiredResponse( G.grid[i] ) / Math.sin( G.grid[i]*Math.PI );
      G.W[i]    = weight( G.grid[i] ) * Math.sin( G.grid[i]*Math.PI );
    }

    G.containsZero = false;
    G.containsPi   = false;
  }
   
  
  
  /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#interpretCoefficients(float[])
   */
  float[] interpretCoefficients( float[] coefficients ) {
    float[] retval = new float[ Nc ];
    Sequence.circularShift( coefficients, N-1 );
    retval[0] = -0.5f*coefficients[0];
    retval[1] = -0.5f*coefficients[1];
    for ( int i = 2;  i < Nc-2;  i++ ) {
      retval[i] = 0.5f*( coefficients[i-2] - coefficients[i] );
    }
    retval[ Nc-2 ] = 0.5f*coefficients[ Nc-4 ];
    retval[ Nc-1 ] = 0.5f*coefficients[ Nc-3 ];
    return retval;
  }


}
