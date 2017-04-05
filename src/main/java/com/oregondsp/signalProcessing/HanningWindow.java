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

package com.oregondsp.signalProcessing;


/**
 * Designs and implements Hanning windows (See Oppenheim and Schafer, 1975).
 * 
 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
public class HanningWindow extends Window {
  
  /**
   * Instantiates a new Hanning window of length N samples.
   *
   * @param N    int specifying the window length in samples.
   */
  public HanningWindow( int N ) {
    
    super(N);
    
    for ( int i = 0;  i < N;  i++ ) {
      w[i] = (float) (0.5 + 0.5*Math.cos( -Math.PI + i*2*Math.PI/(N-1) ) );
    }
    
  }

}
