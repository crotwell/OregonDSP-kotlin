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

package com.oregondsp.signalProcessing.filter.fir;

import com.oregondsp.signalProcessing.Sequence;
import com.oregondsp.signalProcessing.filter.fir.equiripple.CenteredHilbertTransform;


/**
 * Enables a complex analytic signal to be constructed from a real signal.
 * 
 * This class uses the CenteredHilbertTransform class to construct the complex analytic counterpart
 * of a real signal.  The class is perhaps most useful for obtaining the envelope of a signal and
 * a method is supplied for this purpose.  This class is intended to manipulate finite duration 
 * signals in one piece, not to process continuous streams in consecutive, contiguous blocks.
 * 
 *  @author David B. Harris,   Deschutes Signal Processing LLC
 */
public class ComplexAnalyticSignal {
  
  /** The real part of the signal. */
  float[] realPart;
  
  /** The imaginary part of the signal. */
  float[] imagPart;
  
  
  
  /**
   * Instantiates a new complex analytic signal.
   *
   * @param realSignal   float[] containing the original real signal.
   */
  public ComplexAnalyticSignal( float[] realSignal ) {
    realPart = realSignal.clone();
    CenteredHilbertTransform transformer = new CenteredHilbertTransform( 50, 0.03, 0.97 );
    float[] tmp = transformer.filter( realPart );
    Sequence.zeroShift( tmp, -50 );
    imagPart = new float[ realPart.length ];
    System.arraycopy( tmp, 0, imagPart, 0, realPart.length );
  }
  
  
  
  /**
   * Computes and returns the envelope of the signal.
   *
   * @return     float[] containing the signal envelope.
   */
  public float[] getEnvelope() {
    float[] retval = new float[ realPart.length ];
    for ( int i = 0;  i < realPart.length;  i++ ) {
      retval[i] = (float) Math.sqrt( realPart[i]*realPart[i] + imagPart[i]*imagPart[i] );
    }
    
    return retval;
  }
  
  
  
  /**
   * Accessor for the real part of the signal.
   *
   * @return     float[] containing the real part of the signal.
   */
  float[] getRealPart() {
    return realPart.clone();
  }
  
  
  
  /**
   * Accessor for the imaginary part of the signal.
   *
   * @return     float[] containing the imaginary part of the signal.
   */
  float[] getImagPart() {
    return imagPart.clone();
  }
  
}
