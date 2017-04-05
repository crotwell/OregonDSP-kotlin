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



/**
 * Half band prototype filter used by EquirippleHalfBand for interpolation by a factor of 2.
 * 
 * <p>This class is intended be used by EquirippleHalfBand to design filters suitable for
 * interpolating sequences by a factor of two.  It is used in conjuction with the half-band "trick"
 * described in:</p>
 * 
 * <p>A &#8220;TRICK&#8221; for the Design of FIR Half-Band Filters, P. P. VAIDYANATHAN AND TRUONG Q. NGUYEN (1987),
 * IEEE TRANSACTIONS ON CIRCUITS AND SYSTEMS, VOL. CAS-34, NO. 3, pp. 297-300.</p>
 * 
 *  @author David B. Harris,   Deschutes Signal Processing LLC
 */
class EquirippleHalfBandPrototype extends FIRTypeII {


  /**
   * Instantiates a new equiripple half band prototype.
   *
   * @param N         int specifying the design order of the filter.
   * @param OmegaP    double specifying the upper band edge of the single band used in this filter type.
   */
  EquirippleHalfBandPrototype( int N, double OmegaP ) {
    
    super( 1, N );
    
    if ( OmegaP <= 0.0  ||  OmegaP >= 1.0 ) 
      throw new IllegalArgumentException( "OmegaP: " + OmegaP + " out of bounds (0.0 < OmegaP < 1.0)" );

    
    bands[0][0] = 0.0;
    bands[0][1] = OmegaP;
    
    generateCoefficients();
  }


  
  /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#desiredResponse(double)
   */
  double desiredResponse( double Omega ) {
    
    double retval = 0.0;
    if ( LTE( bands[0][0], Omega)  &&  LTE( Omega, bands[0][1] ) )  retval = 1.0;
      
    return retval;
  }



  /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#weight(double)
   */
  double weight( double Omega ) {
    
    double retval = 0.0;
    
    if (  LTE( bands[0][0], Omega)  &&  LTE( Omega, bands[0][1] ) ) 
      retval = 1.0;
    
    return retval;
  }
  
}
