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

package com.oregondsp.signalProcessing.filter.fir.equiripple


/**
 * Implements a centered FIR lowpass filter - the point of symmetry falls on a sample.

 *
 * This class uses the Remez exchange algorithm to design a lowpass filter of length 2N+1.
 * The impulse response is a symmetric sequence about point N (counting from 0). The filter is linear
 * phase, with group delay a constant equal to N.  The design parameters are the _order (N) specifying
 * the number (N+1) of approximating functions in the Remez algorithm and four parameters controlling
 * the cutoffs and design weights of the passband and stopband.  The design problem is performed on
 * a discrete-time frequency axis normalized to range between 0 and 1 (the folding frequency).  The
 * pass band is the interval [0, OmegaP] and the stop band is the interval [OmegaS, 1].  Note that
 * OmegaP < OmegaS and the two bands must have non-zero width.  There also is a transition band,
 * the open interval (OmegaP, OmegaS), that must have non-zero width.  The narrower any of these bands,
 * the larger the _order N must be to obtain a reasonable frequency response.  Weights are specified
 * for each band to control the relative size of the maximum error between bands.
 * For details on the design algorithm and characteristics of the filter response, see

 *
 * A Unified Approach to the Design of Optimum Linear Phase FIR Digital Filters,
 * James H. McClellan and Thomas W. Parks (1973), IEEE Transactions on Circuit Theory, Vol. CT-20,
 * No. 6, pp. 697-701.

 *
 * and

 *
 * FIR Digital Filter Design Techniques Using Weighted Chebyshev Approximation,
 * Lawrence R. Rabiner, James H. McClellan and Thomas W. Parks (1975) PROCEEDINGS OF THE IEEE,
 * VOL. 63, NO. 4, pp. 595-610.

 *
 * and for _order selection, consult:

 *
 * Approximate Design Relationships for Low-Pass FIR Digital Filters, Lawrence R. Rabiner (1973),
 * IEEE TRANSACTIONS ON AUDIO AND ELECTROACOUSTICS, VOL. AU-21, NO. 5, pp. 456-460.


 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
class EquirippleLowpass
/**
 * Instantiates a new equiripple lowpass filter object.

 * @param N        int specifying the design _order of the filter.
 * *
 * @param OmegaP   double specifying the passband upper cutoff frequency.
 * *
 * @param Wp the   double specifying the passband weight.
 * *
 * @param OmegaS   double specifying the stopband lower cutoff frequency.
 * *
 * @param Ws       double specifying the stopband weight.
 */
(N: Int, OmegaP: Double,
 /** double specifying the passband weight.  */
 private val Wp: Double, OmegaS: Double,
 /** double specifying the stopband weight.  */
 private val Ws: Double) : FIRTypeI(2, N) {

    init {

        if (OmegaP >= OmegaS) throw IllegalArgumentException("OmegaP >= OmegaS ")
        if (OmegaP <= 0.0 || OmegaP >= 1.0)
            throw IllegalArgumentException("OmegaP: $OmegaP out of bounds (0.0 < OmegaP < 1.0)")
        if (OmegaS <= 0.0 || OmegaS >= 1.0)
            throw IllegalArgumentException("OmegaS: $OmegaS out of bounds (0.0 < OmegaS < 1.0)")

        bands[0][0] = 0.0
        bands[0][1] = OmegaP
        bands[1][0] = OmegaS
        bands[1][1] = 1.0

        generateCoefficients()
    }


    /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#desiredResponse(double)
   */
    internal override fun desiredResponse(Omega: Double): Double {

        var retval = 0.0
        if (LTE(bands[0][0], Omega) && LTE(Omega, bands[0][1])) retval = 1.0

        return retval
    }


    /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#weight(double)
   */
    internal override fun weight(Omega: Double): Double {

        var retval = 0.0

        if (LTE(bands[0][0], Omega) && LTE(Omega, bands[0][1]))
            retval = Wp
        else if (LTE(bands[1][0], Omega) && LTE(Omega, bands[1][1]))
            retval = Ws

        return retval
    }

}
