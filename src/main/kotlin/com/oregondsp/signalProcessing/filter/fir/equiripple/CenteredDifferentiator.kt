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

import kotlin.math.*


/**
 * Implements a centered equiripple differentiator - the point of symmetry falls on a sample.

 *
 * This class uses the Remez exchange algorithm to design a differentiator of length 2N+1.
 * The Nth sample (counting from 0) is zero, and the impulse response is an anti-symmetric
 * sequence about this point.  The filter is linear phase, with group delay a constant equal
 * to N.  The design parameters are the _order (N) specifying the number (N+1) of approximating
 * functions in the Remez algorithm, the intended sampling interval of the data and the passband
 * edge frequency OmegaP.  The design problem is performed on a discrete-time frequency axis
 * normalized to range between 0 and 1 (the folding frequency).  OmegaP, the upper band edge of
 * the passband must be strictly less than 1, and should be in the range 0.8 - 0.95, depending
 * on the specified _order.  The larger the _order, the closer OmegaP can be to 1.0 and yield an
 * acceptable approximation.  For details on the design algorithm and characteristics of the filter
 * response, see

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

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
class CenteredDifferentiator
/**
 * Instantiates a new centered differentiator.

 * @param N       int specifying the order of the filter, specifying the number of approximating functions
 * *                (N+1) and the number of resulting FIR filter coefficients (2N+1).
 * *
 * @param delta   double specifying the intended sampling interval of the data in seconds.
 * *
 * @param OmegaP  double specifying the upper passband cutoff (0 < OmegaP < 1).  It should be in the range
 * *                0.8 - 0.95+ with larger values of N required to obtain good approximants when OmegaP
 * *                approaches 1.
 */
(N: Int,
 /** Intended sampling interval of the data in seconds.  */
 private val delta: Double, OmegaP: Double) : FIRTypeIII(1, N) {


    init {

        if (!(0.0 < OmegaP && OmegaP < 1.0))
            throw IllegalArgumentException("Check 0.0 < OmegaP < 1.0")

        bands[0][0] = 1.0 / (2 * N)
        bands[0][1] = OmegaP

        generateCoefficients()
    }


    /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#desiredResponse(double)
   */
    internal override fun desiredResponse(Omega: Double): Double {

        var retval = 0.0
        if (LTE(bands[0][0], Omega) && LTE(Omega, bands[0][1])) retval = -PI * Omega / delta

        return retval
    }


    /* (non-Javadoc)
   * @see com.oregondsp.signalProcessing.filter.fir.equiripple.EquirippleFIRFilter#weight(double)
   */
    internal override fun weight(Omega: Double): Double {

        var retval = 0.0

        if (LTE(bands[0][0], Omega) && LTE(Omega, bands[0][1]))
            retval = 1.0 / Omega

        return retval
    }

}
