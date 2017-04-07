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
 * Implements a centered equiripple Hilbert transform operator - the point of symmetry falls on a sample.

 *
 * This class uses the Remez exchange algorithm to design a Hilbert transformer of length 2N+1.
 * The Nth sample is zero (counting from 0), and the impulse response is an anti-symmetric sequence
 * about this point.  The filter is linear phase, with group delay a constant equal to N.  The
 * design parameters are the _order (N) specifying the number (N+1) of approximating functions in the
 * Remez algorithm, and two parameters specifying the band edge frequencies.  The design problem
 * is performed on a discrete-time frequency axis normalized to range between 0 and 1 (the folding
 * frequency).  Omega1, the lower band edge of the passband must be greater than 0 and less than
 * Omega2, the upper band edge.  Omega2 must be strictly less than 1.  A tradeoff exists between the
 * filter _order N and band edge frequencies.  As Omega1 approaches 0 or Omega2 approaches 1, the _order
 * must be increased to obtain an acceptable design.
 *
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

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
class CenteredHilbertTransform
/**
 * Instantiates a new centered Hilbert transform operator.

 * @param N       int specifying the number (N+1) of approximating functions in the Remez design
 * *                  algorithm and the resulting number of FIR filter coefficients (2N+1).
 * *
 * @param Omega1  double specifying the low passband edge of the filter.  Omega1 > 0
 * *
 * @param Omega2  double specifying the high passband edge of the filter. Omega1 < Omega2 < 1.
 */
(N: Int, Omega1: Double, Omega2: Double) : FIRTypeIII(1, N) {


    init {

        if (!(0 < Omega1 && Omega1 < Omega2 && Omega2 < 1.0))
            throw IllegalArgumentException("Check 0.0 < Omega1 < Omega2 < 1.0")

        bands[0][0] = Omega1
        bands[0][1] = Omega2

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
            retval = 1.0

        return retval
    }

}
