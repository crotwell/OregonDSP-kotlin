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
 * Implements an even-length differentiator - the point of symmetry falls between samples.

 *
 * This class implements a differentiator with impulse response that exhibits odd symmetry on a
 * staggered grid, i.e. the point of symmetry falls "between" samples.  This filter has an accurate
 * differentiator response that extends from 0 to the folding frequency.  This advantage may be
 * offset by a non-integer group delay ( (2N-1)/2 ) which shifts the waveform by a fractional sample.
 * The quality of the design is controlled by a single parameter (order N), which specifies the number
 * of approximating basis functions in the Remez algorithm.

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
@JsExport
class StaggeredDifferentiator
/**
 * Instantiates a new differentiator.

 * @param N       int specifying the filter order design parameter.  The larger this value, the
 * *                  more accurate the differentiator response.
 * *
 * @param delta   double specifying the sampling interval of the data to be differentiated.
 */
(N: Int,
 /** double representing the sampling interval of the data to be differentiated.  */
 private val delta: Double) : FIRTypeIV(1, N) {


    init {

        bands[0][0] = 1.0 / (2 * N)
        bands[0][1] = 1.0

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
