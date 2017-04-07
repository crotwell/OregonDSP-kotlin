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
 * Designs a half-band FIR equiripple filter using the "half-band trick" and the Remez algorithm.

 *
 * This class designs a half-band filter (a filter suitable for interpolating data by a factor of 2).
 * It uses the "half-band trick" described in

 *
 * A &#8220;TRICK&#8221; for the Design of FIR Half-Band Filters, P. P. VAIDYANATHAN AND TRUONG Q. NGUYEN (1987),
 * IEEE TRANSACTIONS ON CIRCUITS AND SYSTEMS, VOL. CAS-34, NO. 3, pp. 297-300.

 *
 * The filter is obtained as a transformation of a EquirippleHalfBandPrototype, which is designed
 * with the Remez exchange algorithm.  As with other filters, the prototype is specified by a design
 * order parameter (N) and a band edge parameter (OmegaP).  The resulting FIR filter has 2N+1 coefficients
 * and is evenly symmetric about coefficient N, counting from 0.  The band edge parameter should be close
 * to 0.5, though slightly less: 0 < OmegaP < 0.5.  A value of 0.45 is reasonable, and the closer OmegaP
 * is to 0.5, the larger N must be to obtain a reasonable response.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
class EquirippleHalfBand
/**
 * Instantiates a new equiripple half band filter.

 * @param N         int specifying the design order.
 * *
 * @param OmegaP    double specifying the upper passband cutoff.
 */
(N: Int, OmegaP: Double) {

    /** float[] containing the FIR filter coefficients.  */
    private val coefficients: FloatArray


    init {

        val EHBP = EquirippleHalfBandPrototype(N, 2 * OmegaP)

        val c = EHBP.getCoefficients()

        coefficients = FloatArray(2 * c.size - 1)
        for (i in c.indices) {
            coefficients[2 * i] = 0.5f * c[i]
        }
        coefficients[c.size - 1] = 0.5f

    }


    /**
     * Accessor for the FIR filter coefficients.

     * @return    float[] containing the filter coefficients.
     */
    fun getCoefficients(): FloatArray {
        return coefficients.copyOf()
    }

}
