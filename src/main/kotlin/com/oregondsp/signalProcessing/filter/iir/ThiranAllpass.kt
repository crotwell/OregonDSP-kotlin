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

package com.oregondsp.signalProcessing.filter.iir


import com.oregondsp.signalProcessing.filter.Polynomial
import kotlin.js.Math


/**
 * Designs and implements Thiran allpass filters.

 * Thiran allpass filters are used to interpolate signals by a fractional sample.
 * They have unit amplitude response, thus have no amplitude distortion, and approximate
 * a flat group delay specified by D.  The group delay function is maximally flat at 0 Hz.

 * @author David B. Harris   Deschutes Signal Processing LLC
 */
class ThiranAllpass
/**
 * constructs a Thiran allpass filter.

 * @param N     the order of the allpass filter, typically 3 or 4
 * *
 * @param D     the delay, in samples, best between N-1 and N
 */
(N: Int, D: Double) : Allpass(N) {


    init {

        val a = DoubleArray(N + 1)

        a[0] = 1.0
        for (i in 1..N) {
            var prod = 1.0
            for (n in 0..N) {
                prod *= (D - N + n).toDouble() / (D - N + i.toDouble() + n.toDouble()).toDouble()
            }
            a[i] = Math.pow(-1.0, i.toDouble()) * (factorial(N) / (factorial(N - i) * factorial(i))).toDouble() * prod
        }

        val P = Polynomial(a)
        k = P.reflectionCoefficients()
        constructRationalRepresentation()

    }


    /**
     * Factorial function required to evaluate the coefficients of the Thiran allpass filter.

     * @param n       int argument of the factorial function.
     * *
     * @return        int n!
     */
    private fun factorial(n: Int): Int {

        var retval = 1
        if (n > 1)
            for (i in 2..n) retval *= i

        return retval
    }

}
