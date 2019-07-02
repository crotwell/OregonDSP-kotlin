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

import kotlin.math.*

/**
 * Class to implement a second _order section - basic unit of an Infinite Impulse Response digital filter.

 *
 * Implements the finite difference equation:

 *
 * y[n] = -a[1]*y[n-1] - a[2]*y[n-2] + b[0]*x[n] + b[1]*x[n-1] + b[2]*x[n-2]

 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
class SecondOrderSection
/**
 * Instantiates a new second _order section, with values for the numerator and denominator coefficients.

 * @param b0         Numerator coefficient b[0].
 * *
 * @param b1         Numerator coefficient b[1].
 * *
 * @param b2         Numerator coefficient b[2].
 * *
 * @param a1         Denominator coefficient a[1].
 * *
 * @param a2         Denominator coefficient a[2].
 */
(
        /** Numerator coefficients  */
        internal var b0: Double, internal var b1: Double, internal var b2: Double,
        /** Denominator coefficients (a0 = 1) by assumption.  */
        internal var a1: Double, internal var a2: Double) {

    /** States required to support processing of a continuous data stream in consecutive, contiguous blocks.  */
    internal var s1: Double = 0.toDouble()
    internal var s2: Double = 0.toDouble()


    init {

        initialize()
    }


    /**
     * Initializes states to zero.
     */
    fun initialize() {
        s1 = 0.0
        s2 = 0.0
    }


    /**
     * Filters a single input sample (single-step filtering).

     * @param x     float containing value of the single input sample.
     * *
     * @return      float result of the filter for one time step.
     */
    @JsName("filterSingle")
    fun filter(x: Float): Float {

        val s0 = x.toDouble() - a1 * s1 - a2 * s2
        val retval = (b0 * s0 + b1 * s1 + b2 * s2).toFloat()

        s2 = s1
        s1 = s0

        return retval
    }


    /**
     * Filters a sequence of input samples.

     * @param x     float[] containing the sequence of input samples.
     * *
     * @param y     float[] containing the filtered result.  May be the same array as x.
     */
    @JsName("filter")
    fun filter(x: FloatArray, y: FloatArray) {

        var s0: Double

        val n = min(x.size, y.size)

        for (i in 0..n - 1) {
            s0 = x[i].toDouble() - a1 * s1 - a2 * s2
            y[i] = (b0 * s0 + b1 * s1 + b2 * s2).toFloat()
            s2 = s1
            s1 = s0
        }
    }


    /**
     * Prints the filter coefficients and states.

     * @param ps       PrintStream to which the filter coefficients and states are printed.
     */
    override fun toString():String {
        var out = "  coefficients: \n"
        out += "    b0: " + b0+'\n'
        out += "    b1: " + b1+'\n'
        out += "    b2: " + b2+'\n'
        out += '\n'
        out += "    a1: " + a1+'\n'
        out += "    a2: " + a2+'\n'
        out += "\n  states:  \n"
        out += "    s1: " + s1+'\n'
        out += "    s2: " + s2+'\n'
        return out
    }

}
