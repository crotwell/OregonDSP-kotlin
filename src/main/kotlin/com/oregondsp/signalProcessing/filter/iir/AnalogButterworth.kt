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
import com.oregondsp.signalProcessing.filter.Rational
import kotlin.js.Math


/**
 * Class to design analog Butterworth prototype filters.

 * This analog prototype is a lowpass filter with a cutoff at 1 radian per second.

 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
class AnalogButterworth
/**
 * Instantiates a new analog Butterworth filter design with the indicated number of poles.

 * @param order   int specifying the number of poles of the filter.
 */
(order: Int) : AnalogPrototype() {

    init {

        val nRealPoles = order - 2 * (order / 2)
        val nComplexPolePairs = order / 2
        val nPoles = nRealPoles + 2 * nComplexPolePairs

        if (nRealPoles == 1) {
            val td = doubleArrayOf(1.0, 1.0)
            addSection(Rational(Polynomial(1.0), Polynomial(td)))
        }

        val dAngle = Math.PI / nPoles

        for (i in 0..nComplexPolePairs - 1) {
            val angle = -Math.PI / 2 + dAngle / 2 * (1 + nRealPoles) + i * dAngle
            val td = doubleArrayOf(1.0, -2 * Math.sin(angle), 1.0)
            addSection(Rational(Polynomial(1.0), Polynomial(td)))
        }

    }


}
