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
 * Class to design analog Chebyshev type I prototype filters.

 * This analog prototype is a lowpass filter with a cutoff at 1 radian per second.  Chebyshev
 * type I filters have ripples in the passband and a sharper transition from passband to
 * stopband than a Butterworth filter.

 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
class AnalogChebyshevI
/**
 * Instantiates a new analog Chebyshev type I prototype filter.

 * @param order    int specifying the number of poles of the filter.
 * *
 * @param epsilon  double parameter controlling the stopband attenuation and passband ripple.
 */
(order: Int, epsilon: Double) : AnalogPrototype() {

    init {

        val alpha = (1.0 + Math.sqrt(1.0 + epsilon * epsilon)) / epsilon
        val p = Math.pow(alpha, 1.0 / order)
        val a = 0.5 * (p - 1 / p)
        val b = 0.5 * (p + 1 / p)

        println("alpha: " + alpha)
        println("p:     " + p)
        println("a:     " + a)
        println("b:     " + b)

        val nRealPoles = order - 2 * (order / 2)
        val nComplexPolePairs = order / 2
        val nPoles = nRealPoles + 2 * nComplexPolePairs

        if (nRealPoles == 1) {
            val td = doubleArrayOf(a, 1.0)
            addSection(Rational(Polynomial(1.0), Polynomial(td)))
        }

        val dAngle = Math.PI / nPoles

        for (i in 0..nComplexPolePairs - 1) {
            val angle = -Math.PI / 2 + dAngle / 2 * (1 + nRealPoles) + i * dAngle
            val pole = Complex(a * Math.sin(angle), b * Math.cos(angle))
            val td = doubleArrayOf(pole.real() * pole.real() + pole.imag() * pole.imag(), -2 * pole.real(), 1.0)
            addSection(Rational(Polynomial(1.0), Polynomial(td)))
        }

        // scale to 1 at s = 0

        sections[0].timesEquals(1.0 / (Math.pow(2.0, (order - 1).toDouble()) * epsilon))

    }

}
