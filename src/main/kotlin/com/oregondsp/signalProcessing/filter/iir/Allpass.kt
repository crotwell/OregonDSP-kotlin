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

import java.io.PrintStream
import java.text.DecimalFormat
import java.util.Arrays

import com.oregondsp.signalProcessing.filter.Polynomial
import com.oregondsp.signalProcessing.filter.Rational


/**
 * Allpass filter class.

 *
 * Implements a digital allpass filter, using a lattice structure implementation.  Allpass
 * filters have unit gain by construction, and are used to alter the phase of signals with
 * altering their magnitude responses.

 *
 * An allpass filter has the form:  z<sup>-N</sup> A<sub>N</sub>(z<sup>-1</sup>) / A<sub>N</sub>(z), where N is the order of the
 * filter and A<sub>N</sub>(z) is a polynomial.  The numerator polynomial is a "reflection" of the
 * denominator polynomial.  If

 *
 * A<sub>N</sub>(z)  =  a<sub>0</sub> + a<sub>1</sub>*z<sup>-1</sup> + a<sub>2</sub>*z<sup>-2</sup> + ... + a<sub>N</sub>*z<sup>-N</sup>


 *
 * then the numerator has the following form:

 *
 * z<sup>-N</sup> A<sub>N</sub>(z<sup>-1</sup>)  =  a<sub>N</sub> + a<sub>N-1</sub>*z<sup>-1</sup>) + ... + a<sub>0</sub>*z<sup>-N</sup>

 *
 * The allpass filter is represented internally as a set of reflection coefficients.  These
 * are numbers k in the open interval (-1, 1) and define A<sub>N</sub>(z) by the recursion:

 *
 * A<sub>i</sub>(z) = A<sub>i-1</sub>(z) + z<sup>-i</sup>*k(i)*A<sub>i-1</sub>(z<sup>-1</sup>)     i = 1, ..., N

 *
 * A<sub>0</sub>(z) = 1


 * @author David B. Harris, Deschutes Signal Processing  LLC
 */
open class Allpass {

    /** double[] containing the reflection coefficients specifying this allpass filter.  */
    protected var k: DoubleArray

    /** int containing the order of the filter.  */
    protected var order: Int = 0

    /** double[] containing the state of the filter.  Used to assure continuity in the filtering operation
     * when making repeated calls to the filter methods.  */
    protected var state: DoubleArray

    /** A Rational object containing the rational response representation of the filter.  */
    protected var T: Rational


    /**
     * Instantiates a new allpass filter of a given order with zero reflection coefficients.

     * @param order     int containing the order of the filter.
     */
    constructor(order: Int) {
        this.order = order
        k = DoubleArray(order)
        state = DoubleArray(order + 1)

        constructRationalRepresentation()
    }


    /**
     * Instantiates a new allpass filter object from a Polynomial representing the allpass denominator.

     * @param A     Polynomial object containing the polynomial coefficient representation for the allpass filter.
     */
    constructor(A: Polynomial) {
        k = A.reflectionCoefficients()
        order = k.size
        state = DoubleArray(order + 1)

        constructRationalRepresentation()
    }


    /**
     * Instantiates a new allpass filter object from a double[] containing the reflection coefficients.

     * @param k    double[] containing the reflection coefficients.
     */
    constructor(k: DoubleArray) {
        this.k = k.clone()
        order = this.k.size
        state = DoubleArray(order + 1)

        constructRationalRepresentation()
    }


    /**
     * Initializes the states of the filter to zero.
     */
    fun initialize() {
        Arrays.fill(state, 0.0)
    }


    /**
     * Filters a single sample of a sequence.

     * Because the filter maintains state between calls to filter(x), a sequence (signal) can be filtered
     * by repeated calls to this method using successive continguous samples of the sequence as the argument.

     * @param x      float containing the single sample of a sequence.
     * *
     * @return       float containing the corresponding single output sample of the filter.
     */
    fun filter(x: Float): Float {
        var x = x

        var stage = order
        while (stage >= 0) {

            if (stage > 0) {
                x -= (k[stage - 1] * state[stage - 1]).toFloat()
                state[stage] = k[stage - 1] * x + state[stage - 1]
            } else {
                state[stage] = x.toDouble()
            }

            stage--
        }

        return state[order].toFloat()
    }


    /**
     * Filters a sequence or a segment of a sequence contained in a float array.

     * Because the filter maintains state between successive calls to this method, it is possible
     * to filter a continuous stream in finite consecutive, contiguous blocks with this method with
     * no startup transients between successive calls.  This feature permits filtering in real time
     * as blocks of data become available, or filtering very long sequences incrementally from large
     * files.

     * @param x      float[] containing the sequence or segment of a sequence upon call.  Contains
     * *               the filtered array following the call.
     */
    fun filter(x: FloatArray) {

        for (i in x.indices)
            x[i] = filter(x[i])

    }


    /**
     * Evaluate the filters response at digital frequency omega, element of [0, pi].

     * Calculates A_N( e^(-j*omega) )

     * @param omega     double specifying the discrete frequency for evaluation.
     * *
     * @return          Complex object containing the value of the filter transfer function at omega.
     */
    fun evaluate(omega: Double): Complex {
        val ejOmega = Complex.exp(Complex(0.0, -omega))
        return T.evaluate(ejOmega)
    }


    /**
     * Evaluates the group delay of the allpass filter at digital frequency Omega, element of [0, pi].

     * @param Omega     double containing the discrete frequency for evaluation of the group delay.
     * *
     * @return          double containing the group delay in samples.
     */
    fun groupDelay(Omega: Double): Double {
        return T.discreteTimeGroupDelay(Omega)
    }


    /**
     * Constructs a Rational object with the transfer function representation of this allpass filter.

     * @return       Rational object containing the transfer function representation.
     */
    protected fun constructRationalRepresentation() {

        val a = DoubleArray(order + 1)
        val b = DoubleArray(order + 1)

        a[0] = 1.0
        for (p in 0..order - 1) {

            Arrays.fill(b, 0.0)

            var i = 0
            while (i <= p) {
                b[i] += a[i]
                b[i + 1] += k[p] * a[p - i]
                i++
            }

            System.arraycopy(b, 0, a, 0, p + 2)
        }

        Arrays.fill(b, 0.0)
        for (i in 0..order) b[i] = a[order - i]

        T = Rational(Polynomial(b), Polynomial(a))
    }


    /**
     * Accessor for the Rational representation of this allpass filter.

     * @return      Rational object containing the rational transfer function for this filter.
     */
    fun rationalRepresentation(): Rational {
        return Rational(T)
    }


    /**
     * Prints the reflection coefficients defining this allpass filter and the state vector.

     * @param ps      PrintStream to which the allpass filter description is printed.
     */
    fun print(ps: PrintStream) {

        val DF = DecimalFormat("0.000000")

        ps.println("Allpass order:  " + order)
        for (i in 0..order - 1) {
            if (i < order) {
                if (k[i] < 0.0)
                    ps.println("  " + DF.format(k[i]) + "  " + state[i])
                else
                    ps.println("   " + DF.format(k[i]) + "  " + state[i])
            } else
                ps.println("             " + state[i])
        }

    }

}
