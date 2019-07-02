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

import com.oregondsp.signalProcessing.filter.Rational
import kotlin.math.*


/**
 * Class to implement an Infinite Impulse Response digital filter.

 * Implements the filter as a cascade of second-_order sections.  The filter is obtained using a
 * bilinear transformation of a prototype analog filter.  This implementation saves internal states
 * from one invocation of the filter methods to the next allowing continuous processing of real time
 * data streams or very large files in consecutive, contiguous blocks.

 * @author David B. Harris, Deschutes Signal Processing LLC
 */
open class IIRFilter
/**
 * Instantiates a new IIR filter.

 * @param baseFilter            The AnalogPrototype for this digital filter.
 * *
 * @param type                  PassbandType object specifying lowpass, highpass or bandpass type response.
 * *
 * @param f1                    double specifying the low cutoff frequency - used by highpass and bandpass types.
 * *
 * @param f2                    double specifying the high cutoff frequency - used by lowpass and bandpass types.
 * *
 * @param delta                 double specifying the sampling interval for which the filter is designed.
 */
(baseFilter: AnalogPrototype, type: PassbandType, f1: Double, f2: Double, delta: Double) {

    /** An ArrayList of second _order sections.  */
    protected var sections: ArrayList<SecondOrderSection>

    /** Rational object containing the transfer function of the filter.  */
    protected var T: Rational


    init {

        val prototype: AnalogPrototype

        when (type) {

            PassbandType.LOWPASS -> prototype = baseFilter.lptolp(warp(f2, delta))

            PassbandType.BANDPASS -> prototype = baseFilter.lptobp(warp(f1, delta), warp(f2, delta))

            PassbandType.HIGHPASS -> prototype = baseFilter.lptohp(warp(f1, delta))

            else -> throw IllegalStateException("Undefined passband type")
        }

        val tn = DoubleArray(2)
        val td = DoubleArray(2)
        tn[0] = 1.0
        tn[1] = -1.0
        td[0] = 1.0
        td[1] = 1.0
        val S = Rational(tn, td)

        T = Rational(1.0)

        sections = ArrayList<SecondOrderSection>()

        for (i in 0..prototype.nSections() - 1) {
            val R = prototype.getSection(i).map(S)
            T.timesEquals(R)
            val cn = R.numerator().coefficients()
            val cd = R.denominator().coefficients()
            var s = 1.0
            if (cd[0] != 0.0) s = cd[0]

            val b0 = cn[0] / s
            var b1 = 0.0
            if (cn.size >= 2) b1 = cn[1] / s
            var b2 = 0.0
            if (cn.size >= 3) b2 = cn[2] / s
            var a1 = 0.0
            if (cd.size >= 2) a1 = cd[1] / s
            var a2 = 0.0
            if (cd.size >= 3) a2 = cd[2] / s
            sections.add(SecondOrderSection(b0, b1, b2, a1, a2))
        }

    }


    /**
     * Initializes the states of the filter, i.e. of each of the second-_order sections.
     */
    fun initialize() {
        for (i in sections.indices) {
            sections[i].initialize()
        }
    }


    /**
     * Filters a single sample of a sequence.

     * @param x       float containing the sequence sample.
     * *
     * @return        float value of the resulting filtered sequence.
     */
    @JsName("filterNextSample")
    fun filter(x: Float): Float {
        var retval = sections[0].filter(x)
        for (i in 1..sections.size - 1)
            retval = sections[i].filter(retval)

        return retval
    }


    /**
     * Filters an array of sequence samples.

     * Suitable for use in filtering a long file or continuous data stream broken into
     * consecutive, contiguous blocks.  Maintains state between invocations, allowing
     * continuous processing.

     * @param x    float[] containing samples of the sequence to be filtered.
     * *
     * @param y    float[] containing samples of the resulting filtered sequence.
     */
    @JsName("filter")
    fun filter(x: FloatArray, y: FloatArray) {
        //Arrays.fill(y, 0.0f)
        for (i in y.indices)
            y[i] = 0.0F
        sections[0].filter(x, y)

        for (i in 1..sections.size - 1) {
            sections[i].filter(y, y)
        }
    }


    /**
     * Filters an array of sequence samples in-place.

     * In this implementation, the source and destination arrays are identical, conserving
     * storage.

     * @param x     float[] contains samples of the sequence to be filtered upon call and the filtered
     * *              samples following execution
     */
    @JsName("filterInPlace")
    fun filter(x: FloatArray) {
        for (section in sections) {
            section.filter(x, x)
        }
    }


    /**
     * Evaluates the transfer function of this IIR filter at a specified discrete time frequency.

     * @param Omega      double containing the discrete frequency (in [0, pi]) for evaluation of the transfer function.
     * *
     * @return           Complex object containing the value of the transfer function at frequency Omega.
     */
    @JsName("evaluate")
    fun evaluate(Omega: Double): Complex {
        val ejOmega = Complex.exp(Complex(0.0, -Omega))
        return T.evaluate(ejOmega)
    }


    /**
     * Computes the group delay of the IIR filter at a specified discrete time frequency.

     * @param Omega       double containing the discrete frequency (in [0, pi]) for evaluation of the group delay.
     * *
     * @return            double containing the resulting group delay.
     */
    @JsName("groupDelay")
    fun groupDelay(Omega: Double): Double {
        return T.discreteTimeGroupDelay(Omega)
    }


    /**
     * Prints the coefficients and states of this IIR filter, section by section.

     * @param ps      PrintStream object to which this filters coefficients and states are printed.
     */
    override fun toString():String {

        var out = "IIR Filter:\n"
        for (i in sections.indices) {
            out += "\n  Section " + i + "\n"
            out += sections[i]
            out += '\n'
        }
        return out
    }


    // frequency warping for bilinear transformation

    /**
     * Method to prewarp cutoff frequencies to correct for the nonlinear frequency mapping of the bilinear transformation.

     * @param f        double containing the analog frequency (typically a cutoff specification) before the bilinear transform.
     * *
     * @param delta    double specificying the sampling interval of the data.
     * *
     * @return         double containing the prewarped digital frequency correcting for the nonlinearity of the bilinear transform.
     */
    @JsName("warp")
    private fun warp(f: Double, delta: Double): Double {
        return tan(PI * f * delta)
    }


}
