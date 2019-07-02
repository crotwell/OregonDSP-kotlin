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
import kotlin.math.*


/**
 * Base class, with partial implementation, for analog prototype filters.

 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
open class AnalogPrototype {

    /** Data structure for second-_order sections comprising the filter implementation  */
    protected var sections: ArrayList<Rational>


    private var _T: Rational? = null
    /** Overall transfer function of the filter, represented as a rational object  */
    protected val T: Rational
        get() {
           if (_T == null) {
               _T = computeTransferFunction()
           }
            return _T ?: throw RuntimeException("SHould not happen, _T is null")
        }

    /**
     * Default constructor for a new analog prototype.

     * Instantiates an analog prototype with no second _order sections.  This constructor
     * is called by the super() methods in subclasses.
     */
    init {
        sections = ArrayList<Rational>()
    }


    /**
     * Method to add a second _order section to the analog prototype representation.

     * @param R     Rational object containing a second _order section representation.
     */
    @JsName("addSection")
    fun addSection(R: Rational) {
        sections.add(R)
        _T = null
    }


    /**
     * Returns the number of second _order sections in the analog prototype representation.

     * @return    int containing the number of second _order sections.
     */
    fun nSections(): Int {
        return sections.size
    }


    /**
     * Accessor for second _order sections in the prototype representation.

     * @param index     int specifying the desired second _order section.
     * *
     * @return          Rational object containing the representation of the desired section.
     */
    @JsName("getSection")
    fun getSection(index: Int): Rational {
        return Rational(sections[index])
    }


    // spectral transformations

    //   lowpass to lowpass transformation

    /**
     * Converts a lowpass prototype with cutoff at 1 rad/sec to lowpass with a new cutoff frequency.

     * @param omega0     double specifying the cutoff of the transformed lowpass prototype filter.
     * *
     * @return           AnalogPrototype object containing the transformed filter representation.
     */
    @JsName("lptolp")
    fun lptolp(omega0: Double): AnalogPrototype {

        val tn = doubleArrayOf(0.0, 1.0)
        val td = doubleArrayOf(omega0)

        val T = Rational(tn, td)

        val retval = AnalogPrototype()

        for (i in sections.indices)
            retval.addSection(sections[i].map(T))

        return retval
    }


    //   lowpass to highpass transformation

    /**
     * Converts a lowpass analog prototype with cutoff at 1 rad/sec to a highpass filter with a new cutoff.

     * @param omega0     double specifying the desired new cutoff frequency - now a low cutoff.
     * *
     * @return           AnalogPrototype object containing the transformed filter representation.
     */
    @JsName("lptohp")
    fun lptohp(omega0: Double): AnalogPrototype {

        val tn = doubleArrayOf(omega0)
        val td = doubleArrayOf(0.0, 1.0)

        val T = Rational(tn, td)

        val retval = AnalogPrototype()

        for (i in sections.indices)
            retval.addSection(sections[i].map(T))

        return retval
    }


    //   lowpass to bandpass transformation

    /**
     * Converts a lowpass analog prototype with cutoff at 1 rad/sec to a bandpass filter with specified cutoffs.

     * @param omega1     double containing the low cutoff frequency in radians/sec.
     * *
     * @param omega2     double containing the high cutoff frequency in radians/sec.
     * *
     * @return           AnalogPrototype object containing the transformed filter representation.
     */
    @JsName("lptobp")
    fun lptobp(omega1: Double, omega2: Double): AnalogPrototype {

        val BW = omega2 - omega1
        val prod = omega1 * omega2

        val tn = doubleArrayOf(prod, 0.0, 1.0)
        val td = doubleArrayOf(0.0, BW)

        val T = Rational(tn, td)

        val retval = AnalogPrototype()

        var A = 1.0

        for (i in sections.indices) {

            val section = sections[i]
            val Tsection = section.map(T)
            A *= Tsection.canonicalForm()

            val order = section.order()

            if (order[0] < 2 && order[1] < 2)
                retval.addSection(Tsection)
            else if (order[1] == 2) {

                val DT = lptobpFactors(section.denominator(), BW, prod)
                val t1 = doubleArrayOf(0.0, 1.0)

                if (order[0] == 0) {
                    retval.addSection(Rational(Polynomial(t1), DT[0]))
                    retval.addSection(Rational(Polynomial(t1), DT[1]))
                } else if (order[0] == 1) {
                    retval.addSection(Rational(Polynomial(t1), DT[0]))
                    val t2 = DoubleArray(3)
                    val tc = Tsection.numerator().coefficients()
                    for (j in 0..2) t2[j] = tc[j + 1]
                    retval.addSection(Rational(Polynomial(t2), DT[1]))
                } else if (order[0] == 2) {
                    val NT = lptobpFactors(section.numerator(), BW, prod)
                    retval.addSection(Rational(NT[0], DT[0]))
                    retval.addSection(Rational(NT[1], DT[1]))
                }

            }

        }

        retval.sections[0].timesEquals(A)

        return retval
    }


    /**
     * Method to compute polynomial factors for bandpass transformed quadratic polynomials in a second-_order section.

     * @param P       Polynomial object to be transformed.
     * *
     * @param BW      Bandwidth parameter of the transform.
     * *
     * @param prod    Product parameter of the transform.
     * *
     * @return        Array of Polynomial factors (there will be two for each quadratic in a second _order section).
     */
    @JsName("lptobpFactors")
    private fun lptobpFactors(P: Polynomial, BW: Double, prod: Double): Array<Polynomial> {

        //val retval = arrayOfNulls<Polynomial>(2)

        val p = P.coefficients()
        val c = p[0] / p[2]
        val b = p[1] / p[2]
        val discriminant = b * b - 4 * c
        var t0: Polynomial
        var t1: Polynomial

        if (discriminant >= 0.0) {
            var root = (-b + sqrt(discriminant)) / 2.0
            var f1 = root * BW / 2.0
            var f2 = f1 * f1 - prod
            var C = Complex(f1).plus(Complex.sqrt(Complex(f2)))
            t0 = Polynomial(doubleArrayOf(C.conjugate().times(C).real(), -2.0 * C.real(), 1.0))
            //retval[0] = Polynomial(t0)

            root = (-b - sqrt(discriminant)) / 2.0
            f1 = root * BW / 2.0
            f2 = f1 * f1 - prod
            C = Complex(f1).plus(Complex.sqrt(Complex(f2)))
            t1 = Polynomial(doubleArrayOf(C.conjugate().times(C).real(), -2.0 * C.real(), 1.0))
            //retval[1] = Polynomial(t1)

        } else {
            val root = Complex(-b / 2.0, sqrt(-discriminant) / 2.0)

            val f1 = root.times(BW / 2.0)
            val f2 = f1.times(f1).minus(prod)
            var C = f1.plus(Complex.sqrt(f2))
            t0 = Polynomial(doubleArrayOf(C.conjugate().times(C).real(), -2.0 * C.real(), 1.0))
            //retval[0] = Polynomial(t0)

            C = f1.minus(Complex.sqrt(f2))
            t1 = Polynomial(doubleArrayOf(C.conjugate().times(C).real(), -2.0 * C.real(), 1.0))
            //retval[1] = Polynomial(t1)
        }

        return arrayOf(Polynomial(t0), Polynomial(t1))
    }


    /**
     * Computes the transfer function representation of the filter as a product of second-_order section transfer fuctions.

     * @return     Rational object containing the resulting transfer function representation.
     */
    protected fun computeTransferFunction(): Rational {

        var T = Rational(1.0)

        for (i in sections.indices)
            T.timesEquals(sections[i])
        return T
    }


    /**
     * Accessor for the transfer function representation for the filter.

     * @return      Rational object containing the transfer function representation for the filter.
     */
    val transferFunction: Rational
        get() {
            return Rational(T)
        }


    /**
     * Evaluates the filter transfer function at analog frequency omega.

     * @param omega          double containing the analog frequency for evaluation of the transfer function.
     * *
     * @return               Complex object containing the value of the transfer function at this frequency.
     */
    @JsName("evaluate")
    fun evaluate(omega: Double): Complex {

        return T.evaluate(Complex(0.0, omega))
    }


    /**
     * Evaluates the filter's group delay at analog frequency omega.

     * @param omega          double containing the analog frequency for evaluation of the group delay.
     * *
     * @return               double containing the group delay at this frequency.
     */
    @JsName("groupDelay")
    fun groupDelay(omega: Double): Double {

        return T.groupDelay(omega)
    }


    /**
     * Prints the coefficients of the second-_order section factors of this analog prototype filter.

     * @param ps     PrintStream to which the representation is printed.
     */
    override fun toString():String {

        var out = "AnalogPrototype: \n"

        for (i in sections.indices) {
            out += "  section $i:"+'\n'
            out += sections[i]
        }
        return out
    }


}
