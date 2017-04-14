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

import kotlin.js.Math


/**
 * Class to represent complex numbers and certain basic functions of complex numbers.
 * @author David B. Harris    Deschutes Signal Processing  LLC
 */
class Complex {

    /** Real part of the complex number.  */
    private var real: Double = 0.toDouble()

    /** Imaginary part of the complex number.  */
    private var imag: Double = 0.toDouble()

    // constructors

    /**
     * Instantiates a new complex number object.

     * @param real           double specifying the real part.
     * *
     * @param imag           double specifying the imaginary part.
     */
    constructor(real: Double, imag: Double) {
        this.real = real
        this.imag = imag
    }


    /**
     * Instantiates a new complex number object from a real number (imaginary part is zero).

     * @param real       double specifying the real part.
     */
    constructor(real: Double) {
        this.real = real
        this.imag = 0.0
    }


    // other methods


    /**
     * Returns the real part of a complex number.

     * @return        double containing the real part of a Complex number object.
     */
    fun real(): Double {
        return real
    }


    /**
     * Returns the imaginary part of a complex number.

     * @return        double containing the imaginary part of a Complex number object.
     */
    fun imag(): Double {
        return imag
    }


    /**
     * Computes the absolute value of this Complex number.

     * @return     double containing the absolute value of this Complex number.
     */
    fun abs(): Double {
        return abs(this)
    }


    /**
     * Computes the phase angle of this Complex number.

     * @return     double containing the phase angle of this Complex number.
     */
    fun angle(): Double {
        return angle(this)
    }


    /**
     * Multiplies this Complex number by another Complex number.

     * Does not alter the value of this Complex object.

     * @param c     The other Complex factor.
     * *
     * @return      New Complex object containing the product.
     */
    @JsName("timesComplex")
    operator fun times(c: Complex): Complex {
        return multiply(this, c)
    }


    /**
     * Multiplies this Complex number by a real number.

     * Does not alter the value of this Complex object.

     * @param a     The real multiplicand.
     * *
     * @return      New Complex object containing the product.
     */
    @JsName("timesReal")
    operator fun times(a: Double): Complex {
        return multiply(this, a)
    }


    /**
     * Returns the conjugate of this Complex number.

     * Does not alter the value of this Complex object.

     * @return     New Complex object containing the conjugate of this.
     */
    fun conjugate(): Complex {
        return conjugate(this)
    }


    /**
     * Computes the sum of this Complex number and another.

     * Does not alter the value of this Complex object.

     * @param c     Complex object containing the other summand.
     * *
     * @return      New Complex object containing the sum.
     */
    @JsName("plusComplex")
    operator fun plus(c: Complex): Complex {
        return add(this, c)
    }


    /**
     * Computes the sum of this Complex number and a real number.

     * Does not alter the value of this Complex object.

     * @param a     double containing the real summand.
     * *
     * @return      New Complex object containing the sum.
     */
    @JsName("plusReal")
    operator fun plus(a: Double): Complex {
        return add(this, a)
    }


    /**
     * Subtracts a complex number from this complex number.

     * Does not alter the value of this Complex object.

     * @param c     Complex number to be subtracted from this Complex number.
     * *
     * @return      New Complex object containing the difference.
     */
    @JsName("minusComplex")
    operator fun minus(c: Complex): Complex {
        return subtract(this, c)
    }


    /**
     * Subtracts a real number from this Complex number.

     * Does not alter the value of this Complex object.

     * @param a    double containing the real number to be subtracted from this Complex number.
     * *
     * @return     New Complex object containing the complex difference.
     */
    @JsName("minusReal")
    operator fun minus(a: Double): Complex {
        return subtract(this, a)
    }


    /**
     * Divides this Complex number by a real number.

     * Does not alter the value of this Complex object.

     * @param a    double containing the real divisor.
     * *
     * @return     New Complex object containing the result of division.
     */
    @JsName("overReal")
    fun over(a: Double): Complex {
        return divide(this, a)
    }


    /**
     * Divides this Complex number by another Complex number.

     * Does not alter the value of this Complex object.

     * @param c     The Complex divisor.
     * *
     * @return      New Complex object containint the result of division.
     */
    @JsName("overComplex")
    fun over(c: Complex): Complex {
        return divide(this, c)
    }


    /**
     * Adds a real number to this Complex number.

     * Alters the value of this Complex object.

     * @param a   double containing the other summand.
     */
    @JsName("plusEqualsReal")
    fun plusEquals(a: Double) {
        real += a
    }


    /**
     * Adds a Complex number to this Complex number.

     * Alters the value of this Complex object.

     * @param c   Complex object containing the other summand.
     */
    @JsName("plusEqualsComplex")
    fun plusEquals(c: Complex) {
        real += c.real
        imag += c.imag
    }


    /**
     * Subtracts a real number from this Complex number.

     * Alters the value of this Complex object.

     * @param a    double containing the real number to be subtracted from this Complex number.
     */
    @JsName("minusEqualsReal")
    fun minusEquals(a: Double) {
        real -= a
    }


    /**
     * Subtracts another Complex number from this Complex number.

     * Alters the value of this Complex object.

     * @param c    The other Complex number to be subtracted from this Complex number.
     */
    @JsName("minusEqualsComplex")
    fun minusEquals(c: Complex) {
        real -= c.real
        imag -= c.imag
    }


    /**
     * Multiplies this Complex number by a real number.

     * Alters the value of this Complex object.

     * @param a     double containing the multiplicand.
     */
    @JsName("timesEqualsReal")
    fun timesEquals(a: Double) {
        real *= a
        imag *= a
    }


    /**
     * Multiplies this Complex number by another Complex number.

     * Alters the value of this Complex object.

     * @param c     Complex object containing the other multiplicand.
     */
    @JsName("timesEqualsComplex")
    fun timesEquals(c: Complex) {
        val tmp = real * c.real - imag * c.imag
        imag = real * c.imag + imag * c.real
        real = tmp
    }


    /**
     * Divides this Complex number by a real number.

     * Alters the value of this Complex object.

     * @param a     double containing the real divisor.
     */
    @JsName("divideEqualsReal")
    fun divideEquals(a: Double) {
        real /= a
        imag /= a
    }


    /**
     * Divides this Complex number by another Complex number.

     * Alters the value of this Complex object.

     * @param c     Complex object containing the divisor.
     */
    @JsName("divideEqualsComplex")
    fun divideEquals(c: Complex) {
        val scale = c.real * c.real + c.imag * c.imag
        val tmp = c.real * real + c.imag * imag
        imag = c.real * imag - c.imag * real
        real = tmp
        this.divideEquals(scale)
    }


    /**
     * Generates a String representation for this Complex number object.
     */
    override fun toString(): String {
        return ""+real + "  +  i * " + imag+"\n"
    }

    companion object {


        // static methods


        /**
         * Instantiates a new complex number object from polar representation parameters.

         * @param r          double specifying the radius (magnitude) of the complex number.
         * *
         * @param phi        double specifying the phase angle of the complex number.
         * *
         * @return           Resulting Complex number object.
         */
        fun ComplexFromPolar(r: Double, phi: Double): Complex {
            return Complex(r * Math.cos(phi), r * Math.sin(phi))
        }


        /**
         * Calculates the sum of a real number and a complex number.

         * @param a         double specifying the real number.
         * *
         * @param c         Complex number object.
         * *
         * @return          New Complex object containing the sum.
         */
        fun add(a: Double, c: Complex): Complex {
            return Complex(a + c.real, c.imag)
        }


        /**
         * Calculates the sum of a complex number and a real number.

         * @param c         Complex number object.
         * *
         * @param a         double specifying the real number.
         * *
         * @return          New Complex object containing the sum.
         */
        fun add(c: Complex, a: Double): Complex {
            return add(a, c)
        }


        /**
         * Calculates the difference of a complex number and a real number.

         * @param c         Complex number object.
         * *
         * @param a         double specifying the real number.
         * *
         * @return          New Complex object containing the difference.
         */
        fun subtract(c: Complex, a: Double): Complex {
            return Complex(c.real - a, c.imag)
        }


        /**
         * Calculates the difference of a real number and a complex number.

         * @param a         double specifying the real number.
         * *
         * @param c         Complex number object.
         * *
         * @return          New Complex object containing the difference.
         */
        fun subtract(a: Double, c: Complex): Complex {
            return Complex(a - c.real, c.imag)
        }


        /**
         * Unary minus - negates a complex number.

         * @param c         Complex number to be negated.
         * *
         * @return          New Complex object containing the negative of the operand.
         */
        fun unaryMinus(c: Complex): Complex {
            return Complex(-c.real, -c.imag)
        }


        /**
         * Multiplies a real and a complex number.

         * @param a           double specifying the real factor.
         * *
         * @param c           Complex object specifying the complex factor.
         * *
         * @return            New Complex object containing the product.
         */
        fun multiply(a: Double, c: Complex): Complex {
            return Complex(a * c.real, a * c.imag)
        }


        /**
         * Multiplies a real and a complex number.

         * @param c           Complex object specifying the complex factor.
         * *
         * @param a           double specifying the real factor.
         * *
         * @return            New Complex object containing the product.
         */
        fun multiply(c: Complex, a: Double): Complex {
            return multiply(a, c)
        }


        /**
         * Adds two complex numbers.

         * @param c1     First Complex summand.
         * *
         * @param c2     Second Complex summand.
         * *
         * @return       New Complex object containing the sum.
         */
        fun add(c1: Complex, c2: Complex): Complex {
            return Complex(c1.real + c2.real, c1.imag + c2.imag)
        }


        /**
         * Subtracts one complex number from another.

         * @param c1      First Complex number.
         * *
         * @param c2      Second Complex number to be subtracted from the first.
         * *
         * @return        New Complex object containing the difference.
         */
        fun subtract(c1: Complex, c2: Complex): Complex {
            return Complex(c1.real - c2.real, c1.imag - c2.imag)
        }


        /**
         * Multiplies two complex numbers.

         * @param c1      First Complex factor.
         * *
         * @param c2      Second Complex factor.
         * *
         * @return        New Complex object containing the product.
         */
        fun multiply(c1: Complex, c2: Complex): Complex {
            return Complex(c1.real * c2.real - c1.imag * c2.imag, c1.real * c2.imag + c1.imag * c2.real)
        }


        /**
         * Divides a complex number by a real number.

         * @param c     The Complex number.
         * *
         * @param a     double containing the real divisor.
         * *
         * @return      New Complex object containing the result of division.
         */
        fun divide(c: Complex, a: Double): Complex {
            return Complex(c.real / a, c.imag / a)
        }


        /**
         * Divide a real number by a complex number.

         * @param a     double containing the real number.
         * *
         * @param c     Complex divisor.
         * *
         * @return      New Complex object containing the result of division.
         */
        fun divide(a: Double, c: Complex): Complex {
            val scale = c.real * c.real + c.imag * c.imag
            return Complex(a*c.real / scale, a*(-c.imag) / scale)
        }


        /**
         * Divides one complex number by another.

         * @param c1       The first Complex number.
         * *
         * @param c2       The Complex divisor.
         * *
         * @return         New Complex object containing the result of division.
         */
        fun divide(c1: Complex, c2: Complex): Complex {   // c1/c2 = conjg(c2)*c1/( conjg(c2)*c2 )
            val scale = c2.real * c2.real + c2.imag * c2.imag
            return Complex((c1.real * c2.real + c1.imag * c2.imag) / scale, (c1.imag * c2.real - c1.real * c2.imag) / scale)
        }


        /**
         * Computes the square root of a complex number.

         * @param c     Complex argument of the square root function.
         * *
         * @return      New Complex object containing the square root of the argument.
         */
        fun sqrt(c: Complex): Complex {
            return ComplexFromPolar(Math.sqrt(abs(c)), angle(c) / 2.0)
        }


        /**
         * Computes the absolute value of a complex number.

         * @param c     Complex argument of the absolute value operator.
         * *
         * @return      double containing the absolute value of the argument.
         */
        fun abs(c: Complex): Double {
            return Math.sqrt(c.real * c.real + c.imag * c.imag)
        }


        /**
         * Computes the phase angle of a complex number.

         * @param c      Complex argument of the phase function.
         * *
         * @return       double containing the phase of the argument.
         */
        fun angle(c: Complex): Double {
            return Math.atan2(c.imag, c.real)
        }


        /**
         * Computes the complex exponential function of a complex number.

         * @param c      Complex argument to the exponential function.
         * *
         * @return       New Complex object containing the complex exponential of the argument.
         */
        fun exp(c: Complex): Complex {
            val r = Math.exp(c.real)
            return Complex(r * Math.cos(c.imag), r * Math.sin(c.imag))
        }


        /**
         * Conjugates a complex number.

         * @param c       Complex argument.
         * *
         * @return        New Complex object containing the conjugate of the argument.
         */
        fun conjugate(c: Complex): Complex {
            return Complex(c.real, -c.imag)
        }
    }


}
