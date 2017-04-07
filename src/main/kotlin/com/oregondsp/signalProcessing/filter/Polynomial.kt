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

package com.oregondsp.signalProcessing.filter


import com.oregondsp.signalProcessing.filter.iir.Complex
import kotlin.js.Math


/**
 * Class to implement polynomial functions

 *
 * This class represents polynomials by their coefficients, stored in a double array, with the
 * lowest coefficient (constant) stored at the 0 location.  Consequently the polynomial is of the
 * form:

 *
 * A(x) = a[0] + a[1]*x + a[2]*x<sup>2</sup> + ... + a[N]*x<sup>N</sup>

 *
 * The class implements basic polynomial arithmetic and other functions useful in designing
 * analog and digital filters.  An example of the latter includes the calculation of reflection
 * coefficients for use in allpass filters.

 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
class Polynomial {

    /** double array containing the coefficients of the polynomial.  The low order coefficient is
     * in a[0] and the high order coefficient is in a[order].
     */
    var a: DoubleArray

    /** integer containing the order (degree:  N) of the polynomial.  */
    var order: Int = 0


    /**
     * Instantiates a new polynomial from a double array containing the coefficients.

     * @param a    double[] containing the polynomial coefficients.
     */
    constructor(a: DoubleArray) {
        order = a.size - 1
        this.a = a.copyOf()
    }


    /**
     * Instantiates (copies) a new polynomial from an existing Polynomial instance.

     * @param B       Polynomial to be copied into the new Instance.
     */
    constructor(B: Polynomial) {
        order = B.order
        this.a = B.a.copyOf()
    }


    /**
     * Instantiates a new zero polynomial.

     * @param order     int containing the order (degree) of the polynomial.
     */
    constructor(order: Int) {
        this.order = order
        this.a = DoubleArray(order + 1)
        //Arrays.fill(a, 0.0)
    }


    /**
     * Instantiates a new constant polynomial.

     * @param c        double containing the constant.
     */
    constructor(c: Double) {
        order = 0
        a = DoubleArray(1)
        a[0] = c
    }


    /**
     * Removes leading zero coefficients in a polynomial.

     * This method is used by the Rational class when new polynomials have been constructed in a rational
     * mapping operation.
     */
    fun trim() {

        var i = order
        var n = 0
        while (a[i] == 0.0) {
            n++
            i--
        }

        if (n > 0) {
            val b = DoubleArray(order + 1 - n)
            for (j in b.indices) {
                b[j] = a[j]
            }
            a = b
            order -= n
        }

    }


    /**
     * Returns the order (degree) of the polynomal.

     * @return      int containing the degree of the polynomial.
     */
    fun order(): Int {
        return order
    }


    /**
     * Returns the polynomial coefficients.

     * @return      double[] containing a copy of the coefficients.
     */
    fun coefficients(): DoubleArray {
        return a.copyOf()
    }


    /**
     * Returns a new Polynomial object containing the sum polynomial of this and a constant.

     * The original Polynomial is unchanged.

     * @param c        double constant to be added to this polynomial.
     * *
     * @return         The new polynomial containing the sum of this polynomial and the constant.
     */
    operator fun plus(c: Double): Polynomial {
        val retval = Polynomial(order)
        retval.a = a.copyOf()
        retval.a[0] += c
        return retval
    }


    /**
     * Adds a constant to this polynomial.  Alters this to contain the sum.

     * @param c       double constant to be added to this polynomial.
     */
    fun plusEquals(c: Double) {
        a[0] += c
    }


    /**
     * Returns a new Polynomial object containing the sum of this and the argument polynomial.

     * The original Polynomial is unchanged.

     * @param B       Polynomial to be added to this polynomial.
     * *
     * @return        Polynomial object containing the new sum.
     */
    operator fun plus(B: Polynomial): Polynomial {
        val retval = Polynomial(Math.max(order, B.order))
        for (i in 0..order) retval.a[i] = a[i]
        for (i in 0..B.order) retval.a[i] += B.a[i]
        return retval
    }


    /**
     * Adds a polynomial to this polynomial.  Alters this to contain the sum.

     * @param B       Polynomial object containing the polynomial to be added to this polynomial.
     */
    fun plusEquals(B: Polynomial) {
        val A = DoubleArray(Math.max(order, B.order))
        for (i in 0..order) A[i] = a[i]
        for (i in 0..B.order) A[i] += B.a[i]
        a = A
        order = A.size - 1
    }


    /**
     * Returns a new Polynomial containing the difference between this polynomial and a constant.

     * This polynomial is unchanged.

     * @param c        double containing the constant to be subtracted from this polynomial.
     * *
     * @return         The new Polynomial object containing the difference polynomial.
     */
    operator fun minus(c: Double): Polynomial {
        return plus(-c)
    }


    /**
     * Subtracts a constant from this polynomial.  Alters this to contain the difference.

     * @param c        double constant to be subtracted from this polynomial.
     */
    fun minusEquals(c: Double) {
        plusEquals(-c)
    }


    /**
     * Subtracts a polynomial from this polynomial.  Returns the difference as a new Polynomial.

     * This polynomial is unchanged.

     * @param B            Polynomial to be subtracted from this Polynomial.
     * *
     * @return             Polynomial containing the difference.
     */
    operator fun minus(B: Polynomial): Polynomial {
        val retval = Polynomial(Math.max(order, B.order))
        for (i in 0..order) retval.a[i] = a[i]
        for (i in 0..B.order) retval.a[i] -= B.a[i]
        return retval
    }


    /**
     * Subtracts a polynomial from this polynomial.  Alters this to contain the difference.

     * @param B        Polynomial to be subtracted from this polynomial.
     */
    fun minusEquals(B: Polynomial) {
        val A = DoubleArray(Math.max(order, B.order))
        for (i in 0..order) A[i] = a[i]
        for (i in 0..B.order) A[i] -= B.a[i]
        a = A
        order = A.size - 1
    }


    /**
     * Computes the product of a constant and this polynomial.  The product is returned as a new Polynomial.

     * This polynomial is unchanged.

     * @param c        The double constant factor multiplying this polynomial.
     * *
     * @return         The resulting product polynomial.
     */
    operator fun times(c: Double): Polynomial {
        val retval = Polynomial(order)
        for (i in 0..order) retval.a[i] = c * a[i]
        return retval
    }


    /**
     * Multiplies (scales) this polynomial by a constant.  This polynomial is changed to contain the product.

     * @param c        The constant multiplicative factor.
     */
    fun timesEquals(c: Double) {
        for (i in 0..order) a[i] *= c
    }


    /**
     * Computes the product of this polynomial with another polynomial.  The product is returned in a new Polynomial.

     * @param B        Polynomial object containing the multiplicative factor.
     * *
     * @return         New Polynomial object containing the product.
     */
    operator fun times(B: Polynomial): Polynomial {

        val b = B.a
        val prod = DoubleArray(order + B.order + 1)
        //Arrays.fill(prod, 0.0)

        for (i in 0..B.order) {
            for (j in 0..order) {
                prod[i + j] += b[i] * a[j]
            }
        }

        return Polynomial(prod)
    }


    /**
     * Multiplies this by a Polynomial factor.  Alters this polynomial to contain the product.

     * @param B        Polynomial object containing the multiplicative factor.
     */
    fun timesEquals(B: Polynomial) {

        val b = B.a
        val prod = DoubleArray(order + B.order + 1)
        //Arrays.fill(prod, 0.0)

        for (i in 0..B.order) {
            for (j in 0..order) {
                prod[i + j] += b[i] * a[j]
            }
        }

        a = prod
        order += B.order
    }


    /**
     * Divides this polynomial by a constant and returns the result of division in a new Polynomial object.

     * This polynomial is unchanged by this method.

     * @param c        The double divisor.
     * *
     * @return         New Polynomial object containing the result of division.
     */
    fun over(c: Double): Polynomial {
        val tmp = DoubleArray(order + 1)
        for (i in 0..order + 1 - 1)
            tmp[i] = a[i] / c

        return Polynomial(tmp)
    }


    /**
     * Divides this polynomial by a constant.  This polynomial is altered to contain the result of division.

     * @param c          The double divisor.
     */
    fun overEquals(c: Double) {
        for (i in 0..order + 1 - 1)
            a[i] /= c
    }


    /**
     * Divides this polynomial by another polynomial.  Returns the result of division in a new Rational object.

     * This polynomial is unchanged by this operation.

     * @param B         The Polynomial divisor.
     * *
     * @return          New Rational object containing the result of division.
     */
    fun over(B: Polynomial): Rational {
        return Rational(this, B)
    }


    /**
     * Computes the derivative of a polynomial.  Returns the derivative as a new Polynomial object.

     * @return          New Polynomial object containing the derivative of this polynomial.
     */
    fun derivative(): Polynomial {
        val tmp = DoubleArray(order)
        for (i in 0..order - 1) {
            tmp[i] = (i + 1) * a[i + 1]
        }

        return Polynomial(tmp)
    }


    /**
     * Evaluates this polynomial for a real argument.

     * @param x        double containing the argument of the polynomial.
     * *
     * @return         double containing the value of the polynomial at this argument.
     */
    fun evaluate(x: Double): Double {

        var retval = a[order]

        for (i in order - 1 downTo 0) {
            retval = x * retval + a[i]
        }

        return retval
    }


    /**
     * Evaluates this polynomial for a complex argument.

     * @param c        Complex object containing the argument of the polynomial.
     * *
     * @return         Complex object containing the value of the polynomial at this argument.
     */
    fun evaluate(c: Complex): Complex {

        var retval = Complex(a[order])

        for (i in order - 1 downTo 0) {
            retval = retval.times(c).plus(a[i])
        }

        return retval
    }


    /**
     * Evaluates the group delay of an analog filter transfer function specified by this polynomial.

     * @param omega    double specifying the radial frequency (2*pi*f) at which the group delay is evaluated.
     * *
     * @return         double containing the resulting group delay in seconds.
     */
    fun groupDelay(omega: Double): Double {

        if (order == 0)
            return 0.0
        else {
            val c = Complex(0.0, omega)
            val N = derivative().evaluate(c)
            val D = evaluate(c)

            return -N.over(D).real()
        }

    }


    /**
     * Evaluates the group delay of a discrete time filter transfer function specified by this polynomial.

     * @param Omega         double specifying the value discrete-time frequency [0 pi] at which the group delay is evaluated.
     * *
     * @return              double containing the resulting group delay in samples.
     */
    fun discreteTimeGroupDelay(Omega: Double): Double {

        val c = Complex.exp(Complex(0.0, -Omega))

        var N = Complex(a[order] * order)
        for (i in order - 1 downTo 0) {
            N = N.times(c).plus(a[i] * i)
        }

        val D = evaluate(c)

        return N.over(D).real()
    }


    /**
     * Computes reflection coefficients for this polynomial.

     * @return         double[] containing the reflection coefficient representation for this polynomial.
     */
    fun reflectionCoefficients(): DoubleArray {

        val k = DoubleArray(order)

        // assure that polynomial is monic

        val b = DoubleArray(order + 1)
        b[0] = 1.0
        for (i in 0..order - 1) b[i + 1] = a[i + 1] / a[0]

        // recursion to calculate reflection coefficients


        for (i in order downTo 1) {

            k[i - 1] = b[i]

            val scale = 1.0 - k[i - 1] * k[i - 1]

            //Arrays.fill(c, 0.0)
            val c = DoubleArray(order)

            for (j in 0..i - 1) {
                c[j] = (b[j] - k[i - 1] * b[i - j]) / scale
            }

            for (j in 0..i - 1) {
                b[j] = c[j]
            }
        }

        return k
    }


    /**
     * Prints the coefficients of this polynomial to a PrintStream.

     * @param ps            PrintStream object to which this polynomial's coefficients are printed.
     */
    fun print(ps: PrintStream) {
        for (i in 0..order) {
            if (i >= 0 && i < 10)
                ps.println(i.toString() + "    " + a[i])
            else if (i >= 10 && i <= 100)
                ps.println(i.toString() + "   " + a[i])
        }
    }

}
