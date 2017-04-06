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

import java.io.PrintStream

import com.oregondsp.signalProcessing.filter.iir.Complex


/**
 * Class implementing rational functions.  Used to design and represent analog and digital filters.

 * The rational function is of the form H(s) = N(s)/D(s), where N and D are polynomials.

 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
class Rational {

    /** Numerator polynomial.  */
    private val N: Polynomial

    /** Denominator polynomial.  */
    private val D: Polynomial


    /**
     * Instantiates a new rational function, given arrays specifying the numerator and denominator polynomials.

     * @param num      double[] specifying coefficients of the numerator polynomial.
     * *
     * @param denom    double[] specifying coefficients of the denominator polynomial.
     */
    constructor(num: DoubleArray, denom: DoubleArray) {
        N = Polynomial(num)
        D = Polynomial(denom)
    }


    /**
     * Instantiates a new rational function, given Polynomial objects for the numerator and denominator.

     * @param N        Polynomial instance specifying the numerator polynomial.
     * *
     * @param D        Polynomial instance specifying the denominator polynomial.
     */
    constructor(N: Polynomial, D: Polynomial) {
        this.N = Polynomial(N)
        this.D = Polynomial(D)
    }


    /**
     * Constructs a new rational function by copying an existing rational function object.

     * @param R       Rational function object to be copied.
     */
    constructor(R: Rational) {
        this.N = Polynomial(R.N)
        this.D = Polynomial(R.D)
    }


    /**
     * Instantiates a new rational function object from a constant.

     * @param c      double specifying the constant of the numerator.
     */
    constructor(c: Double) {
        N = Polynomial(c)
        D = Polynomial(1.0)
    }


    /**
     * Returns the orders of the numerator and denominator polynomials.

     * @return  int[] - first element is the order of the numerator and second element is the order of the denominator.
     */
    fun order(): IntArray {
        val retval = intArrayOf(N.order(), D.order())
        return retval
    }


    /**
     * Returns a copy of the numerator polynomial.

     * @return    Polynomial copy of the numerator.
     */
    fun numerator(): Polynomial {
        return Polynomial(N)
    }


    /**
     * Returns a copy of the denominator polynomial.

     * @return    Polynomial copy of the denominator.
     */
    fun denominator(): Polynomial {
        return Polynomial(D)
    }


    /**
     * Puts the rational function representation in canonical form.

     * Normalizes the numerator and denominator polynomials to have unit lead coefficients.  Returns the
     * gain factor required to perform the normalization.

     * @return   double specifying the gain
     */
    fun canonicalForm(): Double {

        val scaleN = N.a[N.order]
        for (i in N.a.indices) N.a[i] /= scaleN
        val scaleD = D.a[D.order]
        for (i in D.a.indices) D.a[i] /= scaleD

        return scaleN / scaleD
    }


    /**
     * Scales (in-place) a rational function by a constant.


     * @param A    double specifying the scale factor.
     */
    fun timesEquals(A: Double) {
        N.timesEquals(A)
    }


    /**
     * Multiplies (in-place) a rational function by a polynomial.

     * @param P    Polynomial object specifying the multiplicative polynomial factor.
     */
    fun timesEquals(P: Polynomial) {
        N.timesEquals(P)
    }


    /**
     * Multiplies(in-place) a rational function by another rational function.

     * @param R    Rational object specifying the multiplicative rational factor.
     */
    fun timesEquals(R: Rational) {
        N.timesEquals(R.N)
        D.timesEquals(R.D)
    }


    /**
     * Evaluates the rational function for a real argument.

     * @param x     double specifying the argument to the rational function.
     * *
     * @return      double specifying the resulting value of the rational function.
     */
    fun evaluate(x: Double): Double {
        var retval = 0.0
        val num = N.evaluate(x)
        val denom = D.evaluate(x)
        if (denom != 0.0) retval = num / denom

        return retval
    }


    /**
     * Evaluates a rational function for a complex argument.

     * @param c    Complex object specifying the complex argument.
     * *
     * @return     Complex object specifying the resulting complex value of the rational function.
     */
    fun evaluate(c: Complex): Complex {
        var retval = Complex(0.0, 0.0)
        val num = N.evaluate(c)
        val denom = D.evaluate(c)
        if (denom.abs() != 0.0) retval = num.over(denom)

        return retval
    }


    /**
     * Maps a rational function to a new rational function by substitution of a rational function for the argument.

     * Uses a modification of Horner's scheme to perform the mapping.

     * @param S     Rational object specifying the mapping function.
     * *
     * @return      Rational object specifying the resulting mapped rational function.
     */
    fun map(S: Rational): Rational {

        //  Modified Horner's scheme evaluation

        //    Numerator

        var P = Polynomial(N.a[N.order])
        var T = Polynomial(1.0)
        for (i in N.order - 1 downTo 0) {
            T = T.times(S.D)
            P = P.times(S.N).plus(T.times(N.a[i]))
        }

        //    Denominator

        var Q = Polynomial(D.a[D.order])
        T = Polynomial(1.0)
        for (i in D.order - 1 downTo 0) {
            T = T.times(S.D)
            Q = Q.times(S.N).plus(T.times(D.a[i]))
        }

        if (D.order > N.order) {
            for (i in 0..D.order - N.order - 1)
                P = P.times(S.D)
        } else if (N.order > D.order) {
            for (i in 0..N.order - D.order - 1)
                Q = Q.times(S.D)
        }

        P.trim()
        Q.trim()

        return Rational(P, Q)
    }


    /**
     * Calculates the residue of a real pole of the rational function.

     * Uses L'Hopital's rule to calculate the residue of a real pole.  Potentially useful to find parallel implementations
     * of IIR digital filters.  Suitable only for simple poles.

     * @param pole     double specifying the real pole (root of the denominator).
     * *
     * @return         double specifying the residue.
     */
    fun residue(pole: Double): Double {

        // using L'Hopital's rule - assumes single pole

        return N.evaluate(pole) / D.derivative().evaluate(pole)
    }


    /**
     * Calculates the residue of a complex pole of the rational function.

     * Uses L'Hopital's rule to calculate the residue of a complex pole.  Potentially useful to find parallel implementations
     * of IIR digital filters.  Suitable only for simple poles.

     * @param pole     Complex object specifying the pole.
     * *
     * @return         Complex object specifying the residue.
     */
    fun residue(pole: Complex): Complex {

        // using L'Hopital's rule - assumes single pole

        return N.evaluate(pole).over(D.derivative().evaluate(pole))
    }


    /**
     * Calculates the group delay of an analog filter with transfer function specified by this rational function.

     * @param omega      double specifying the radial frequency (2*pi*f) at which the group delay is evaluated.
     * *
     * @return           double specifying the group delay in seconds.
     */
    fun groupDelay(omega: Double): Double {
        return N.groupDelay(omega) - D.groupDelay(omega)
    }


    /**
     * Calculates the group delay of a digital filter with transfer function specified by this rational function.

     * For this evaluation, the numerator and denominator are assumed to be polynomials in z^-1.

     * @param Omega    double specifying the digital frequency [0 pi] at which the group delay is evaluated.
     * *
     * @return         double specifying the group delay in samples.
     */
    fun discreteTimeGroupDelay(Omega: Double): Double {
        return N.discreteTimeGroupDelay(Omega) - D.discreteTimeGroupDelay(Omega)
    }


    /**
     * Prints the coefficients of the rational function.

     * @param ps      Printstream to which the rational function coefficients are printed.
     */
    fun print(ps: PrintStream) {
        ps.println("Numerator: ")
        N.print(ps)
        ps.println("Denominator: ")
        D.print(ps)
    }

    companion object {


        @JvmStatic fun main(args: Array<String>) {
            val a = DoubleArray(4)
            a[0] = 1.0
            a[1] = 2.0
            a[2] = 2.0
            a[3] = 1.0
            val b = DoubleArray(1)
            b[0] = 1.0
            val R = Rational(b, a)

            for (i in 0..99) {
                val omega = i / 25.0
                println(omega.toString() + "  " + R.evaluate(Complex(0.0, omega)).abs() + "   " + R.groupDelay(omega))
            }
        }
    }


}
