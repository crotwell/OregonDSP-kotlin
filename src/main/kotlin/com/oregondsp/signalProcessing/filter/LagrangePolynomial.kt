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


/**
 * Implements Lagrange Polynomials.

 * This class principally supports the design of equiripple FIR digital filters with the
 * Remez exchange algorithm (package equiripple).

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
class LagrangePolynomial
/**
 * Instantiates a new Lagrange polynomial given the set of ordinate and matching abscissa values that this
 * polynomial interpolates.

 * @param x    double[] containing the ordinates of the points interpolated by this polynomial.
 * *
 * @param y    double[] containing the abscissas of the points interpolated by this polynomial.
 */
(x: DoubleArray, y: DoubleArray) {


    /** Order of the polynomial.  */
    private val order: Int

    /** The ordinates of the points interpolated by this Lagrange polynomial.  */
    private val x: DoubleArray

    /** The abscissas of the points interpolated by this Lagrange polynomial.  */
    private val y: DoubleArray

    /** Barycentric weights for this Lagrange polynomial.  */
    private val weights: DoubleArray                 // Barycentric weights

    init {

        if (x.size != y.size) throw IllegalArgumentException("Lengths of x and y arrays do not match")
        this.x = x.clone()
        this.y = y.clone()

        order = x.size - 1

        // calculate Barycentric weights

        weights = BarycentricWeights(x)

    }


    /**
     * Accessor - returns the order of the polynomial.

     * @return the int
     */
    fun order(): Int {
        return order
    }


    /**
     * Evaluates the Lagrange polynomial at real value xp.

     * Evaluates the Lagrange polynomial using the barycentric formula.

     * @param xp    double containing the real value for evaluation of the polynomial.
     * *
     * @return      double containing the value of the polynomial at xp.
     */
    fun evaluate(xp: Double): Double {

        var num = 0.0
        var denom = 0.0
        for (j in 0..order) {
            if (xp == x[j]) {
                num = y[j]
                denom = 1.0
                break
            }
            val term = weights[j] / (xp - x[j])
            num += term * y[j]
            denom += term
        }

        return num / denom
    }


    /**
     * Computes Chebyshev nodes for approximation of a function on interval [a, b].

     * The Chebyshev nodes are a particularly good set of ordinate values for
     * approximation of a function on the interval [a, b] by interpolation with
     * a Lagrange polynomial.

     * @param a    double containing the starting point of the interval.
     * *
     * @param b    double containing the ending point of the interval.
     * *
     * @param n    int specifying the number of Chebyshev nodes desired.
     * *
     * @return     double[] containing the Chebyshev nodes.
     */
    fun ChebyshevNodes(a: Double, b: Double, n: Int): DoubleArray {

        val t0 = (a + b) / 2.0
        val t1 = (b - 1) / 2.0
        val retval = DoubleArray(n)
        for (i in 0..n - 1)
            retval[i] = t0 + t1 * Math.cos((2 * i + 1) / (2 * n) * Math.PI)

        return retval
    }

    companion object {


        /**
         * Calculates barycentric weights for a collection of abscissa values.

         * @param z     double[] containing the ordinate values for which the barycentric weights are computed.
         * *
         * @return      double[] containing the resulting barycentric weights.
         */
        fun BarycentricWeights(z: DoubleArray): DoubleArray {

            val n = z.size

            val retval = DoubleArray(n)

            for (j in 0..n - 1) {
                var w = 1.0
                for (i in 0..n - 1) {
                    if (i != j) w *= z[j] - z[i]
                }
                retval[j] = 1.0 / w
            }

            return retval
        }


        /**
         * The main method.

         * @param args the arguments
         */
        @JvmStatic fun main(args: Array<String>) {

            val p = DoubleArray(3)
            p[0] = 6.0
            p[1] = -11.0
            p[2] = 6.0
            val P = Polynomial(p)

            val x = DoubleArray(3)
            x[0] = 1.0
            x[1] = 2.0
            x[2] = 3.0

            val f = DoubleArray(3)
            f[0] = 1.0
            f[1] = 8.0
            f[2] = 27.0
            val L = LagrangePolynomial(x, f)

            for (i in 0..20) {
                val z = 1.0 + i * 0.1
                println(P.evaluate(z).toString() + "  " + L.evaluate(z))
            }

        }
    }

}
