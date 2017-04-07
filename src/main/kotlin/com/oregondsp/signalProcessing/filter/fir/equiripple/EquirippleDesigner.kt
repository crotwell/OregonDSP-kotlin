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

package com.oregondsp.signalProcessing.filter.fir.equiripple

import java.util.ArrayList


import com.oregondsp.signalProcessing.fft.RDFT
import com.oregondsp.signalProcessing.filter.LagrangePolynomial
import kotlin.js.Math


/**
 * Implements the Parks-McClellan algorithm for designing equiripple FIR digital filters.

 *
 * See

 *
 * Chebyshev Approximation for Nonrecursive Digital Filters with Linear Phase,
 * Thomas W. Parks and James H. McClellan (1972), IEEE Transactions on Circuit Theory, Vol. CT-19, no. 2,
 * pp. 184-194.

 *
 * A Unified Approach to the Design of Optimum Linear Phase FIR Digital Filters,
 * James H. McClellan and Thomas W. Parks (1973), IEEE Transactions on Circuit Theory, Vol. CT-20,
 * No. 6, pp. 697-701.

 *
 * A Program for the Design of Linear Phase Finite Impulse Response Digital Filters,
 * Thomas W. Parks and James H. McClellan (1972), IEEE Transactions on Audio and Electroacoustics,
 * Vol. AU-20 no. 3, pp. 195-199.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
internal object EquirippleDesigner {

    /** Constant specifying the maximum number of iterations of the Remez exchange algorithm  */
    private val MAXITER = 25


    /**
     * Remez exchange algorithm.

     * @param G     DesignGrid object containing the finite frequency-sampling grid used by the Remez
     * *                exchange algorithm.
     */
    fun remez(G: DesignGrid) {

        val nextrema = G.extremaIndices!!.size

        val newExtrema = ArrayList<Int>()
        val E = DoubleArray(G.gridSize)
        val GA = DoubleArray(G.gridSize)


        var niter = 0

        do {

            val delta = computeDelta(G)

            println("delta: " + delta)

            val LP = constructInterpolatingPolynomial(G, delta)

            //  Compute current approximant (GA) and error function (E) on grid

            for (i in 0..G.gridSize - 1) {
                GA[i] = LP.evaluate(G.X!![i])
                E[i] = GA[i] - G.H!![i]
            }

            // Search for new extrema starting from old extrema

            newExtrema.clear()

            var change = 0

            for (currentExtremum in 0..nextrema - 1) {

                val currentGridPt = G.extremaIndices!![currentExtremum]
                val s = sgn(E[currentGridPt])

                // search forward

                var ptr = currentGridPt + 1
                if (ptr < G.gridSize) {
                    while (sgn(E[ptr] - E[ptr - 1]) == s) {
                        ptr++
                        if (ptr == G.gridSize) break
                    }
                }
                ptr--

                if (ptr == currentGridPt) {

                    // forward search failed, try backward search

                    ptr = currentGridPt - 1
                    if (ptr >= 0) {
                        while (sgn(E[ptr] - E[ptr + 1]) == s) {
                            ptr--
                            if (ptr < 0) break
                        }
                    }
                    ptr++

                }

                newExtrema.add(ptr)
                if (ptr != currentGridPt) change++
            }

            // test for exchanges at 0 and pi


            if (G.containsZero && G.containsPi) {

                val gridPi = G.gridSize - 1

                if (newExtrema.contains(0)) {

                    if (!newExtrema.contains(gridPi)) {
                        if (sgn(E[gridPi]) != sgn(E[G.extremaIndices!![nextrema - 1]])) {
                            if (Math.abs(E[gridPi]) > Math.abs(E[0])) {
                                newExtrema.removeAt(0)
                                newExtrema.add(gridPi)
                                change++
                            }
                        }
                    }
                } else {

                    if (newExtrema.contains(gridPi)) {

                        if (sgn(E[0]) != sgn(E[G.extremaIndices!![0]])) {
                            if (Math.abs(E[0]) > Math.abs(E[gridPi])) {
                                newExtrema.removeAt(newExtrema.size - 1)
                                newExtrema.add(0, 0)
                                change++
                            }
                        }
                    }

                }

            }

            if (change == 0) break

            // exchange extrema

            for (i in 0..nextrema - 1) {
                G.extremaIndices[i] = newExtrema[i]
            }

            niter++
        } while (niter < MAXITER)

    }


    /**
     * Method to compute the Linfinity norm best approximation error on the current set of extrema.

     * @param G      DesignGrid instance, which contains a list of current extrema.
     * *
     * @return       double containing the error on the current set of extrema.
     */
    fun computeDelta(G: DesignGrid): Double {

        val nextrema = G.extremaIndices!!.size
        val extrema = DoubleArray(nextrema)
        for (i in 0..nextrema - 1) {
            extrema[i] = G.X!![G.extremaIndices!![i]]
        }
        val gamma = LagrangePolynomial.BarycentricWeights(extrema)

        var num = 0.0
        var denom = 0.0
        var s = 1.0
        for (i in 0..nextrema - 1) {
            val j = G.extremaIndices!![i]
            num += gamma[i] * G.H!![j]
            denom += s * gamma[i] / G.W!![j]
            s = -s
        }

        return num / denom
    }


    /**
     * Constructs the Lagrange polynomial on the set of extrema.

     * @param  G        DesignGrid instance which contains the current set of extrema.
     * *
     * @param  delta    Current deviation on that set of extrema.
     * *
     * @return          Lagrange polynomial instance that interpolates the extrema.
     */
    fun constructInterpolatingPolynomial(G: DesignGrid, delta: Double): LagrangePolynomial {

        val extremaSubset = DoubleArray(G.extremaIndices!!.size - 1)
        val n = extremaSubset.size
        val x = DoubleArray(n)
        val f = DoubleArray(n)
        var s = 1.0
        for (i in 0..n - 1) {
            val j = G.extremaIndices!![i]
            x[i] = G.X!![j]
            f[i] = G.H!![j] - s * delta / G.W!![j]
            s = -s
        }

        return LagrangePolynomial(x, f)
    }


    /**
     * Calculates coefficients of the best Chebyshev approximation out of a cosine basis.

     * @param   G    DesignGrid instance that contains a list of the current extrema.
     * *
     * @param   Nc   The number of coefficients of the corresponding FIR filter.
     * *
     * @return       float[] containing the coefficients of the FIR filter.
     */
    fun calculateCoefficients(G: DesignGrid, Nc: Int): FloatArray {

        val LP = constructInterpolatingPolynomial(G, computeDelta(G))

        var log2nfft = 6
        var nfft = 64
        while (nfft < Nc) {
            nfft *= 2
            log2nfft++
        }
        val X = FloatArray(nfft)
        val x = FloatArray(nfft)
        for (i in 0..nfft / 2) {
            X[i] = LP.evaluate(Math.cos(2.0 * Math.PI * i.toDouble() / nfft)).toFloat()
        }

        val dft = RDFT(log2nfft)
        dft.evaluateInverse(X, x)

        return x
    }


    /**
     * Method to compute the sign of a double.

     * @param x    the double in question
     * *
     * @return     int ( = 1 if x > 0, = 0 if x = 0, = -1 if x < 0 )
     */
    fun sgn(x: Double): Int {
        if (x > 0.0)
            return 1
        else if (x < 0.0)
            return -1
        else
            return 0
    }


}
