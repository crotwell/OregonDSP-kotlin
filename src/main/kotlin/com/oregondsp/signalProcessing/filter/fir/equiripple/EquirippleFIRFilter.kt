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


import com.oregondsp.signalProcessing.fft.RDFT
import com.oregondsp.signalProcessing.filter.fir.OverlapAdd
import kotlin.math.*
import kotlin.js.Math.random


/**
 * Base class for equiripple FIR filters designed with the Parks-McClellan algorithm.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
abstract class EquirippleFIRFilter
/**
 * Instantiates a new equiripple FIR filter.

 * @param numBands    int specifying the number of pass and stop bands.
 * *
 * @param N           int specifying the design _order of the filter.
 * *
 * @param Nc          int specifying the number of FIR filter coefficients
 */
(
        /** int specifying the number of bands  */
        protected var numBands: Int          // number of bands
        ,
        /** int specifying the number of approximating functions in the Remez algorithm  */
        protected var N: Int                 // number of approximating degrees of freedom
        ,
        /** int containing the number of filter coefficients  */
        protected var Nc: Int                // number of coefficients
) {

    /** double[][] specifying band edge information  */
    protected var bands: Array<DoubleArray>             // specifies band edges

    /** float[] containing the FIR filter coefficients  */
    var _coefficients: FloatArray?  = null      // filter coefficients

    fun getCoefficients(): FloatArray {
        return _coefficients?.copyOf() ?: throw RuntimeException("Should not happen, access to coefficients before initialized.")
    }



    /** An OverlapAdd instance that can be used to filter data with the filter.  */
    protected var implementation: OverlapAdd? = null


    init {
        bands = Array(numBands) { DoubleArray(2) }
    }


    /**
     * Method to create the design grid.

     * @return    DesignGrid object used by the Remez exchange algorithm
     */
    protected fun createGrid(): DesignGrid {

        val G = DesignGrid()

        //  initial guess for extreme points - need N + 1 approximately equally spaced among pass and stop bands
        //    include band edges, 0 and pi

        val nextrema = IntArray(numBands)  // defines allocation of extrema among bands

        var totalBandwidth = 0.0
        for (ib in 0..numBands - 1) totalBandwidth += bands[ib][1] - bands[ib][0]

        val m = N + 1 - 2 * numBands
        var np = 0
        var largestBand = 0
        var nmax = 0
        for (ib in 0..numBands - 1) {
            val B = bands[ib][1] - bands[ib][0]
            nextrema[ib] = round(m * B / totalBandwidth).toInt() + 2
            if (nextrema[ib] > nmax) {
                nmax = nextrema[ib]
                largestBand = ib
            }
            np += nextrema[ib]
        }

        // Add or delete extrema to largest band to assure exactly N+1 extrema are specified
        //   Because of rounding the last step may not produce the correct number

        while (np < N + 1) {
            nextrema[largestBand]++
            np++
        }
        while (np > N + 1) {
            nextrema[largestBand]--
            np--
        }

        // set up dense grid and initial estimate of extrema

        G.bandEdgeIndices = IntArray(numBands * 2)
        G.extremaIndices = IntArray(N + 1)

        //   grid sampling proportional to band widths

        val gridArray = ArrayList<Double>()
        var gridpt = 0
        var extremum = 0
        var bandEdgeCount = 0
        var perturbation: Int
        for (ib in 0..numBands - 1) {
            val B = bands[ib][1] - bands[ib][0]
            val n = 1 + (nextrema[ib] - 1) * DesignGrid.GRIDDENSITY
            val dB = B / (n - 1)
            val base = bands[ib][0]
            for (i in 0..n - 1) {

                val Omega = base + dB * i
                gridArray.add(Omega)

                if (i % DesignGrid.GRIDDENSITY == 0) {
                    if (i != 0 && i != n - 1)
                        perturbation = (floor(random()*3)).roundToInt() - 1
                    else
                        perturbation = 0
                    G.extremaIndices[extremum++] = gridpt + perturbation
                }
                if (i == 0 || i == n - 1) {
                    G.bandEdgeIndices[bandEdgeCount] = gridpt
                    bandEdgeCount++
                }

                gridpt++
            }
        }
        G.gridSize = gridArray.size
        G.grid = DoubleArray(G.gridSize)
        G.X = DoubleArray(G.gridSize)
        G.H = DoubleArray(G.gridSize)
        G.W = DoubleArray(G.gridSize)
        for (i in 0..G.gridSize - 1) {
            G.grid[i] = gridArray[i]
            G.X[i] = cos(G.grid[i] * PI)
        }

        return G
    }


    /**
     * Method made concrete by TypeI-IV filters to populate the unified weight and objective functions.

     * @param G      DesignGrid object used by the Remez exchange algorithm.
     */
    @JsName("populateGrid")
    internal abstract fun populateGrid(G: DesignGrid)


    /**
     * Method made concrete by specific filter classes (e.g. EquirippleLowpass) to specify the objective function.

     * @param Omega      double specifying the normalized frequency at which the desired response is evaluated.
     * *
     * @return           double containing the desired response at this frequency.
     */
    @JsName("desiredResponse")
    internal abstract fun desiredResponse(Omega: Double): Double


    /**
     * Weight.

     * @param Omega      double specifying the normalized frequency at which the weight function is evaluated.
     * *
     * @return           double containing the weight function at this frequency.
     */
    @JsName("weight")
    internal abstract fun weight(Omega: Double): Double


    /**
     * Method to interpret cosine basis coefficients as TypeI-TypeIV FIR filter coefficients.

     * @param coefficients     float[] containing the cosine sequence coefficients
     * *
     * @return                 float[] containing the corresponding FIR filter coefficients.
     */
    @JsName("interpretCoefficients")
    internal abstract fun interpretCoefficients(coefficients: FloatArray): FloatArray


    /**
     * Method to generate cosine basis coefficients from response function on a dense grid.
     */
    fun generateCoefficients() {
        val G = createGrid()
        populateGrid(G)
        EquirippleDesigner.remez(G)
        _coefficients = interpretCoefficients(EquirippleDesigner.calculateCoefficients(G, Nc))
    }




    /**
     * Method to provide a new OverlapAdd instance to implement the filter.

     * @param blockSize the block size
     * *
     * @return the implementation
     */
    @JsName("getImplementation")
    fun getImplementation(blockSize: Int): OverlapAdd {
        return OverlapAdd(getCoefficients(), blockSize)
    }


    /**
     * Method to filter a fixed-length sequence with this filter.

     * @param x       float[] containing the input sequence.
     * *
     * @return        float[] containing the resulting filtered sequence.
     */
    @JsName("filter")
    fun filter(x: FloatArray): FloatArray {

        var nfft = 16
        var log2nfft = 4
        val coefficients = getCoefficients()
        val n = x.size + coefficients.size - 1
        while (nfft < n) {
            nfft *= 2
            log2nfft++
        }

        val fft = RDFT(log2nfft)
        val tmp = FloatArray(nfft)
        val transform = FloatArray(nfft)
        var kernel = FloatArray(nfft)

        for (i in x.indices) {
            tmp[i] = x[i]
        }
        fft.evaluate(tmp, transform)

        //Arrays.fill(tmp, 0.0f)
        for (i in tmp.indices)
            tmp[i] = 0.0F
        for (i in coefficients.indices) {
            tmp[i] = coefficients[i]
        }
        fft.evaluate(tmp, kernel)

        RDFT.dftProduct(kernel, transform, 1.0f)
        fft.evaluateInverse(transform, tmp)

        // trim off trailing zeros

        kernel = FloatArray(n)
        for (i in 0..n) {
            kernel[i] = tmp[i]
        }

        return kernel
    }


    /**
     * Method to determine whether one double is close to another.

     * @param x     the first double
     * *
     * @param y     the second double
     * *
     * @return true, if the numbers are close, false otherwise
     */
    @JsName("LTE")
    protected fun LTE(x: Double, y: Double): Boolean {
        var retval = false

        if (x < y) retval = true

        if (abs(x - y) < MACHINETOLERANCE) retval = true

        return retval
    }

    companion object {

        /** Constant specifying a tolerance for checking band edge inclusion in the design grid.  */
        protected var MACHINETOLERANCE = 1.0E-6
    }

}
