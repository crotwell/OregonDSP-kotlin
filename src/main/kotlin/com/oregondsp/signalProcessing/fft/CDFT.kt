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


package com.oregondsp.signalProcessing.fft


/**
 * Class to calculate the complex discrete Fourier transform of a complex sequence and its inverse using the split-radix algorithm.

 *
 * This class is designed for efficient calculation of many discrete Fourier transforms of the
 * same length.  It is limited to transform lengths that are powers of two and greater than or
 * equal to 32.  The class recursively constructs and links smaller DFTs with hard-wired array indices
 * to minimize index calculations during the overall DFT evaluation.  This approach may produce large run-time
 * images for very large DFTs (> 32768).  Special hand-coded implementations of length 8 and 16 DFTs eliminate
 * many unnecessary calculations.  The code uses precomputed sine and cosine tables and does not implement
 * in-place calculations in order to eliminate the bit reversal step.  Consequently, this implementation
 * trades memory for speed.

 *
 *  Example of use:
 *
 *
 * <font face="courier">
 * int N &nbsp&nbsp&nbsp&nbsp&nbsp= 16384;<BR></BR>
 * int log2N &nbsp= 14;<BR></BR>
 * float[] xr = new float[N];<BR></BR>
 * float[] xi = new float[N];<BR></BR>
 * float[] outXr = new float[N];<BR></BR>
 * float[] outXi = new float[N];<BR></BR>
 * CDFT Xfm = new CDFT( log2N );<BR></BR>
 * <BR></BR>
 * // load data<BR></BR>
 * for ( int i = 0;  i < N;  i++ ) {<BR></BR>
 * &nbsp xr[i] = ...<BR></BR>
 * &nbsp xi[i] = ...<BR></BR>
 * }<BR></BR>
 * <BR></BR>
 * // evaluate transform of data<BR></BR>
 * Xfm.evaluate( xr, xi, outXr, outXi );<BR></BR>
</font> *
 *

 *
 * The real and imaginary parts of the transform are stored in outXr and outXi in natural order, with the zeroth
 * discrete frequency value in outXr(0) and outXi(0), and the N-1st value ( 2*pi*(N-1)/N ) in outXr(N-1) and outXi(N-1).
 *

 *
 * As long as the transform size does not change, the CDFT object does not need to be reinstantiated.
 * Consequently, the data arrays can be reloaded and the evaluate method invoked to compute additional
 * DFTs without incurring the cost of CDFT object instantiation.

 *
 * It may happen in some applications that the array arguments in the evaluate() and evaluateInverse()
 * methods never change, i.e. the same arrays are used repeatedly.  Since this implementation is recursive,
 * the input and output arrays are recursively linked down the chain of smaller DFTs that implement the full
 * DFT.  This linking operation can be avoided when the arguments to evaluate() and evaluateInverse() never vary.
 * For this circumstance an alternative constructor is provided, that links the input and output arrays at
 * construction time (for a slight performance improvement).  To avoid relinking arrays, this constructor should
 * be paired with the evaluate() and evaluateInverse() methods that have NO arguments.  Example:

 *
 *
 * <font face="courier">
 * CDFT Xfm = new CDFT( xr, xi, outXr, outXi, log2N );<BR></BR>
 * <BR></BR>
 * // load data<BR></BR>
 * for ( int i = 0;  i < N;  i++ ) {<BR></BR>
 * &nbsp xr[i] = ...<BR></BR>
 * &nbsp xi[i] = ...<BR></BR>
 * }<BR></BR>
 * <BR></BR>
 * // evaluate transform of data<BR></BR>
 * Xfm.evaluate();<BR></BR>
</font> *
 *

 *
 * For the inverse transform in this usage, the roles of (xr,xi) and (outXr,outXi) are reversed.  The pair
 * (xr,xi) contains the transform real and imaginary parts in natural order, and upon execution of
 * evaluateInverse(), the pair (outXr,outXi) contains the real and imaginary parts of the corresponding sequence
 * (inverse transform).
 *

 *
 * See "On Computing the Split-Radix FFT", Sorensen, H. V., Heideman, M. T. and Burrus, C. S.
 * IEEE TRANSACTIONS ON ACOUSTICS, SPEECH, AND SIGNAL PROCESSING, VOL. ASSP-34, NO. 1,
 * FEBRUARY, 1986, pp. 152-156.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
class CDFT {

    private lateinit var yr: FloatArray
    private lateinit var yi: FloatArray
    private var arraysUnlinked: Boolean = false

    private lateinit var c: FloatArray
    private lateinit var c3: FloatArray
    private lateinit var s: FloatArray
    private lateinit var s3: FloatArray

    internal val N: Int
    internal val log2N: Int

    private val dft: CDFTsr


    /**
     * Constructs a CDFT instance without references to sequence and transform arrays
     * @param log2N       base-2 logarithm of the length of the transform
     */
    constructor(log2N: Int) {

        if (log2N < 3) throw IllegalArgumentException("DFT size must be >= 8")
        arraysUnlinked = true

        this.log2N = log2N
        N = 1 shl log2N

        createTable()

        if (log2N == 3)
            dft = CDFTsr8(0, 1, 0)
        else if (log2N == 4)
            dft = CDFTsr16(0, 1, 0)
        else if (log2N >= 5) {
            dft = CDFTsr(log2N, c, c3, s, s3)
        } else
            throw IllegalArgumentException("unknown log2N size, must be >=3 but was: "+log2N);
    }


    /**
     * evaluates the DFT with specified sequence and transform arrays
     * @param xr          float array containing sequence real part
     * *
     * @param xi          float array containing sequence imaginary part
     * *
     * @param Xr          float array containing transform real part
     * *
     * @param Xi          float array containing transform imaginary part
     */
    fun evaluate(xr: FloatArray, xi: FloatArray, Xr: FloatArray, Xi: FloatArray) {
        this.yr = Xr
        this.yi = Xi
        dft.link(xr, xi, Xr, Xi)
        arraysUnlinked = false
        dft.evaluate()
    }


    /**
     * evaluates the inverse DFT with specified transform and sequence arrays
     * @param Xr          float array containing transform real part
     * *
     * @param Xi          float array containing transform imaginary part
     * *
     * @param xr          float array containing sequence real part
     * *
     * @param xi          float array containing sequence imaginary part
     */
    fun evaluateInverse(Xr: FloatArray, Xi: FloatArray, xr: FloatArray, xi: FloatArray) {
        this.yr = xr
        this.yi = xi
        dft.link(Xr, Xi, xr, xi)
        arraysUnlinked = false
        evaluateInverse()
    }


    /**
     * constructs a CDFT instance with references to sequence and transform arrays
     * @param xr          float array containing sequence real part on forward evaluation,
     * *                    transform real part on inverse evaluation
     * *
     * @param xi          float array containing sequence imaginary part on forward evaluation,
     * *                    transform imaginary part on inverse evaluation
     * *
     * @param yr          float array containing transform real part on forward evaluation,
     * *                    sequence real part on inverse evaluation
     * *
     * @param yi          float array containing transform imaginary part on forward evaluation,
     * *                    sequence imaginary part on inverse evaluation
     * *
     * @param log2N       base-2 logarithm of the length of the transform
     */
    constructor(xr: FloatArray, xi: FloatArray, yr: FloatArray, yi: FloatArray, log2N: Int) :this(log2N) {

        this.yr = yr
        this.yi = yi

        dft.link(xr, xi, yr, yi)
        arraysUnlinked = false

    }


    /**
     * evaluates the DFT assuming sequence and transformed arrays have been linked at construction time

     */
    fun evaluate() {
        if (arraysUnlinked)
            throw IllegalStateException("Sequence and transform arrays are not linked")
        dft.evaluate()
    }


    /**
     * evaluates the inverse DFT assuming the sequence and transform arrays have been linked at construction time
     */
    fun evaluateInverse() {

        if (arraysUnlinked)
            throw IllegalStateException("Sequence and transform arrays are not linked")

        dft.evaluate()

        val scale = 1.0f / N.toFloat()
        val N2 = N / 2

        yr[0] *= scale
        yi[0] *= scale
        yr[N2] *= scale
        yi[N2] *= scale

        var i = 1
        var j = N - 1

        var tmp: Float

        while (i < j) {
            tmp = yr[i]
            yr[i] = yr[j] * scale
            yr[j] = tmp * scale
            tmp = yi[i]
            yi[i] = yi[j] * scale
            yi[j] = tmp * scale

            i++
            j--
        }

    }


    private fun createTable() {

        val N8 = N / 8

        c = FloatArray(N8)
        c3 = FloatArray(N8)
        s = FloatArray(N8)
        s3 = FloatArray(N8)

        for (i in 0..N8 - 1) {
            c[i] = Math.cos(2.0 * Math.PI * i.toDouble() / N).toFloat()
            c3[i] = Math.cos(2.0 * Math.PI * 3.0 * i.toDouble() / N).toFloat()
            s[i] = -Math.sin(2.0 * Math.PI * i.toDouble() / N).toFloat()
            s3[i] = -Math.sin(2.0 * Math.PI * 3.0 * i.toDouble() / N).toFloat()
        }

    }

    companion object {


        /**
         * Convenience method to multiply two complex transforms of the same size.
         * @param Xr     float array containing the real part of the first transform
         * *
         * @param Xi     float array containing the imaginary part of the first transform
         * *
         * @param Yr     float array containing the real part of the second transform before call, real part of the product after call
         * *
         * @param Yi     float array containing the imaginary part of the second transform before call, imaginary part of the product after call
         * *
         * @param sign   +1 for convolution type product, -1 for correlation type product
         */
        fun dftProduct(Xr: FloatArray, Xi: FloatArray, Yr: FloatArray, Yi: FloatArray, sign: Float) {

            if (Xr.size != Yr.size || Xi.size != Yi.size || Xr.size != Xi.size)
                throw IllegalArgumentException("Transform array lengths are not equal")

            var tmp: Float
            for (i in Xr.indices) {
                tmp = Xr[i] * Yr[i] - sign * Xi[i] * Yi[i]
                Yi[i] = Xr[i] * Yi[i] + sign * Xi[i] * Yr[i]
                Yr[i] = tmp
            }

        }
    }

}
