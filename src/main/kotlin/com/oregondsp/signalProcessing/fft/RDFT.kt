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

import kotlin.js.Math


/**

 * Class to calculate the discrete Fourier transform of a real sequence using symmetries and a complex DFT of half the length.

 *
 * This class is designed for efficient calculation of many discrete Fourier transforms of the
 * same length.  It is limited to transform lengths that are powers of two and are greater than or
 * equal to 32.  It uses the identity:

 *
 *
 * x[2n] + j*x[2n+1]  <->  X[k] + X[k+N/2] + j*W<sup>k</sup>*( X[k] - X[k+N/2] )
 *
 *
 *
 * and the fact that, for real sequences x[n], the transform exhibits conjugate symmetry:
 *
 *
 * X[k] = conjg( X[N-k] )
 *
 *
 * to calculate a real DFT with a complex DFT of half the length.
 *

 *
 *  Example of use:
 *
 *
 * <font face="courier">
 * int N &nbsp&nbsp&nbsp&nbsp= 16384;<BR></BR>
 * int log2N = 14;<BR></BR>
 * float[] x = new float[N];<BR></BR>
 * float[] X = new float[N];<BR></BR>
 * RDFT Xfm = new RDFT( log2N );<BR></BR>
 * <BR></BR>
 * // load data<BR></BR>
 * for ( int i = 0;  i < N;  i++ ) {<BR></BR>
 * &nbsp x[i] = ...<BR></BR>
 * }<BR></BR>
 * <BR></BR>
 * // evaluate transform of data<BR></BR>
 * Xfm.evaluate( x, X );<BR></BR>
</font> *
 *
 *
 *
 * The input sequence, x[n], is in natural order, and the output transform is in the packed form used
 * by Sorensen et al. (1987) for conjugate symmetric transforms:
 *
 *

 *
 *
 * <font face="courier">
 * __0_____1_____2_____3_____..._______N/2-1______N/2_______N/2+1______N/2+2____...____N-1 <BR></BR>
 * outXr(0)_Xr(1)_Xr(2)_Xr(3)___..._____Xr(N/2-1)__Xr(N/2)___Xi(N/2-1)__Xi(N/2-2)__...___Xi(1)
</font> *
 *
 *
 *  where outXr is the real part of the transform and outXi is the imaginary part.

 *
 * As long as the transform size does not change, the RDFT object does not need to be reinstantiated.
 * Consequently, the data arrays can be reloaded and the evaluate method invoked to compute additional
 * DFTs without incurring the cost of RDFT object instantiation.

 *
 *
 * The inverse DFT is calculated with a call to evaluateInverse():
 *
 *
 *
 * <font face="courier">
 * Xfm.evaluateInverse( X, x );
</font> *
 *

 *
 *
 * The input transform, X, is in conjugate symmetric packed form and the output sequence, x, is in natural order.
 *
 *
 *
 * This class also contains a convenience method to support convolution and correlation via the FFT.
 *

 *
 *  See "Real-valued Fast Fourier Transform Algorithms", Sorensen, H. V., et al., IEEE TRANSACTIONS ON
 * ACOUSTICS, SPEECH, AND SIGNAL PROCESSING, VOL. ASSP-35, NO. 6, JUNE 1987, pp. 849-863.

 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
class RDFT(log2N: Int) {

    private val N: Int
    private val N2: Int
    private val N4: Int

    private val xr: FloatArray
    private val xi: FloatArray
    private val Xr: FloatArray
    private val Xi: FloatArray

    private val dft: CDFT

    private val c: FloatArray
    private val s: FloatArray


    init {

        if (log2N < 4) throw IllegalArgumentException("DFT size must be >= 16")

        N = 1 shl log2N
        N2 = N / 2
        N4 = N / 4

        xr = FloatArray(N2)
        xi = FloatArray(N2)
        Xr = FloatArray(N2)
        Xi = FloatArray(N2)

        s = FloatArray(N4)
        c = FloatArray(N4)

        for (i in 0..N4 - 1) {
            s[i] = Math.sin(2.0 * Math.PI / N * i).toFloat()
            c[i] = Math.cos(2.0 * Math.PI / N * i).toFloat()
        }


        dft = CDFT(log2N - 1)

    }


    /**
     * Evaluates the DFT of a real sequence x.
     * @param x     float[] containing the real sequence in natural order.
     * *
     * @param X     float[] containing the transform of the sequence in conjugate symmetric packed form.
     */
    fun evaluate(x: FloatArray, X: FloatArray) {

        // Uses symmetries to perform the real length-N DFT with a special length-N set of butterflies
        // and one length-N/2 complex DFT.
        //
        // uses, specifically, the identity:  x[2n] + j*x[2n+1]  <->  X[k] + X[k+N/2] + j*W^k*( X[k] - X[k+N/2] )
        // and the fact that, for real sequences, X[k] = conjg( X[N-k] )

        for (i in 0..N2 - 1) {
            var j = i shl 1
            xr[i] = x[j++]
            xi[i] = x[j]
        }

        dft.evaluate(xr, xi, Xr, Xi)

        // special case at k = 0

        X[0] = Xr[0] + Xi[0]
        X[N2] = Xr[0] - Xi[0]

        // 1 <= k < N/4

        var N2pk = N2 + 1
        var N2mk = N2 - 1
        var Nmk = N - 1
        for (k in 1..N4 - 1) {

            val Xrk = Xr[k]
            val Xik = Xi[k]
            val XrN2mk = Xr[N2mk]
            val XiN2mk = Xi[N2mk]

            val Sr = (Xrk + XrN2mk) / 2
            val Si = (Xik - XiN2mk) / 2

            var Dr = (Xik + XiN2mk) / 2
            var Di = (XrN2mk - Xrk) / 2

            val tmp = c[k] * Dr + s[k] * Di
            Di = c[k] * Di - s[k] * Dr
            Dr = tmp

            X[k] = Sr + Dr
            X[Nmk] = Si + Di

            X[N2mk] = Sr - Dr
            X[N2pk] = Di - Si

            N2pk++
            N2mk--
            Nmk--
        }

        // special case at k = N/4

        //  cos( 2*pi/N * k ) = cos( pi/2 ) = 0
        //  sin( 2*pi/N * k ) = sin( pi/2 ) = 1

        X[N4] = Xr[N4]
        X[N2 + N4] = -Xi[N4]

    }


    /**
     * Evaluates the inverse DFT of a conjugate symmetric transform.
     * @param X     float[] containing the input transform of the sequence in conjugate symmetric packed form.
     * *
     * @param x     float[] containing the output real sequence in natural order.
     */
    fun evaluateInverse(X: FloatArray, x: FloatArray) {

        // Assumed input storage:
        //   0     1     2     3     ...       N/2-1      N/2       N/2+1      N/2+2    ...    N-1
        // outXr(0) outXr(1) outXr(2) outXr(3)   ...     outXr(N/2-1)  outXr(N/2)   outXi(N/2-1)  outXi(N/2-2)        outXi(1)


        // Uses symmetries to perform the real length-N inverse DFT with a special length-N set of butterflies
        // and one length-N/2 complex DFT.
        //
        // uses, specifically, the identity:  x[2n] + j*x[2n+1]  <->  X[k] + X[k+N/2] + j*W^k*( X[k] - X[k+N/2] )
        // and the fact that, for real sequences, X[k] = conjg( X[N-k] )


        // special case at k = 0

        Xr[0] = X[0] + X[N2]
        Xi[0] = X[0] - X[N2]

        // 1 <= k < N/4

        var N2pk = N2 + 1
        var N2mk = N2 - 1
        var Nmk = N - 1
        for (k in 1..N4 - 1) {

            val Xrk = X[k]
            val Xik = X[Nmk]
            val XrkpN2 = X[N2mk]
            val XikpN2 = -X[N2pk]

            val Dr = Xrk - XrkpN2
            val Di = Xik - XikpN2

            Xr[k] = Xrk + XrkpN2 - s[k] * Dr - c[k] * Di
            Xi[k] = Xik + XikpN2 + c[k] * Dr - s[k] * Di

            N2pk++
            N2mk--
            Nmk--
        }

        // special case at k = N/4

        //  cos( 2*pi/N * k ) = cos( pi/2 ) = 0
        //  sin( 2*pi/N * k ) = sin( pi/2 ) = 1

        Xr[N4] = 2.0f * X[N4]
        Xi[N4] = -2.0f * X[N2 + N4]

        // N/4 + 1  <=  k  <  N/2 - 1;

        //  cos( 2*pi/N * (N/4+m) ) = cos( 2*pi/N*m + pi/2 ) = -cos( 2*pi/N*(N/4-m) )
        //  sin( 2*pi/N * (N/4+m) ) = sin( 2*pi/N * (N/4-m) )

        N2pk = N2 + N4 + 1
        N2mk = N4 - 1
        Nmk = N - N4 - 1
        var reflect = N4 - 1
        for (k in N4 + 1..N2 - 1) {

            val Xrk = X[k]
            val Xik = X[Nmk]
            val XrkpN2 = X[N2mk]
            val XikpN2 = -X[N2pk]

            val Dr = Xrk - XrkpN2
            val Di = Xik - XikpN2

            Xr[k] = Xrk + XrkpN2 - s[reflect] * Dr + c[reflect] * Di
            Xi[k] = Xik + XikpN2 - c[reflect] * Dr - s[reflect] * Di

            N2pk++
            N2mk--
            Nmk--
            reflect--
        }

        dft.evaluate(Xr, Xi, xr, xi)

        x[0] = xr[0] / N
        x[1] = xi[0] / N

        var j = N2 - 1
        for (k in 1..N2 - 1) {
            var i = k shl 1
            x[i++] = xr[j] / N
            x[i] = xi[j] / N
            j--
        }

    }

    companion object {


        /**
         * Calculates the product of two conjugate symmetric dfts of the same length and stores the result in the second dft.

         *
         * This is a convenience method to support convolution and correlation operations using the Fast Fourier Transform.
         * Example of use to compute the convolution of two sequences:
         *
         *
         * <font face="courier">
         * int N &nbsp&nbsp&nbsp&nbsp= 1024;<BR></BR>
         * int log2N = 10;<BR></BR>
         * float[] x = new float[N];<BR></BR>
         * float[] y = new float[N];<BR></BR>
         * float[] X = new float[N];<BR></BR>
         * float[] Y = new float[N];<BR></BR>
         * RDFT Xfm = new RDFT( log2N );<BR></BR>
         * <BR></BR>
         * // load sequences<BR></BR>
         * for ( int i = 0;  i < N;  i++ ) {<BR></BR>
         * &nbsp x[i] = ...<BR></BR>
         * &nbsp y[i] = ...<BR></BR>
         * }<BR></BR>
         * <BR></BR>
         * // evaluate transforms of sequences<BR></BR>
         * Xfm.evaluate( x, X );<BR></BR>
         * Xfm.evaluate( y, Y );<BR></BR>
         * <BR></BR>
         * // product of transforms <BR></BR>
         * RDFT.dftProduct( X, Y, 1 );<BR></BR>
         * <BR></BR>
         * // inverse transform to obtain convolution<BR></BR>
         * float[] xy = new float[1024];<BR></BR>
         * RDFT.evaluateInverse( Y, xy );<BR></BR>
        </font> *
         *

         * @param kernel       first DFT
         * *
         * @param transform    second DFT before call, contains product after the call
         * *
         * @param sign         +1 if a convolution type product, -1 if a correlation type product
         */
        fun dftProduct(kernel: FloatArray, transform: FloatArray, sign: Float) {

            if (kernel.size != transform.size)
                throw IllegalArgumentException("kernel and transform arrays must have the same size")

            val n = kernel.size
            val half = n / 2
            transform[0] *= kernel[0]
            transform[half] *= kernel[half]

            var tmp: Float
            for (i in 1..half - 1) {
                val im = n - i
                tmp = kernel[i] * transform[i] - sign * kernel[im] * transform[im]
                transform[im] = kernel[i] * transform[im] + sign * kernel[im] * transform[i]
                transform[i] = tmp
            }

        }
    }

}
