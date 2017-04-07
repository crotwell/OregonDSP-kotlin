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
 * Package-private class implementing a length-8 complex DFT with a split-radix algorithm.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
internal class CDFTsr8
/**
 * Instantiates a new CDFTsr8.

 * @param xoffset  int specifying offset into the top-level length-N sequence array.
 * *
 * @param xstride  int specifying the stride of butterflies into the top-level length-N sequence array.
 * *
 * @param Xoffset  int specifying the offset into the length-N transform array.
 */
(
        /** Input sequence array indices.  */
        private val n0: Int, xstride: Int,
        /** Output transform array indices.  */
        private val m0: Int) : CDFTsr() {
    private val n1: Int
    private val n2: Int
    private val n3: Int
    private val n4: Int
    private val n5: Int
    private val n6: Int
    private val n7: Int
    private val m1: Int
    private val m2: Int
    private val m3: Int
    private val m4: Int
    private val m5: Int
    private val m6: Int
    private val m7: Int


    init {

        m = 3
        N = 8
        this.xoffset = n0
        this.xstride = xstride
        this.outXoffset = m0
        n1 = n0 + xstride
        n2 = n1 + xstride
        n3 = n2 + xstride
        n4 = n3 + xstride
        n5 = n4 + xstride
        n6 = n5 + xstride
        n7 = n6 + xstride
        m1 = m0 + 1
        m2 = m1 + 1
        m3 = m2 + 1
        m4 = m3 + 1
        m5 = m4 + 1
        m6 = m5 + 1
        m7 = m6 + 1

    }


    /**
     * Links the user-supplied input sequence and output transform arrays.

     * @param xr  float[] containing the input sequence real part.
     * *
     * @param xi  float[] containing the input sequence imaginary part.
     * *
     * @param Xr  float[] containing the output sequence real part.
     * *
     * @param Xi  float[] containing the output sequence imaginary part.
     */
    internal override fun link(xr: FloatArray, xi: FloatArray, Xr: FloatArray, Xi: FloatArray) {
        this.xr = xr
        this.xi = xi
        this.outXr = Xr
        this.outXi = Xi
    }


    /**
     * Evaluates the length-8 complex DFT.
     */
    internal override fun evaluate() {

        val T1r: Float
        val T1i: Float
        val T3r: Float
        val T3i: Float
        var Rr: Float
        var Ri: Float
        var Sr: Float
        var Si: Float

        // Length 2 DFT

        outXr[m0] = xr[n0] + xr[n4]
        outXi[m0] = xi[n0] + xi[n4]
        outXr[m1] = xr[n0] - xr[n4]
        outXi[m1] = xi[n0] - xi[n4]

        // length 4 dft

        // k = 0 butterfly

        Rr = xr[n2] + xr[n6]
        Ri = xi[n2] + xi[n6]
        Sr = xi[n6] - xi[n2]
        Si = xr[n2] - xr[n6]

        outXr[m2] = outXr[m0] - Rr
        outXi[m2] = outXi[m0] - Ri
        outXr[m3] = outXr[m1] + Sr
        outXi[m3] = outXi[m1] + Si

        outXr[m0] += Rr
        outXi[m0] += Ri
        outXr[m1] -= Sr
        outXi[m1] -= Si

        // Length 2 DFT

        outXr[m4] = xr[n1] + xr[n5]
        outXi[m4] = xi[n1] + xi[n5]
        outXr[m5] = xr[n1] - xr[n5]
        outXi[m5] = xi[n1] - xi[n5]

        // Length 2 DFT

        outXr[m6] = xr[n3] + xr[n7]
        outXi[m6] = xi[n3] + xi[n7]
        outXr[m7] = xr[n3] - xr[n7]
        outXi[m7] = xi[n3] - xi[n7]


        // length 8 dft


        // k = 0 butterfly

        Rr = outXr[m4] + outXr[m6]
        Ri = outXi[m4] + outXi[m6]
        Sr = outXi[m6] - outXi[m4]
        Si = outXr[m4] - outXr[m6]

        outXr[m4] = outXr[m0] - Rr
        outXi[m4] = outXi[m0] - Ri
        outXr[m6] = outXr[m2] + Sr
        outXi[m6] = outXi[m2] + Si

        outXr[m0] += Rr
        outXi[m0] += Ri
        outXr[m2] -= Sr
        outXi[m2] -= Si


        // k = 1 butterfly

        // T1 = Wk*O1
        // T3 = W3k*O3

        T1r = SQRT2BY2 * (outXr[m5] + outXi[m5])
        T1i = SQRT2BY2 * (outXi[m5] - outXr[m5])
        T3r = SQRT2BY2 * (outXi[m7] - outXr[m7])
        T3i = -SQRT2BY2 * (outXi[m7] + outXr[m7])

        // R = T1 + T3
        // S = i*(T1 - T3)

        Rr = T1r + T3r
        Ri = T1i + T3i
        Sr = T3i - T1i
        Si = T1r - T3r

        outXr[m5] = outXr[m1] - Rr
        outXi[m5] = outXi[m1] - Ri
        outXr[m7] = outXr[m3] + Sr
        outXi[m7] = outXi[m3] + Si

        outXr[m1] += Rr
        outXi[m1] += Ri
        outXr[m3] -= Sr
        outXi[m3] -= Si

    }

    companion object {

        /** Constant twiddle factor.  */
        val SQRT2BY2 = (Math.sqrt(2.0) / 2.0).toFloat()
    }

}
