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
 * Package-private class implementing a length-16 complex DFT with a split-radix algorithm.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
internal class CDFTsr16
/**
 * Instantiates a new CDFTsr16 instance.

 * @param dataOffset       int specifying offset into the top-level length-N sequence array.
 * *
 * @param dataStride       int specifying the stride of butterflies into the top-level length-N sequence array.
 * *
 * @param transformOffset  int specifying the offset into the length-N transform array.
 */
(dataOffset: Int, dataStride: Int, transformOffset: Int) : CDFTsr() {


    /** Input sequence indices  */
    private val n0: Int
    private val n1: Int
    private val n2: Int
    private val n3: Int
    private val n4: Int
    private val n5: Int
    private val n6: Int
    private val n7: Int
    private val n8: Int
    private val n9: Int
    private val n10: Int
    private val n11: Int
    private val n12: Int
    private val n13: Int
    private val n14: Int
    private val n15: Int

    /** Output transform indices  */
    private val m0: Int
    private val m1: Int
    private val m2: Int
    private val m3: Int
    private val m4: Int
    private val m5: Int
    private val m6: Int
    private val m7: Int
    private val m8: Int
    private val m9: Int
    private val m10: Int
    private val m11: Int
    private val m12: Int
    private val m13: Int
    private val m14: Int
    private val m15: Int


    init {

        m = 4
        N = 16
        xoffset = dataOffset
        xstride = dataStride
        outXoffset = transformOffset

        n0 = xoffset
        n1 = n0 + xstride
        n2 = n1 + xstride
        n3 = n2 + xstride
        n4 = n3 + xstride
        n5 = n4 + xstride
        n6 = n5 + xstride
        n7 = n6 + xstride
        n8 = n7 + xstride
        n9 = n8 + xstride
        n10 = n9 + xstride
        n11 = n10 + xstride
        n12 = n11 + xstride
        n13 = n12 + xstride
        n14 = n13 + xstride
        n15 = n14 + xstride

        m0 = outXoffset
        m1 = m0 + 1
        m2 = m1 + 1
        m3 = m2 + 1
        m4 = m3 + 1
        m5 = m4 + 1
        m6 = m5 + 1
        m7 = m6 + 1
        m8 = m7 + 1
        m9 = m8 + 1
        m10 = m9 + 1
        m11 = m10 + 1
        m12 = m11 + 1
        m13 = m12 + 1
        m14 = m13 + 1
        m15 = m14 + 1

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
     * Evaluates the length-16 complex DFT.
     */
    internal override fun evaluate() {

        var T1r: Float
        var T1i: Float
        var T3r: Float
        var T3i: Float
        var Rr: Float
        var Ri: Float
        var Sr: Float
        var Si: Float


        // Length 2 DFT

        outXr[m0] = xr[n0] + xr[n8]
        outXi[m0] = xi[n0] + xi[n8]
        outXr[m1] = xr[n0] - xr[n8]
        outXi[m1] = xi[n0] - xi[n8]

        // length 4 dft


        // k = 0 butterfly

        Rr = xr[n4] + xr[n12]
        Ri = xi[n4] + xi[n12]
        Sr = xi[n12] - xi[n4]
        Si = xr[n4] - xr[n12]

        outXr[m2] = outXr[m0] - Rr
        outXi[m2] = outXi[m0] - Ri
        outXr[m3] = outXr[m1] + Sr
        outXi[m3] = outXi[m1] + Si

        outXr[m0] += Rr
        outXi[m0] += Ri
        outXr[m1] -= Sr
        outXi[m1] -= Si


        // Length 2 DFT

        outXr[m4] = xr[n2] + xr[n10]
        outXi[m4] = xi[n2] + xi[n10]
        outXr[m5] = xr[n2] - xr[n10]
        outXi[m5] = xi[n2] - xi[n10]

        // Length 2 DFT

        outXr[m6] = xr[n6] + xr[n14]
        outXi[m6] = xi[n6] + xi[n14]
        outXr[m7] = xr[n6] - xr[n14]
        outXi[m7] = xi[n6] - xi[n14]


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


        // all other butterflies


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


        // Length 2 DFT

        outXr[m8] = xr[n1] + xr[n9]
        outXi[m8] = xi[n1] + xi[n9]
        outXr[m9] = xr[n1] - xr[n9]
        outXi[m9] = xi[n1] - xi[n9]

        // length 4 dft


        // k = 0 butterfly

        Rr = xr[n5] + xr[n13]
        Ri = xi[n5] + xi[n13]
        Sr = xi[n13] - xi[n5]
        Si = xr[n5] - xr[n13]

        outXr[m10] = outXr[m8] - Rr
        outXi[m10] = outXi[m8] - Ri
        outXr[m11] = outXr[m9] + Sr
        outXi[m11] = outXi[m9] + Si

        outXr[m8] += Rr
        outXi[m8] += Ri
        outXr[m9] -= Sr
        outXi[m9] -= Si


        // Length 2 DFT

        outXr[m12] = xr[n3] + xr[n11]
        outXi[m12] = xi[n3] + xi[n11]
        outXr[m13] = xr[n3] - xr[n11]
        outXi[m13] = xi[n3] - xi[n11]


        // length 4 dft


        // k = 0 butterfly

        Rr = xr[n7] + xr[n15]
        Ri = xi[n7] + xi[n15]
        Sr = xi[n15] - xi[n7]
        Si = xr[n7] - xr[n15]

        outXr[m14] = outXr[m12] - Rr
        outXi[m14] = outXi[m12] - Ri
        outXr[m15] = outXr[m13] + Sr
        outXi[m15] = outXi[m13] + Si

        outXr[m12] += Rr
        outXi[m12] += Ri
        outXr[m13] -= Sr
        outXi[m13] -= Si


        // length 16 dft


        // k = 0 butterfly

        Rr = outXr[m8] + outXr[m12]
        Ri = outXi[m8] + outXi[m12]
        Sr = outXi[m12] - outXi[m8]
        Si = outXr[m8] - outXr[m12]

        outXr[m8] = outXr[m0] - Rr
        outXi[m8] = outXi[m0] - Ri
        outXr[m12] = outXr[m4] + Sr
        outXi[m12] = outXi[m4] + Si

        outXr[m0] += Rr
        outXi[m0] += Ri
        outXr[m4] -= Sr
        outXi[m4] -= Si


        // all other butterflies

        // k = 1
        // T1 = Wk*O1
        // T3 = W3k*O3

        T1r = C_1_16 * outXr[m9] + C_3_16 * outXi[m9]
        T1i = C_1_16 * outXi[m9] - C_3_16 * outXr[m9]
        T3r = C_3_16 * outXr[m13] + C_1_16 * outXi[m13]
        T3i = C_3_16 * outXi[m13] - C_1_16 * outXr[m13]

        // R = T1 + T3
        // S = i*(T1 - T3)

        Rr = T1r + T3r
        Ri = T1i + T3i
        Sr = T3i - T1i
        Si = T1r - T3r

        outXr[m9] = outXr[m1] - Rr
        outXi[m9] = outXi[m1] - Ri
        outXr[m13] = outXr[m5] + Sr
        outXi[m13] = outXi[m5] + Si

        outXr[m1] += Rr
        outXi[m1] += Ri
        outXr[m5] -= Sr
        outXi[m5] -= Si

        // k = 2
        // T1 = Wk*O1
        // T3 = W3k*O3

        T1r = SQRT2BY2 * (outXr[m10] + outXi[m10])
        T1i = SQRT2BY2 * (outXi[m10] - outXr[m10])
        T3r = SQRT2BY2 * (outXi[m14] - outXr[m14])
        T3i = -SQRT2BY2 * (outXi[m14] + outXr[m14])

        // R = T1 + T3
        // S = i*(T1 - T3)

        Rr = T1r + T3r
        Ri = T1i + T3i
        Sr = T3i - T1i
        Si = T1r - T3r

        outXr[m10] = outXr[m2] - Rr
        outXi[m10] = outXi[m2] - Ri
        outXr[m14] = outXr[m6] + Sr
        outXi[m14] = outXi[m6] + Si

        outXr[m2] += Rr
        outXi[m2] += Ri
        outXr[m6] -= Sr
        outXi[m6] -= Si

        // k = 3
        // T1 = Wk*O1
        // T3 = W3k*O3

        T1r = C_3_16 * outXr[m11] + C_1_16 * outXi[m11]
        T1i = C_3_16 * outXi[m11] - C_1_16 * outXr[m11]
        T3r = -C_1_16 * outXr[m15] - C_3_16 * outXi[m15]
        T3i = -C_1_16 * outXi[m15] + C_3_16 * outXr[m15]

        // R = T1 + T3
        // S = i*(T1 - T3)

        Rr = T1r + T3r
        Ri = T1i + T3i
        Sr = T3i - T1i
        Si = T1r - T3r

        outXr[m11] = outXr[m3] - Rr
        outXi[m11] = outXi[m3] - Ri
        outXr[m15] = outXr[m7] + Sr
        outXi[m15] = outXi[m7] + Si

        outXr[m3] += Rr
        outXi[m3] += Ri
        outXr[m7] -= Sr
        outXi[m7] -= Si

    }

    companion object {

        /** Constant twiddle factor   */
        val C_1_16 = Math.cos(2.0 * Math.PI / 16).toFloat()

        /** Constant twiddle factor  */
        val C_3_16 = Math.cos(2.0 * Math.PI * 3.0 / 16).toFloat()

        /** Constant twiddle factor  */
        val SQRT2BY2 = (Math.sqrt(2.0) / 2.0).toFloat()
    }

}
