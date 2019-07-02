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

import kotlin.math.*


/**
 * Package-private class implementing an arbitrary power-of-two length complex DFT with the split radix algorithm.

 * Creates smaller CDFTsr instances recursively and calls these in the evaluation.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
internal open class CDFTsr {

    /** float[] pair specifying complex input sequence.  */
    protected lateinit var xr: FloatArray
    protected lateinit var xi: FloatArray

    /** float[] pair specifying complex output transform.  */
    protected lateinit var outXr: FloatArray
    protected lateinit var outXi: FloatArray

    /** int specifying offset into the top-level length-N sequence arrays.  */
    protected var xoffset: Int = 0

    /** int specifying the stride of the butterflies at this level.  */
    protected var xstride: Int = 0

    /** int specifying the offset into the top-level length-N transform arrays.  */
    protected var outXoffset: Int = 0

    /** Log base 2 of the length of this DFT.  */
    protected var m: Int = 0

    /** int specifying the length of this transform.  */
    protected var N: Int = 0

    /** N/8.  */
    protected var Ndiv8: Int = 0

    /** N/4.  */
    protected var Ndiv4: Int = 0

    /** CDFTsr instances created recursively from this level.  */
    private var dft1: CDFTsr? = null
    private var dft2: CDFTsr? = null
    private var dft3: CDFTsr? = null

    /** float[] containing a reference to the cos(x) table.  */
    private var c: FloatArray? = null

    /** float[] containing a reference to the cos(3*x) table.  */
    private var c3: FloatArray? = null

    /** float[] containing a reference to the sin(x) table.  */
    private var s: FloatArray? = null

    /** float[] containing a reference to the sin(3*x) table.  */
    private var s3: FloatArray? = null

    /** int factor specifying stride in cos/sin tables.  */
    private var f: Int = 0

    /** int specifying offset into cos/sin tables.  */
    private var reflect: Int = 0


    /**
     * Instantiates a top-level CDFTsr.

     * @param m    int specifying the log base 2 of the length of this CDFTsr.
     * *
     * @param c    reference to the cos(x) table.
     * *
     * @param c3   reference to the cos(3*x) table.
     * *
     * @param s    reference to the sin(x) table.
     * *
     * @param s3   reference to the sin(3*x) table.
     */
    constructor(m: Int, c: FloatArray, c3: FloatArray, s: FloatArray, s3: FloatArray) {

        this.m = m
        N = 1 shl m
        Ndiv8 = N / 8
        Ndiv4 = N / 4
        xoffset = 0
        xstride = 1
        outXoffset = 0

        this.c = c
        this.c3 = c3
        this.s = s
        this.s3 = s3

        f = 1
        reflect = 2 * c.size

        if (m > 6) {
            dft1 = CDFTsr(this, 0, 2, 0, m - 1)
            dft2 = CDFTsr(this, 1, 4, N / 2, m - 2)
            dft3 = CDFTsr(this, 3, 4, 3 * N / 4, m - 2)
        } else if (m == 6) {
            dft1 = CDFTsr(this, 0, 2, 0, 5)
            dft2 = CDFTsr16(1, 4, N / 2)
            dft3 = CDFTsr16(3, 4, 3 * N / 4)
        } else if (m == 5) {
            dft1 = CDFTsr16(0, 2, 0)
            dft2 = CDFTsr8(1, 4, N / 2)
            dft3 = CDFTsr8(3, 4, 3 * N / 4)
        }

    }


    /**
     * Default constructor.
     */
    protected constructor() {
        dft1 = null
        dft2 = null
        dft3 = null
    }


    /**
     * Instantiates a new CDFTsr - called by a parent CDFTsr.

     * @param parent          CDFTsr reference to the parent
     * *
     * @param dataOffset      int specifying the offset into the top-level data arrays.
     * *
     * @param dataStride      int specifying the stride of butterflies in the top-level data arrays.
     * *
     * @param transformOffset int specifying the offset of this transform into the top-level transform arrays.
     * *
     * @param m               int specifying the log base 2 of the length of this DFT.
     */
    protected constructor(parent: CDFTsr, dataOffset: Int, dataStride: Int, transformOffset: Int, m: Int) {

        c = parent.c
        c3 = parent.c3
        s = parent.s
        s3 = parent.s3

        this.m = m
        N = 1 shl m
        Ndiv8 = N / 8
        Ndiv4 = N / 4
        this.xoffset = dataOffset
        this.xstride = dataStride
        this.outXoffset = transformOffset

        f = c!!.size / Ndiv8
        reflect = 2 * c!!.size

        if (m > 6) {
            dft1 = CDFTsr(this, dataOffset, dataStride * 2, transformOffset, m - 1)
            dft2 = CDFTsr(this, dataOffset + dataStride, dataStride * 4, transformOffset + N / 2, m - 2)
            dft3 = CDFTsr(this, dataOffset + 3 * dataStride, dataStride * 4, transformOffset + 3 * N / 4, m - 2)
        } else if (m == 6) {
            dft1 = CDFTsr(this, dataOffset, dataStride * 2, transformOffset, 5)
            dft2 = CDFTsr16(dataOffset + dataStride, dataStride * 4, transformOffset + N / 2)
            dft3 = CDFTsr16(dataOffset + 3 * dataStride, dataStride * 4, transformOffset + 3 * N / 4)
        } else if (m == 5) {
            dft1 = CDFTsr16(dataOffset, dataStride * 2, transformOffset)
            dft2 = CDFTsr8(dataOffset + dataStride, dataStride * 4, transformOffset + N / 2)
            dft3 = CDFTsr8(dataOffset + 3 * dataStride, dataStride * 4, transformOffset + 3 * N / 4)
        }

    }


    /**
     * Links the user-supplied input sequence and output transform arrays.  Propagates the links to child DFTs.

     * @param xr  float[] containing the input sequence real part.
     * *
     * @param xi  float[] containing the input sequence imaginary part.
     * *
     * @param Xr  float[] containing the output sequence real part.
     * *
     * @param Xi  float[] containing the output sequence imaginary part.
     */
    @JsName("link")
    internal open fun link(xr: FloatArray, xi: FloatArray, Xr: FloatArray, Xi: FloatArray) {
        this.xr = xr
        this.xi = xi
        this.outXr = Xr
        this.outXi = Xi
        dft1!!.link(xr, xi, Xr, Xi)
        dft2!!.link(xr, xi, Xr, Xi)
        dft3!!.link(xr, xi, Xr, Xi)
    }


    /**
     * Evaluates the complex DFT.
     */
    internal open fun evaluate() {

        var T1r: Float
        var T1i: Float
        var T3r: Float
        var T3i: Float
        var Rr: Float
        var Ri: Float
        var Sr: Float
        var Si: Float
        var Wr: Float
        var Wi: Float

        dft1!!.evaluate()
        dft2!!.evaluate()
        dft3!!.evaluate()

        // k = 0 butterfly

        var kp = outXoffset
        var kpN4 = kp + Ndiv4
        var kpN2 = kpN4 + Ndiv4
        var kp3N4 = kpN2 + Ndiv4

        Rr = outXr[kpN2] + outXr[kp3N4]
        Ri = outXi[kpN2] + outXi[kp3N4]
        Sr = outXi[kp3N4] - outXi[kpN2]
        Si = outXr[kpN2] - outXr[kp3N4]

        outXr[kpN2] = outXr[kp] - Rr
        outXi[kpN2] = outXi[kp] - Ri
        outXr[kp3N4] = outXr[kpN4] + Sr
        outXi[kp3N4] = outXi[kpN4] + Si

        outXr[kp] += Rr
        outXi[kp] += Ri
        outXr[kpN4] -= Sr
        outXi[kpN4] -= Si

        // k = 1 through Ndiv8-1 butterflies

        var fk: Int

        for (k in 1..Ndiv8 - 1) {

            fk = f * k
            kp = k + outXoffset
            kpN4 = kp + Ndiv4
            kpN2 = kpN4 + Ndiv4
            kp3N4 = kpN2 + Ndiv4

            // T1 = Wk*O1
            // T3 = W3k*O3

            Wr = c!![fk]
            Wi = s!![fk]
            T1r = Wr * outXr[kpN2] - Wi * outXi[kpN2]
            T1i = Wr * outXi[kpN2] + Wi * outXr[kpN2]
            Wr = c3!![fk]
            Wi = s3!![fk]
            T3r = Wr * outXr[kp3N4] - Wi * outXi[kp3N4]
            T3i = Wr * outXi[kp3N4] + Wi * outXr[kp3N4]

            // R = T1 + T3
            // S = i*(T1 - T3)

            Rr = T1r + T3r
            Ri = T1i + T3i
            Sr = T3i - T1i
            Si = T1r - T3r

            outXr[kpN2] = outXr[kp] - Rr
            outXi[kpN2] = outXi[kp] - Ri
            outXr[kp3N4] = outXr[kpN4] + Sr
            outXi[kp3N4] = outXi[kpN4] + Si

            outXr[kp] += Rr
            outXi[kp] += Ri
            outXr[kpN4] -= Sr
            outXi[kpN4] -= Si
        }

        // k = N/8 butterfly

        kp = Ndiv8 + outXoffset
        kpN4 = kp + Ndiv4
        kpN2 = kpN4 + Ndiv4
        kp3N4 = kpN2 + Ndiv4

        // T1 = Wk*O1
        // T3 = W3k*O3

        T1r = SQRT2BY2 * (outXr[kpN2] + outXi[kpN2])
        T1i = SQRT2BY2 * (outXi[kpN2] - outXr[kpN2])

        T3r = SQRT2BY2 * (outXi[kp3N4] - outXr[kp3N4])
        T3i = -SQRT2BY2 * (outXi[kp3N4] + outXr[kp3N4])

        // R = T1 + T3
        // S = i*(T1 - T3)

        Rr = T1r + T3r
        Ri = T1i + T3i
        Sr = T3i - T1i
        Si = T1r - T3r

        outXr[kpN2] = outXr[kp] - Rr
        outXi[kpN2] = outXi[kp] - Ri
        outXr[kp3N4] = outXr[kpN4] + Sr
        outXi[kp3N4] = outXi[kpN4] + Si

        outXr[kp] += Rr
        outXi[kp] += Ri
        outXr[kpN4] -= Sr
        outXi[kpN4] -= Si

        // k = N/8+1 through N/4-1 butterflies

        for (k in Ndiv8 + 1..Ndiv4 - 1) {

            fk = reflect - f * k
            kp = k + outXoffset
            kpN4 = kp + Ndiv4
            kpN2 = kpN4 + Ndiv4
            kp3N4 = kpN2 + Ndiv4

            // T1 = Wk*O1
            // T3 = W3k*O3

            Wr = -s!![fk]
            Wi = -c!![fk]
            T1r = Wr * outXr[kpN2] - Wi * outXi[kpN2]
            T1i = Wr * outXi[kpN2] + Wi * outXr[kpN2]
            Wr = s3!![fk]
            Wi = c3!![fk]
            T3r = Wr * outXr[kp3N4] - Wi * outXi[kp3N4]
            T3i = Wr * outXi[kp3N4] + Wi * outXr[kp3N4]

            // R = T1 + T3
            // S = i*(T1 - T3)

            Rr = T1r + T3r
            Ri = T1i + T3i
            Sr = T3i - T1i
            Si = T1r - T3r

            outXr[kpN2] = outXr[kp] - Rr
            outXi[kpN2] = outXi[kp] - Ri
            outXr[kp3N4] = outXr[kpN4] + Sr
            outXi[kp3N4] = outXi[kpN4] + Si

            outXr[kp] += Rr
            outXi[kp] += Ri
            outXr[kpN4] -= Sr
            outXi[kpN4] -= Si
        }

    }

    companion object {

        /** Constant twiddle factor for N/2 butterfly.  */
        private val SQRT2BY2 = (sqrt(2.0) / 2.0).toFloat()
    }

}
