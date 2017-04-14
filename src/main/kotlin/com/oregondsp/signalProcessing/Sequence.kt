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

package com.oregondsp.signalProcessing

import kotlin.js.Math


/**
 * Class to implement basic signal processing operations on scalar sequences used by other classes in this package.

 *
 * This class can be used in two ways.  Objects of this class may be instantiated to represent sequences and in that
 * case will serve as a containers for sequence values (float precision).  For this use, methods are supplied
 * that alter or operate on the values of the contained sequence.  Alternatively, the methods (where it makes
 * sense) are supplied in static form.  The class may be used to operate on the user's float arrays without
 * need to instantiate a Sequence object.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
class Sequence {


    /** float[] containing sequence values  */
    /**
     * Accessor for the entire array of sequence values.  Allows the sequence to be modified.

     * @return      float[] reference to the sequence.  Does not return a copy.
     */
    var array: FloatArray
        protected set


    /**
     * Instantiates a new sequence from an array of float values.

     * @param x    float[] containing the Sequence values.
     */
    constructor(x: FloatArray) {
        this.array = FloatArray(x.size)
        //System.arraycopy(x, 0, this.array, 0, x.size)
        this.array = x.copyOf();
    }


    /**
     * Instantiates a new sequence of all zeros, of length N samples.

     * @param N    int specifying the length of the sequence.
     */
    constructor(N: Int) {
        array = FloatArray(N)
    }


    /**
     * Aliases the current sequence into a smaller sequence.  Modifies this Sequence.

     * @param N    New length of the sequence, and alias modulus.
     */
    fun alias(N: Int) {
        val newx = FloatArray(N)
        alias(array, newx)
        array = newx
    }


    /**
     * Accessor for an individual value of the sequence.

     * @param index     int containing the index of the desired sequence value.
     * *
     * @return          Sequence value at index.
     */
    operator fun get(index: Int): Float {
        var retval = 0.0f
        if (index >= 0 && index < array.size) retval = array[index]
        return retval
    }


    /**
     * Reverses this sequence in-place.
     */
    fun reverse() {
        reverse(array)
    }


    /**
     * Removes the mean of this sequence in-place.
     */
    fun rmean() {
        rmean(array)
    }


    /**
     * Performs a circular shift on this sequence, in-place

     * @param shift     int specifying the size and direction of the shift.  A negative number
     * *                  specifies a left shift and a positive number, a right shift.  A zero shift
     * *                  leaves the sequence unchanged.
     */
    fun circularShift(shift: Int) {
        circularShift(array, shift)
    }


    /**
     * Performs a shift on this sequence with zero-fill.

     * @param shift   int specifying the direction and size of the shift.  A negative shift is to the left.
     * *                Zeros are fed in from the right in that case.  A positive shift is to the right.  Zeros
     * *                are fed in from the left in that case.  A zero shift leaves the sequence unchanged.
     */
    fun zeroShift(shift: Int) {
        zeroShift(array, shift)
    }


    /**
     * Decimates this sequence in-place.

     * @param decrate   int specifying the decimation rate.
     */
    fun decimate(decrate: Int) {
        val tmp = FloatArray(array.size / decrate)
        decimate(array, decrate, tmp)
        array = tmp
    }


    /**
     * Stretches this sequence by a specified rate, in place.

     * This operation spreads the sequence values out and fills between them with interstitial zeros.
     * It is a basic operation needed for interpolation by an integer rate.

     * @param rate     int containing the stretch rate (factor).
     */
    fun stretch(rate: Int) {
        val tmp = FloatArray(array.size * rate)
        stretch(array, rate, tmp)
        array = tmp
    }


    /**
     * Multiplies this sequence by a constant, in-place.

     * @param f   float specifying the multiplicative constant.
     */
    fun timesEquals(f: Float) {
        timesEquals(array, f)
    }


    /**
     * Pads this sequence to length n, by zero filling on right if n > length of this sequence, no-op otherwise.

     * @param n           int specifying desired length of padded sequence
     */
    fun pad(n: Int) {
        if (n > array.size) {
            val tmp = FloatArray(n)
            pad(array, tmp)
            array = tmp
        }
    }

    companion object {


        /**
         * Method to alias a source sequence into a destination sequence.

         * Source value src[n] is added to dst[ mod(n,dst.length) ].

         * @param src     float[] containing the source sequence.
         * *
         * @param dst     float[] containing the destination sequence.
         */
        @JsName("aliasArray")
        fun alias(src: FloatArray, dst: FloatArray) {

            val slength = src.size
            val dlength = dst.size

            //Arrays.fill(dst, 0.0f)
            for (i in 0..dlength)
                dst[i] = 0.0F

            for (i in 0..slength - 1)
                dst[i % dlength] += src[i]

        }


        /**
         * Reverses a sequence in place.

         * @param  y   float[] containing the sequence to be reversed.
         */
        @JsName("reverseArray")
        fun reverse(y: FloatArray) {
            var i = 0
            var j = y.size - 1
            while (i < j) {
                val tmp = y[i]
                y[i] = y[j]
                y[j] = tmp
                i++
                j--
            }
        }


        /**
         * Removes the mean of a sequence.

         * @param y    float[] specifying the sequence to be demeaned.
         */
        @JsName("rmeanArray")
        fun rmean(y: FloatArray) {
            var mean = 0.0f
            for (i in y.indices) mean += y[i]
            mean /= y.size.toFloat()
            for (i in y.indices) y[i] -= mean
        }


        /**
         * Performs a circular shift of a sequence.

         * @param y       float[] containing sequence to be shifted.
         * *
         * @param shift   int specifying the size and direction of the shift.  A negative
         * *                number specifies a left shift, a positive number, a right shift.
         * *                A zero shift leaves the sequence unchanged.
         */
        @JsName("circularShiftArray")
        fun circularShift(y: FloatArray, shift: Int) {

            val N = y.size
            var s = shift % N

            // minimize shift - consider alternative shift in opposite direction

            if (s > 0 && N - s < s)
                s -= N
            else if (s < 0 && N + s < -s)
                s += N

            // right shift
            if (s<0) s *= -1
            val tmp = FloatArray(s)

            if (s > 0) {
                for (i in 0..s - 1)
                    tmp[i] = y[N - s + i]
                for (i in N - 1 - s downTo 0)
                    y[i + s] = y[i]
                for (i in 0..s - 1)
                    y[i] = tmp[i]
            }

            // left shift

            if (s < 0) {
                for (i in 0..-s - 1)
                    tmp[i] = y[i]
                for (i in -s..N - 1)
                    y[i + s] = y[i]
                for (i in 0..-s - 1)
                    y[N + s + i] = tmp[i]
            }


        }


        /**
         * Shifts a sequence left or right and pads with zeros (unlike the circular shift, sequence values are lost).

         * @param y       float[] containing the sequence to be shifted.
         * *
         * @param shift   int specifying the direction and size of the shift.  A negative shift is to the left.
         * *                Zeros are fed in from the right in that case.  A positive shift is to the right.  Zeros
         * *                are fed in from the left in that case.  A zero shift leaves the sequence unchanged.
         */
        @JsName("zeroShiftArray")
        fun zeroShift(y: FloatArray, shift: Int) {

            //if (Math.abs(shift) >= y.size)
            if (-1*shift >= y.size || shift >= y.size)
                //Arrays.fill(y, 0.0f)
                for (i in 0..y.size)
                    y[i] = 0.0F
            else if (shift > 0) {
                for (i in y.size - 1 downTo shift)
                    y[i] = y[i - shift]
                for (i in 0..shift - 1)
                    y[i] = 0.0f
            } else if (shift < 0) {
                for (i in 0..y.size + shift - 1)
                    y[i] = y[i - shift]
                for (i in y.size + shift..y.size - 1)
                    y[i] = 0.0f
            }

        }


        /**
         * Decimates a sequence by a specified rate.

         * @param y            float[] containing the sequence to be decimated.
         * *
         * @param decrate      int specifying the decimation rate.
         * *
         * @param ydecimated   float[] containing the decimated sequence.
         */

        @JsName("decimateArray")
        fun decimate(y: FloatArray, decrate: Int, ydecimated: FloatArray) {
            val n = Math.min(ydecimated.size, y.size / decrate)
            for (i in 0..n - 1) ydecimated[i] = y[i * decrate]
        }


        /**
         * Stretches a sequence by a specified rate.

         * This operation spreads the sequence values out and fills between them with interstitial zeros.
         * It is a basic operation needed for interpolation by an integer rate.

         * @param y            float[] containing the sequence to be stretched.
         * *
         * @param rate         int containing the stretch rate.
         * *
         * @param ystretched   float[] containing the stretched sequence.
         */

        @JsName("stretchArray")
        fun stretch(y: FloatArray, rate: Int, ystretched: FloatArray) {
            val n = Math.min(y.size, ystretched.size / rate)
            //Arrays.fill(ystretched, 0f)
            for (i in 0..ystretched.size)
                ystretched[i] = 0.0F
            for (i in 0..n - 1) ystretched[i * rate] = y[i]
        }


        /**
         * Multiplies a sequence by a constant.

         * @param y       float[] containing the sequence to be scaled.
         * *
         * @param f       float containing the multiplicative constant.
         */
        @JsName("timesEqualsArray")
        fun timesEquals(y: FloatArray, f: Float) {
            for (i in y.indices) y[i] *= f
        }


        /**
         * Pad a sequence with zeros (on the right)

         * If ypadded.length < y.length, this method performs a truncation of y.

         * @param y           float[] containing original sequence
         * *
         * @param ypadded     float[] containing padded sequence
         */
        @JsName("padArray")
        fun pad(y: FloatArray, ypadded: FloatArray) {
            if (y.size < ypadded.size) {
                //Arrays.fill(ypadded, 0.0f)
                //System.arraycopy(y, 0, ypadded, 0, y.size)
                for (i in 0..y.size)
                    ypadded[i] = y[i]
                for (i in y.size..ypadded.size)
                    ypadded[i] = 0.0F

            } else {
                //System.arraycopy(y, 0, ypadded, 0, ypadded.size)
                for (i in 0..ypadded.size)
                    ypadded[i] = y[i]
            }
        }
    }

}
