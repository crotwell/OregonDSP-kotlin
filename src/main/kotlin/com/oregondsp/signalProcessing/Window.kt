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


/**
 * Base class for implementing windows - partial implementation.

 * This class and its derived classes make it possible to "snip out" a portion of a signal
 * beginning at a specified index and shape the resulting segment with a specified weighting
 * function.  The length of the resulting segment is the same length as the window.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
open class Window {

    /** float[] containing the window coefficients.  */
    protected var w: FloatArray


    /**
     * Instantiates a new Window from a vector of coefficients.

     * @param w     float[] containin the vector of window coefficients.
     */
    constructor(w: FloatArray) {
        this.w = w.copyOf()
    }


    /**
     * Instantiates a new length-N window containing zeros.

     * @param N     int specifying the window length in samples.
     */
    constructor(N: Int) {
        w = FloatArray(N)
    }


    /**
     * Returns the length of the window in samples.

     * @return    int containing the window length in samples.
     */
    fun length(): Int {
        return w.size
    }


    /**
     * Allows a window to be modified in-place by multiplication by another window.

     * @param x    float[] containing the coefficients of the second window, which modifies the first (this) window.
     */
    fun timesEquals(x: FloatArray) {
        if (x.size != w.size) throw IllegalArgumentException("Argument length does not match window length")
        for (i in w.indices) w[i] *= x[i]
    }


    /**
     * Returns a copy of the coefficients of this window.

     * @return     float[] containing window coefficients.
     */
    val array: FloatArray
        get() = w.copyOf()


    /**
     * Windows a sequence and places the result in a specified array.

     * @param x          float[] containing the sequence to be windowed by this Window.
     * *
     * @param index      start point in the input sequence at which this Window is applied.
     * *
     * @param y          float[] containing the resulting windowed sequence.
     */
    fun window(x: FloatArray, index: Int, y: FloatArray) {

        if (y.size != w.size) throw IllegalArgumentException("Destination array length does not match window length")

        for (i in w.indices) {
            val j = index + i
            if (j >= 0 && j < x.size)
                y[i] = w[i] * x[j]
            else
                y[i] = 0.0f
        }

    }

}
