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

package com.oregondsp.signalProcessing.filter.fir

import com.oregondsp.signalProcessing.HammingWindow
import com.oregondsp.signalProcessing.Sequence
import kotlin.js.Math


/**
 * Implements a finite impulse response (FIR) interpolator for resampling sequences up by an integer rate.

 *
 * This implementation uses a windowed sinc function as the interpolation kernel.  The length of the window
 * is determined by the interpolation rate and a design factor which controls interpolation accuracy.  A
 * Hamming window is used to window the sinc function.

 *
 * The window length is N = 2*M + 1, where M = (interpolation rate)*designFactor.  The interpolation kernel
 * is:
 *
 *
 * h[n] = ( 0.54 + 0.46*cos(pi/M * n) ) * sin( n*pi/R ) / ( n*pi/R )
 *
 *
 * R is the interpolation rate.

 *
 * This implementation assumes that the data will be processed in consecutive, contiguous blocks of uniform
 * size.  The block size is specified in the constructor for the Interpolator, and thereafter calls to interpolate()
 * will interpolate successive blocks of that size.  Very important:  the block size must not change.  State
 * information required to assure continuous interpolation from block to block is maintained by the class instance,
 * so that there are no edge effects in the interpolated result sequence (except at the very beginning and end of
 * the data stream).

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
class Interpolator
/**
 * Instantiates a new Interpolator.

 * @param rate                int containing the interpolation rate.
 * *
 * @param designFactor        int containing parameter controlling the accuracy of interpolation
 * *                            through the length of the interpolation filter
 * *
 * @param blockSize           int controlling the size of blocks in which the data will be processed
 */
(
        /** int containing the interpolation rate.  */
        private val rate: Int, designFactor: Int, blockSize: Int) {

    /** OverlapAdd instance that performs the interpolation by convolution with a windowed sinc function.  */
    private val overlapAdd: OverlapAdd

    /** float[] that buffers the low-rate sequence to be interpolated.  */
    private val buffer: FloatArray


    init {

        val half = rate * designFactor
        val N = 2 * half + 1
        val kernel = HammingWindow(N).array
        for (i in 1..half) {
            kernel[half + i] *= (Math.sin(Math.PI * i / rate) / (Math.PI * i / rate)).toFloat()
            kernel[half - i] = kernel[half + i]
        }

        overlapAdd = OverlapAdd(kernel, blockSize * rate)

        buffer = FloatArray(blockSize * rate)
    }


    /**
     * Interpolates blocks of a data stream.



     * @param block                 float[] containing a block of the sequence to be interpolated
     * *
     * @param interpolatedBlock     float[] containing the interpolated block
     */
    fun interpolate(block: FloatArray, interpolatedBlock: FloatArray) {
        Sequence.stretch(block, rate, buffer)
        overlapAdd.filter(buffer, 0, interpolatedBlock, 0)
    }

}
