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

package com.oregondsp.signalProcessing.filter.iir

import java.io.FileNotFoundException
import java.io.FileOutputStream
import java.io.PrintStream
import kotlin.js.Math


/**
 * Class implementing Chebyshev Type II filters - characterized by zeros in the stop band.

 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
class ChebyshevII
/**
 * Instantiates a new Chebyshev type II filter.

 * @param order      int specifying the order (number of poles) of the filter.
 * *
 * @param epsilon    double design parameter specifying the passband ripple and stopband attenuation.
 * *
 * @param type       PassbandType specifying whether the filter is a lowpass, bandpass or highpass filter.
 * *
 * @param f1         double specifying the low cutoff frequency (must always be present, but used only for
 * *                   bandpass and highpass filters).
 * *
 * @param f2         double specifying the high cutoff frequency (must always be present, but used only for
 * *                   bandpass and lowpass filters).
 * *
 * @param delta      double specifying the sampling interval of the data to which this filter will be applied.
 */
(order: Int, epsilon: Double, type: PassbandType, f1: Double, f2: Double, delta: Double) : IIRFilter(AnalogChebyshevII(order, epsilon), type, f1, f2, delta) {
    companion object {


        @JvmStatic fun main(args: Array<String>) {

            val epsilon = 0.01
            val F = ChebyshevII(7, epsilon, PassbandType.LOWPASS, 2.0, 0.0, 0.05)
            F.print(System.out)
            val tmp = FloatArray(201)
            for (i in 0..200) {
                val C = F.evaluate(Math.PI / 200.0 * i)
                tmp[i] = Complex.abs(C).toFloat()
            }

            val x = FloatArray(1001)
            x[200] = 1.0f
            val y = FloatArray(1001)
            F.filter(x, y)

            var ps: PrintStream
            try {
                ps = PrintStream(FileOutputStream("C:\\DATA\\Test\\Response.m"))
                ps.print("R = [ ")
                for (i in 0..199) {
                    ps.println(tmp[i].toString() + "; ...")
                }
                ps.println(tmp[200].toString() + "];")
                ps.close()
            } catch (e: FileNotFoundException) {
                e.printStackTrace()
            }

            try {
                ps = PrintStream(FileOutputStream("C:\\DATA\\Test\\ImpulseResponse.m"))
                ps.print("IR = [ ")
                for (i in 0..999) {
                    ps.println(y[i].toString() + "; ...")
                }
                ps.println(y[1000].toString() + "];")
                ps.close()
            } catch (e: FileNotFoundException) {
                e.printStackTrace()
            }

        }
    }


}
