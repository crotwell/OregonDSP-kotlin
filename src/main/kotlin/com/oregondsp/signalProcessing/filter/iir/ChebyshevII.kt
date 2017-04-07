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
    


}
