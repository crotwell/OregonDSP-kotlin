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
 * Designs and implements Hamming windows (see Oppenheim and Schafer, 1975).

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
class HammingWindow
/**
 * Instantiates a new Hamming window of length N samples.

 * @param N       int specifying the window length.
 */
(N: Int) : Window(N) {


    init {

        for (i in 0..N - 1) {
            w[i] = (0.54 + 0.46 * Math.cos(-Math.PI + i.toDouble() * 2.0 * Math.PI / (N - 1))).toFloat()
        }

    }

}
