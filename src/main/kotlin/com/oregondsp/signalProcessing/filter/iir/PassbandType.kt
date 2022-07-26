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


/**
 * Enum PassbandType used to specify the pass band type (lowpass, highpass, bandpass) of analog and digital filters.

 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
@JsExport
enum class PassbandType {

    /** Specifies a lowpass filter  */
    LOWPASS,
    /** Specifies a bandpass filter  */
    BANDPASS,
    /** Specifies a highpass filter  */
    HIGHPASS
}
