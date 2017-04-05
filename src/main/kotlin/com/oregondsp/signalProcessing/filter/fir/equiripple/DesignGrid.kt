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

package com.oregondsp.signalProcessing.filter.fir.equiripple

import java.io.PrintStream
import java.text.DecimalFormat


/**
 * Contains the finite frequency-sampling grid used by the Remez exchange algorithm.

 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
internal class DesignGrid {

    /** double[] containing the grid samples.  */
    var grid: DoubleArray? = null                 // sample points where designs are evaluated

    /** int specifying the grid size.  */
    var gridSize: Int = 0

    /** double array containing the transformed grid points.  */
    var X: DoubleArray? = null

    /** double[] containing the desired (weighted) response of the filter on the grid.  */
    var H: DoubleArray? = null                    // desired response of filter on grid

    /** double[] containing the weighting function on grid.  */
    var W: DoubleArray? = null                    // weight function on grid

    /** int[] specifying indices of grid points that are band edges.  */
    var bandEdgeIndices: IntArray? = null

    /** int[] specifying indices of grid points that are current extrema in the Remez exchange.  */
    var extremaIndices: IntArray? = null

    /** boolean value specifying whether the grid contains 0 frequency as a sample.  */
    var containsZero: Boolean = false

    /** boolean value specifying whether the grid contains frequency pi (1.0 normalized) as a sample.  */
    var containsPi: Boolean = false


    /**
     * Prints the grid to a PrintStream instance - useful for debugging.

     * @param ps    PrintStream instance to which the grid is printed.
     */
    fun print(ps: PrintStream) {

        val F = DecimalFormat("0.000000")
        val I = DecimalFormat("000")
        var extremum = 0
        var bandEdgeCount = 0
        for (i in 0..gridSize - 1) {
            val Omega = grid!![i]
            var line = I.format(i.toLong()) + "  " + F.format(Omega) + "  " +
                    F.format(X!![i]) + "  " +
                    F.format(H!![i]) + "  " +
                    F.format(W!![i])
            if (bandEdgeIndices!![bandEdgeCount] == i) {
                line = line + "  band edge"
                bandEdgeCount++
            }
            if (Omega == grid!![extremaIndices!![extremum]]) {
                line = line + "  extremum"
                extremum++
            }

            ps.println(line)
        }

    }

    companion object {

        /** Constant GRIDDENSITY partly specifies the number of grid points (along with N).  */
        val GRIDDENSITY = 20
    }

}
