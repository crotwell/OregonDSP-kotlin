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

package com.oregondsp.signalProcessing.filter.iir;

import java.io.PrintStream;
import java.text.DecimalFormat;

/**
 * Class to implement a second order section - basic unit of an Infinite Impulse Response digital filter.
 * 
 * <p>Implements the finite difference equation:</p>
 * 
 * <p>y[n] = -a[1]*y[n-1] - a[2]*y[n-2] + b[0]*x[n] + b[1]*x[n-1] + b[2]*x[n-2]</p>
 * 
 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
public class SecondOrderSection {
	
	/** Numerator coefficients */
	double b0, b1, b2;
	
	/** Denominator coefficients (a0 = 1) by assumption. */
	double a1, a2;
	
	/** States required to support processing of a continuous data stream in consecutive, contiguous blocks. */
	double s1, s2;
	
	
	/**
	 * Instantiates a new second order section, with values for the numerator and denominator coefficients.
	 *
	 * @param b0         Numerator coefficient b[0].
	 * @param b1         Numerator coefficient b[1].
	 * @param b2         Numerator coefficient b[2].
	 * @param a1         Denominator coefficient a[1].
	 * @param a2         Denominator coefficient a[2].
	 */
	public SecondOrderSection( double b0, double b1, double b2, double a1, double a2 ) {
		
		this.b0 = b0;
		this.b1 = b1;
		this.b2 = b2;
		this.a1 = a1;
		this.a2 = a2;
		
		initialize();
	}
	
	
	/**
	 * Initializes states to zero.
	 */
	public void initialize() {
	  s1 = 0.0;
	  s2 = 0.0;
	}
	
	
	/**
	 * Filters a single input sample (single-step filtering).
	 *
	 * @param x     float containing value of the single input sample.
	 * @return      float result of the filter for one time step.
	 */
	public float filter( float x ) {
		
		double s0    = x - a1*s1 - a2*s2;
		float retval = (float) (b0*s0 + b1*s1 + b2*s2);
		
		s2 = s1;
		s1 = s0;
		
		return retval;
	}

	
	/**
	 * Filters a sequence of input samples.
	 *
	 * @param x     float[] containing the sequence of input samples.
	 * @param y     float[] containing the filtered result.  May be the same array as x.
	 */
	public void filter( float[] x, float[] y ) {
		
		double s0;
		
		int n = Math.min( x.length, y.length );
		
		for ( int i = 0;  i < n;  i++ ) {
			s0   = x[i] - a1*s1 -a2*s2;
			y[i] = (float) ( b0*s0 + b1*s1 + b2*s2 );
			s2 = s1;
			s1 = s0;
		}
	}
	
	
	
	/**
	 * Prints the filter coefficients and states.
	 *
	 * @param ps       PrintStream to which the filter coefficients and states are printed.
	 */
	public void print( PrintStream ps ) {
		
		DecimalFormat formatter = new DecimalFormat( "##0.00000" );
		
		ps.println( "  coefficients: \n" );
		ps.println( "    b0: " + formatter.format(b0) );
		ps.println( "    b1: " + formatter.format(b1) );
		ps.println( "    b2: " + formatter.format(b2) );
		ps.println();
		ps.println( "    a1: " + formatter.format(a1) );
		ps.println( "    a2: " + formatter.format(a2) );
		ps.println( "\n  states:  \n" );
		ps.println( "    s1: " + formatter.format(s1) );
		ps.println( "    s2: " + formatter.format(s2) );
	}
	
}
