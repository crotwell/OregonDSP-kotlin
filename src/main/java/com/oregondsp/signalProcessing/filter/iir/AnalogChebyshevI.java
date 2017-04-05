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

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;

import com.oregondsp.signalProcessing.filter.Polynomial;
import com.oregondsp.signalProcessing.filter.Rational;

/**
 * Class to design analog Chebyshev type I prototype filters.
 * 
 * This analog prototype is a lowpass filter with a cutoff at 1 radian per second.  Chebyshev
 * type I filters have ripples in the passband and a sharper transition from passband to 
 * stopband than a Butterworth filter.
 * 
 * @author David B. Harris,  Deschutes Signal Processing LLC
 */
public class AnalogChebyshevI extends AnalogPrototype {
  
  /**
   * Instantiates a new analog Chebyshev type I prototype filter.
   *
   * @param order    int specifying the number of poles of the filter.
   * @param epsilon  double parameter controlling the stopband attenuation and passband ripple. 
   */
  public AnalogChebyshevI ( int order, double epsilon ) {
    
    super();
    
    double alpha = ( 1.0 + Math.sqrt( 1.0 + epsilon*epsilon ) ) / epsilon ;
    double p     = Math.pow( alpha, 1.0/order );
    double a     = 0.5*( p - 1/p );
    double b     = 0.5*( p + 1/p );
    
    System.out.println( "alpha: " + alpha );
    System.out.println( "p:     " + p );
    System.out.println( "a:     " + a );
    System.out.println( "b:     " + b );
    
    int nRealPoles        = order - 2*(order/2);
    int nComplexPolePairs = order/2;
    int nPoles            = nRealPoles + 2*nComplexPolePairs;
    
    if ( nRealPoles == 1 ) {
      double[] td = { a, 1.0};
      addSection( new Rational( new Polynomial(1.0), new Polynomial(td) ) );
    }
    
    double dAngle = Math.PI/nPoles;

    for ( int i = 0;  i < nComplexPolePairs;  i++ ) {
        double   angle = -Math.PI/2  +  dAngle/2 *( 1 + nRealPoles )  +  i*dAngle;
        Complex pole = new Complex( a*Math.sin(angle), b*Math.cos(angle ) );
        double[] td    = { pole.real()*pole.real() + pole.imag()*pole.imag(), -2*pole.real(), 1.0 };
        addSection( new Rational( new Polynomial(1.0), new Polynomial(td) ) );
    }
    
    // scale to 1 at s = 0
    
    sections.get( 0 ).timesEquals( 1.0 / ( Math.pow(2.0, order-1)*epsilon ) );
    
  }
  
  
  
  public static void main( String[] args ) {
    
    AnalogChebyshevI A = new AnalogChebyshevI( 4, 0.50885 );
    AnalogPrototype B = A.lptolp( 0.2*Math.PI );

    float[] tmp = new float[201];
    for ( int i = 0;  i < 201;  i++ ) {
      Complex C = B.evaluate( i*0.02 );
      tmp[i] = (float) Complex.abs( C );
    }
    
    PrintStream ps;
    try {
      ps = new PrintStream( new FileOutputStream( "C:\\DATA\\Test\\AnalogResponse.m" ) );
      ps.print( "R = [ " );
      for ( int i = 0;  i < 200;  i++ ) {
          ps.println( tmp[i] + "; ..." );
      }
      ps.println( tmp[200] + "];" );
      ps.close();
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }

  }

}
