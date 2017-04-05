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

package com.oregondsp.signalProcessing;

import java.util.Arrays;


/**
 * Class to implement basic signal processing operations on scalar sequences used by other classes in this package.
 * 
 * <p>This class can be used in two ways.  Objects of this class may be instantiated to represent sequences and in that
 * case will serve as a containers for sequence values (float precision).  For this use, methods are supplied
 * that alter or operate on the values of the contained sequence.  Alternatively, the methods (where it makes
 * sense) are supplied in static form.  The class may be used to operate on the user's float arrays without
 * need to instantiate a Sequence object.</p>
 * 
 * @author David B. Harris,   Deschutes Signal Processing LLC
 */
public class Sequence {
  
  
  /** float[] containing sequence values */
  protected float[] x;
  
  
  /**
   * Instantiates a new sequence from an array of float values.
   *
   * @param x    float[] containing the Sequence values.
   */
  public Sequence( float[] x ) {
    this.x = new float[ x.length ];
    System.arraycopy( x, 0, this.x, 0, x.length );
  }
  
  
  
  /**
   * Instantiates a new sequence of all zeros, of length N samples.
   *
   * @param N    int specifying the length of the sequence.
   */
  public Sequence( int N ) {
    x = new float[ N ];
  }
  
  
  
  /**
   * Method to alias a source sequence into a destination sequence.
   * 
   * Source value src[n] is added to dst[ mod(n,dst.length) ].
   *
   * @param src     float[] containing the source sequence.
   * @param dst     float[] containing the destination sequence.
   */
  public static void alias( float[] src, float[] dst ) {
    
    int slength = src.length;
    int dlength = dst.length;
    
    Arrays.fill( dst, 0.0f );
    
    for ( int i = 0;  i < slength;  i++ ) 
      dst[ i % dlength ] += src[i];
    
  }
  
  
  
  /**
   * Aliases the current sequence into a smaller sequence.  Modifies this Sequence.
   *
   * @param N    New length of the sequence, and alias modulus.
   */
  public void alias( int N ) {
    float[] newx = new float[N];
    alias( x, newx );
    x = newx;
  }
  
  
  
  /**
   * Accessor for an individual value of the sequence.
   *
   * @param index     int containing the index of the desired sequence value.
   * @return          Sequence value at index.
   */
  public float get( int index ) {
    float retval = 0.0f;
    if ( index >= 0  &&  index < x.length ) retval = x[index];
    return retval;
  }
  
  
  
  /**
   * Accessor for the entire array of sequence values.  Allows the sequence to be modified.
   *
   * @return      float[] reference to the sequence.  Does not return a copy.
   */
  public float[] getArray() {
    return x;
  }
  
  
  
  /**
   * Reverses a sequence in place.
   *
   * @param  y   float[] containing the sequence to be reversed.
   */
  public static void reverse( float[] y ) {
    int i = 0;  
    int j = y.length-1;
    while ( i < j ) {
      float tmp = y[i];
      y[i] = y[j];
      y[j] = tmp;
      i++;
      j--;
    }
  }
  
  
  
  /**
   * Reverses this sequence in-place.
   */
  public void reverse() {
    reverse(x);
  }
  
  
  
  /**
   * Removes the mean of a sequence.
   *
   * @param y    float[] specifying the sequence to be demeaned.
   */
  public static void rmean( float[] y ) {
    float mean = 0.0f;
    for ( int i = 0;  i < y.length;  i++ ) mean += y[i];
    mean /= y.length;
    for ( int i = 0;  i < y.length;  i++ ) y[i] -= mean;
  }
  
  
  
  /**
   * Removes the mean of this sequence in-place.
   */
  public void rmean() {
    rmean(x);
  }
  
  
  
  /**
   * Performs a circular shift of a sequence.
   *
   * @param y       float[] containing sequence to be shifted.
   * @param shift   int specifying the size and direction of the shift.  A negative
   *                number specifies a left shift, a positive number, a right shift.
   *                A zero shift leaves the sequence unchanged.
   */
  public static void circularShift( float[] y, int shift ) {
    
    int N = y.length;
    int s = shift % N;
    
    // minimize shift - consider alternative shift in opposite direction
    
    if ( s > 0  &&  N-s < s ) 
      s -= N;
    else if ( s < 0  &&  N+s < -s ) 
      s += N;
    
    // right shift
    
    float[] tmp = new float[ Math.abs( s ) ];

    if ( s > 0 ) {
      for ( int i = 0;  i < s;  i++ ) 
        tmp[i] = y[N-s+i];
      for ( int i = N-1-s;  i >= 0;  i-- ) 
        y[i+s] = y[i];
      for ( int i = 0;  i < s;  i++ ) 
        y[i] = tmp[i];
    }
    
    // left shift
    
    if ( s < 0 ) {
      for ( int i = 0;  i < -s;  i++ )
        tmp[i] = y[i];
      for ( int i = -s;  i < N;  i++ )
        y[i+s] = y[i];
      for ( int i = 0;  i < -s;  i++ ) 
        y[N+s+i] = tmp[i];
    }
    
    
  }
  
  
  
  /**
   * Performs a circular shift on this sequence, in-place
   *
   * @param shift     int specifying the size and direction of the shift.  A negative number
   *                  specifies a left shift and a positive number, a right shift.  A zero shift
   *                  leaves the sequence unchanged.
   */
  public void circularShift( int shift ) {
    circularShift( x, shift );
  }
  
  
  
  /**
   * Shifts a sequence left or right and pads with zeros (unlike the circular shift, sequence values are lost).
   *
   * @param y       float[] containing the sequence to be shifted.
   * @param shift   int specifying the direction and size of the shift.  A negative shift is to the left.
   *                Zeros are fed in from the right in that case.  A positive shift is to the right.  Zeros
   *                are fed in from the left in that case.  A zero shift leaves the sequence unchanged.
   */
  public static void zeroShift( float[] y, int shift ) {
    
    if ( Math.abs( shift) >= y.length )
      Arrays.fill( y, 0.0f );
    
    else if ( shift > 0 ) {
      for ( int i = y.length-1;  i >= shift;  i-- )
        y[i] = y[i-shift];
      for ( int i = 0;  i < shift;  i++ ) 
        y[i] = 0.0f;
    }
    
    else if ( shift < 0 ) {
      for ( int i = 0;  i < y.length+shift;  i++ ) 
        y[i] = y[i-shift];
      for ( int i = y.length+shift;  i < y.length;  i++ ) 
        y[i] = 0.0f;
    }
    
  }
  
  
  
  /**
   * Performs a shift on this sequence with zero-fill.
   *
   * @param shift   int specifying the direction and size of the shift.  A negative shift is to the left.
   *                Zeros are fed in from the right in that case.  A positive shift is to the right.  Zeros
   *                are fed in from the left in that case.  A zero shift leaves the sequence unchanged.
   */
  public void zeroShift( int shift ) {
    zeroShift( x, shift );
  }
  
  
  
  /**
   * Decimates a sequence by a specified rate.
   *
   * @param y            float[] containing the sequence to be decimated.
   * @param decrate      int specifying the decimation rate.
   * @param ydecimated   float[] containing the decimated sequence.  
   */
  public static void decimate( float[] y, int decrate, float[] ydecimated ) {
    int n = Math.min( ydecimated.length, y.length/decrate );
    for ( int i = 0;  i < n;  i++ ) ydecimated[i] = y[i*decrate];
  }
  
  
  
  /**
   * Decimates this sequence in-place.
   *
   * @param decrate   int specifying the decimation rate.
   */
  public void decimate( int decrate ) {
    float[] tmp = new float[ x.length/decrate ];
    decimate( x, decrate, tmp );
    x = tmp;
  }
  
  
  
  /**
   * Stretches a sequence by a specified rate.
   * 
   * This operation spreads the sequence values out and fills between them with interstitial zeros.
   * It is a basic operation needed for interpolation by an integer rate.
   *
   * @param y            float[] containing the sequence to be stretched.
   * @param rate         int containing the stretch rate.
   * @param ystretched   float[] containing the stretched sequence.
   */
  public static void stretch( float[] y, int rate, float[] ystretched ) {
    int n = Math.min( y.length, ystretched.length/rate );
    Arrays.fill( ystretched, 0.f );
    for ( int i = 0;  i < n;  i++ ) ystretched[i*rate] = y[i];
  }
  
  
  
  /**
   * Stretches this sequence by a specified rate, in place.
   *
   * This operation spreads the sequence values out and fills between them with interstitial zeros.
   * It is a basic operation needed for interpolation by an integer rate.
   * 
   * @param rate     int containing the stretch rate (factor).
   */
  public void stretch( int rate ) {
    float[] tmp = new float[ x.length*rate ];
    stretch( x, rate, tmp );
    x = tmp;
  }
  
  

  /**
   * Multiplies a sequence by a constant.
   *
   * @param y       float[] containing the sequence to be scaled.
   * @param f       float containing the multiplicative constant.
   */
  public static void timesEquals( float[] y, float f ) {
    for ( int i = 0;  i < y.length;  i++ ) y[i] *= f;
  }
  
  
  
  /**
   * Multiplies this sequence by a constant, in-place.
   *
   * @param f   float specifying the multiplicative constant.
   */
  public void timesEquals( float f ) {
    timesEquals( x, f );
  }
  
  
  
  /**
   * Pad a sequence with zeros (on the right)
   * 
   * If ypadded.length < y.length, this method performs a truncation of y.
   * 
   * @param y           float[] containing original sequence
   * @param ypadded     float[] containing padded sequence
   */
  public static void pad( float[] y, float[] ypadded ) {
    if ( y.length < ypadded.length ) {
      Arrays.fill( ypadded, 0.0f );
      System.arraycopy( y, 0, ypadded, 0, y.length );
    }
    else {
      System.arraycopy( y, 0, ypadded, 0, ypadded.length );
    }
  }
  
  
  
  /**
   * Pads this sequence to length n, by zero filling on right if n > length of this sequence, no-op otherwise.
   * 
   * @param n           int specifying desired length of padded sequence
   */
  public void pad( int n ) {
    if ( n > x.length ) {
      float[] tmp = new float[ n ];
      pad( x, tmp );
      x = tmp;
    }
  }
  
}
