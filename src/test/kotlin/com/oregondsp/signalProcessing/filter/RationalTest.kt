package com.oregondsp.signalProcessing.filter

import com.oregondsp.signalProcessing.filter.iir.Complex
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

/**
 * Created by crotwell on 4/7/17.
 */
class RationalTest  : StringSpec() {
    init {
        "strings.length should return size of string" {
            "hello".length shouldBe 5
        }

        "rational from main, not sure what to test" {
            val a = DoubleArray(4)
            a[0] = 1.0
            a[1] = 2.0
            a[2] = 2.0
            a[3] = 1.0
            val b = DoubleArray(1)
            b[0] = 1.0
            val R = Rational(b, a)

            for (i in 0..99) {
                val omega = i / 25.0
                println(omega.toString() + "  " + R.evaluate(Complex(0.0, omega)).abs() + "   " + R.groupDelay(omega))
            }
        }
    }
}