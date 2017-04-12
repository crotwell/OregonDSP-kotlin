package com.oregondsp.signalProcessing.filter

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

/**
 * Created by crotwell on 4/7/17.
 */
class LagrangePolynomialTest  : StringSpec() {
    init {
        "strings.length should return size of string" {
            "hello".length shouldBe 5
        }

        "polynomian and lagrange equivalence" {

            val p = DoubleArray(3)
            p[0] = 6.0
            p[1] = -11.0
            p[2] = 6.0
            val P = Polynomial(p)

            val x = DoubleArray(3)
            x[0] = 1.0
            x[1] = 2.0
            x[2] = 3.0

            val f = DoubleArray(3)
            f[0] = 1.0
            f[1] = 8.0
            f[2] = 27.0
            val L = LagrangePolynomial(x, f)

            for (i in 0..20) {
                val z = 1.0 + i * 0.1
                P.evaluate(z).toString() shouldBe L.evaluate(z)
            }
        }
    }
}