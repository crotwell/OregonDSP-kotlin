package com.oregondsp.signalProcessing.filter.iir

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec
import kotlin.js.Math

/**
 * Created by crotwell on 4/7/17.
 */
class ChebyshevITest : StringSpec() {
    init {
        "strings.length should return size of string" {
            "hello".length shouldBe 5
        }

        "from main" {

            val epsilon = 0.50885
            val F = ChebyshevI(4, epsilon, PassbandType.LOWPASS, 2.0, 0.0, 0.05)
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