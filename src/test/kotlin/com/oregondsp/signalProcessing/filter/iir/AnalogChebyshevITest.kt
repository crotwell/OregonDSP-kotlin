package com.oregondsp.signalProcessing.filter.iir

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec
import kotlin.js.Math

/**
 * Created by crotwell on 4/7/17.
 */
class AnalogChebyshevITest : StringSpec() {
    init {
        "strings.length should return size of string" {
            "hello".length shouldBe 5
        }

        "from main" {

            val A = AnalogChebyshevI(4, 0.50885)
            val B = A.lptolp(0.2 * Math.PI)

            val tmp = FloatArray(201)
            for (i in 0..200) {
                val C = B.evaluate(i * 0.02)
                tmp[i] = Complex.abs(C).toFloat()
            }

            val ps: PrintStream
            try {
                ps = PrintStream(FileOutputStream("C:\\DATA\\Test\\AnalogResponse.m"))
                ps.print("R = [ ")
                for (i in 0..199) {
                    ps.println(tmp[i].toString() + "; ...")
                }
                ps.println(tmp[200].toString() + "];")
                ps.close()
            } catch (e: FileNotFoundException) {
                e.printStackTrace()
            }

        }
        }
    }
}