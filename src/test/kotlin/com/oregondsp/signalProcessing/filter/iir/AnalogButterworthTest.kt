package com.oregondsp.signalProcessing.filter.iir

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec
import kotlin.js.Math

/**
 * Created by crotwell on 4/7/17.
 */
class AnalogButterworthTest  : StringSpec() {
    init {
        "strings.length should return size of string" {
            "hello".length shouldBe 5
        }

        "main from AnalogButterworth" {

            val B = AnalogButterworth(6)

            val A = B.lptobp(2.0 * Math.PI * 2.0, 2.0 * Math.PI * 3.0)

            val tmp = FloatArray(201)
            val gd = FloatArray(201)
            for (i in 0..200) {
                val omega = i.toDouble() * 2.0 * Math.PI / 20.0
                tmp[i] = Complex.abs(A.evaluate(omega)).toFloat()
                gd[i] = A.groupDelay(omega).toFloat()
            }

            var ps: PrintStream
            try {
                ps = PrintStream(FileOutputStream("C:\\DATA\\Test\\AnalogButterworthResponse.m"))
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
                ps = PrintStream(FileOutputStream("C:\\DATA\\Test\\AnalogButterworthGroupDelay.m"))
                ps.print("gd = [ ")
                for (i in 0..199) {
                    ps.println(gd[i].toString() + "; ...")
                }
                ps.println(gd[200].toString() + "];")
                ps.close()
            } catch (e: FileNotFoundException) {
                e.printStackTrace()
            }

        }
    }
}