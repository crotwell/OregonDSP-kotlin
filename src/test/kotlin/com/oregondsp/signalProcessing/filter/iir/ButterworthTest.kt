package com.oregondsp.signalProcessing.filter.iir

import kotlin.js.Math
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

/**
 * Created by crotwell on 4/7/17.
 */
class ButterworthTest : StringSpec() {

    init {
        "from main" {

            val B = Butterworth(3, PassbandType.BANDPASS, 2.0, 5.0, 0.025)
            B.print(System.out)
            val tmp = FloatArray(201)
            for (i in 0..200) {
                val C = B.evaluate(Math.PI / 200.0 * i)
                tmp[i] = Complex.abs(C).toFloat()
            }

            val x = FloatArray(1001)
            x[200] = 1.0f
            val y = FloatArray(1001)
            B.filter(x, y)

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