package com.oregondsp.signalProcessing.filter.iir

/**
 * Created by crotwell on 4/7/17.
 */
class AnalogChebyshevIITest : StringSpec() {
    init {
        "strings.length should return size of string" {
            "hello".length shouldBe 5
        }

        "from main" {

            val A = AnalogChebyshevII(8, 0.1)

            val tmp = FloatArray(201)
            for (i in 0..200) {
                val C = A.evaluate(i * 0.02)
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