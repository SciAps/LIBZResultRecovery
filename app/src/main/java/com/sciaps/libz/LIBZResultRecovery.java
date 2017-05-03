package com.sciaps.libz;

import com.devsmart.IOUtils;
import com.sciaps.common.AtomicElement;
import com.sciaps.common.Global;
import com.sciaps.common.algorithms.GradeMatchRanker;
import com.sciaps.common.calculation.libs.EmpiricalCurveCalc;
import com.sciaps.common.calculation.libs.EmpiricalResult;
import com.sciaps.common.data.SpectrometerCalibration;
import com.sciaps.common.data.flatbuff.*;
import com.sciaps.common.libs.LIBAnalysisResult;
import com.sciaps.common.spectrum.LIBZPixelSpectrum;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Date;

public class LIBZResultRecovery {

    private static final Logger LOGGER = LoggerFactory.getLogger(LIBZResultRecovery.class);

    private static LIBZPixelSpectrum loadSpectrum(PixSpectrum pixSpectrum) {
        double[] knots = new double[pixSpectrum.knotsLength()];
        for (int i = 0; i < pixSpectrum.knotsLength(); i++) {
            knots[i] = pixSpectrum.knots(i);
        }

        double[][] pixels = new double[pixSpectrum.dataLength()][];
        SpectrometerCalibration[] wlcal = new SpectrometerCalibration[pixSpectrum.dataLength()];

        SpectrometerPixData data = new SpectrometerPixData();
        for (int i = 0; i < pixSpectrum.dataLength(); i++) {
            pixSpectrum.data(data, i);

            double[] pixNmCoeff = new double[data.pixToNmLength()];
            for (int q = 0; q < data.pixToNmLength(); q++) {
                pixNmCoeff[q] = data.pixToNm(q);
            }

            wlcal[i] = new SpectrometerCalibration();
            wlcal[i].pixToNm = new PolynomialFunction(pixNmCoeff);

            pixels[i] = Global.ALLOCATOR.alloc(data.pixLength());
            for (int q = 0; q < data.pixLength(); q++) {
                pixels[i][q] = data.pix(q);
            }
        }

        LIBZPixelSpectrum retval = new LIBZPixelSpectrum();
        retval.setData(pixels, wlcal, knots);
        retval.shotNum = pixSpectrum.shotNum();

        for (int i = 0; i < pixels.length; i++) {
            Global.ALLOCATOR.free(pixels[i]);
        }

        return retval;
    }

    private static com.sciaps.common.data.Grade loadGrade(Grade grade) {
        com.sciaps.common.data.Grade retval = new com.sciaps.common.data.Grade();

        retval.name = grade.name();

        ChemRange chemRange = new ChemRange();
        for (int i = 0; i < grade.specLength(); i++) {
            if (grade.spec(chemRange, i) != null) {
                com.sciaps.common.data.Grade.Spec spec = new com.sciaps.common.data.Grade.Spec();
                spec.min = chemRange.min();
                spec.max = chemRange.max();
                retval.spec.put(
                        AtomicElement.getElementByAtomicNum(chemRange.element()),
                        spec);
            }
        }

        return retval;
    }


    public static LIBAnalysisResult loadResult(InputStream in) throws IOException {
        //71k should be about a typical result size
        ByteArrayOutputStream bout = new ByteArrayOutputStream(72704);
        IOUtils.pump(in, bout);

        Test test = Test.getRootAsTest(ByteBuffer.wrap(bout.toByteArray()));

        LIBAnalysisResult retval = new LIBAnalysisResult();
        retval.mTime = new Date(test.unixTime());
        retval.mTitle = test.displayName();
        LOGGER.info("loadResult | retval.mTitle: " + retval.mTitle);

        AlloyResult alloyResult = new AlloyResult();
        test.result(alloyResult);
        retval.mResult = new EmpiricalResult();
        retval.mResult.base = alloyResult.base();

        retval.mResult.spectrum = loadSpectrum(test.shots(0));

        ChemValue chemValue = new ChemValue();
        retval.mResult.chemResults = new ArrayList<EmpiricalCurveCalc.EmpiricalCurveResult>(alloyResult.chemResultsLength());
        for (int i = 0; i < alloyResult.chemResultsLength(); i++) {
            if (alloyResult.chemResults(chemValue, i) != null) {
                EmpiricalCurveCalc.EmpiricalCurveResult chemResult = new EmpiricalCurveCalc.EmpiricalCurveResult();
                chemResult.element = AtomicElement.getElementByAtomicNum(chemValue.element());
                chemResult.percent = chemValue.value();
                chemResult.error = chemValue.stdErr();

                retval.mResult.chemResults.add(chemResult);
            }
        }

        try {
            GradeMatch gradeMatch = new GradeMatch();
            retval.mResult.gradeRanks = new ArrayList<GradeMatchRanker.GradeRank>(alloyResult.gradeMatchesLength());
            final int numGradeMatches = alloyResult.gradeMatchesLength();
            for (int i = 0; i < numGradeMatches; i++) {
                if (alloyResult.gradeMatches(gradeMatch, i) != null) {
                    GradeMatchRanker.GradeRank gradeResult = new GradeMatchRanker.GradeRank();
                    gradeResult.grade = loadGrade(gradeMatch.grade());
                    gradeResult.matchNumber = gradeMatch.match();

                    retval.mResult.gradeRanks.add(gradeResult);
                }
            }
        } catch (Exception e) {
            LOGGER.error(e.toString(), e);
        }

        return retval;
    }


    public static void main(String[] args) {

    }
}
