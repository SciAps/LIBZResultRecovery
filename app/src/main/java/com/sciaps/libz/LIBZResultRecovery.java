package com.sciaps.libz;

import au.com.bytecode.opencsv.CSVWriter;
import com.devsmart.IOUtils;
import com.google.common.collect.ImmutableMap;
import com.sciaps.common.AtomicElement;
import com.sciaps.common.Global;
import com.sciaps.common.algorithms.GradeMatchRanker;
import com.sciaps.common.calculation.libs.EmpiricalCurveCalc;
import com.sciaps.common.calculation.libs.EmpiricalResult;
import com.sciaps.common.data.ChemResult;
import com.sciaps.common.data.SpectrometerCalibration;
import com.sciaps.common.data.flatbuff.*;
import com.sciaps.common.libs.LIBAnalysisResult;
import com.sciaps.common.spectrum.LIBZPixelSpectrum;
import org.apache.commons.cli.*;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.nio.ByteBuffer;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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


    private static final Pattern REGEX_RESULTFILENAME = Pattern.compile("alloyresult([0-9]+)\\.json");

    private static int getResultNum(File f) {
        final String filename = f.getName();
        Matcher m = REGEX_RESULTFILENAME.matcher(filename);
        if(m.find()) {
            return Integer.parseInt(m.group(1));
        } else {
            return -1;
        }
    }

    private static final int COLUMN_DATE = 0;
    private static final int COLUMN_TESTNAME = 1;
    private static final int COLUMN_TESTNUM = 2;
    private static final int COLUMN_MATCH1 = 3;
    private static final int COLUMN_MATCH2 = 4;
    private static final int COLUMN_MATCH3 = 5;
    private static final int COLUMN_MATCH4 = 6;
    private static final int COLUMN_MATCH5 = 7;
    private static final Map<AtomicElement, Integer> COLUMN_CHEM_MAP;
    private static final DateFormat DATE_FORMAT = SimpleDateFormat.getDateTimeInstance();

    static {
        ImmutableMap.Builder<AtomicElement, Integer> builder = new ImmutableMap.Builder<AtomicElement, Integer>();

        int column = 8;
        builder.put(AtomicElement.Silicon, column++);
        builder.put(AtomicElement.Iron, column++);
        builder.put(AtomicElement.Copper, column++);
        builder.put(AtomicElement.Manganese, column++);
        builder.put(AtomicElement.Magnesium, column++);
        builder.put(AtomicElement.Chromium, column++);
        builder.put(AtomicElement.Nickel, column++);
        builder.put(AtomicElement.Zinc, column++);
        builder.put(AtomicElement.Titanium, column++);
        builder.put(AtomicElement.getElementBySymbol("Ag"), column++);
        builder.put(AtomicElement.getElementBySymbol("Bi"), column++);
        builder.put(AtomicElement.getElementBySymbol("Li"), column++);
        builder.put(AtomicElement.getElementBySymbol("Pb"), column++);
        builder.put(AtomicElement.getElementBySymbol("Sn"), column++);
        builder.put(AtomicElement.getElementBySymbol("V"), column++);
        builder.put(AtomicElement.getElementBySymbol("Zr"), column++);
        builder.put(AtomicElement.getElementBySymbol("Al"), column++);

        COLUMN_CHEM_MAP = builder.build();

    }

    private static void doRecovery(final File rootDir, File output_csv) throws IOException {
        CSVWriter csvWriter = new CSVWriter(new FileWriter(output_csv));
        String[] line = new String[8 + COLUMN_CHEM_MAP.size()];
        line[COLUMN_DATE] = "Date";
        line[COLUMN_TESTNAME] = "Test Name";
        line[COLUMN_TESTNUM] = "Test #";
        line[COLUMN_MATCH1] = "Match #1";
        line[COLUMN_MATCH2] = "Match #2";
        line[COLUMN_MATCH3] = "Match #3";
        line[COLUMN_MATCH4] = "Match #4";
        line[COLUMN_MATCH5] = "Match #5";
        for(Map.Entry<AtomicElement, Integer> e : COLUMN_CHEM_MAP.entrySet()) {
            line[e.getValue()] = String.format("%s (%%)", e.getKey().symbol);
        }

        csvWriter.writeNext(line);



        int maxResultNum = 0;
        for(File f : rootDir.listFiles()) {
            maxResultNum = Math.max(maxResultNum, getResultNum(f));
        }

        for(int resultNum=0;resultNum<=maxResultNum;resultNum++) {
            File f = new File(rootDir, String.format("alloyresult%d.json", resultNum));
            if(f.exists()) {
                try {

                    for(int i=0;i<line.length;i++){
                        line[i] = "";
                    }

                    InputStream in = new FileInputStream(f);
                    LIBAnalysisResult r = loadResult(in);
                    in.close();

                    line[COLUMN_DATE] = DATE_FORMAT.format(r.mTime);
                    line[COLUMN_TESTNAME] = r.mTitle;
                    line[COLUMN_TESTNUM] = Integer.toString(resultNum);

                    if(r.mResult.gradeRanks != null && r.mResult.gradeRanks.size() >= 1) {
                        line[COLUMN_MATCH1] = r.mResult.gradeRanks.get(0).grade.getDisplayName();
                    }

                    if(r.mResult.gradeRanks != null && r.mResult.gradeRanks.size() >= 2) {
                        line[COLUMN_MATCH2] = r.mResult.gradeRanks.get(1).grade.getDisplayName();
                    }

                    if(r.mResult.gradeRanks != null && r.mResult.gradeRanks.size() >= 3) {
                        line[COLUMN_MATCH3] = r.mResult.gradeRanks.get(2).grade.getDisplayName();
                    }

                    if(r.mResult.gradeRanks != null && r.mResult.gradeRanks.size() >= 4) {
                        line[COLUMN_MATCH4] = r.mResult.gradeRanks.get(3).grade.getDisplayName();
                    }

                    if(r.mResult.gradeRanks != null && r.mResult.gradeRanks.size() >= 5) {
                        line[COLUMN_MATCH5] = r.mResult.gradeRanks.get(4).grade.getDisplayName();
                    }

                    if(r.mResult.chemResults != null) {
                        for(EmpiricalCurveCalc.EmpiricalCurveResult chemR : r.mResult.chemResults) {
                            Integer columnIndex = COLUMN_CHEM_MAP.get(chemR.element);
                            if(columnIndex != null) {
                                if(chemR.percent < 0 || chemR.type == ChemResult.TYPE_LESSLOD) {
                                    line[columnIndex] = "0";
                                } else {
                                    line[columnIndex] = String.format("%f", chemR.percent);
                                }

                            }
                        }
                    }

                    csvWriter.writeNext(line);


                } catch (Exception e) {
                    LOGGER.error("", e);
                }
            }
        }

        csvWriter.close();
    }

    public static void main(String[] args) {

        Options options = new Options();
        options.addOption(Option.builder("d")
                .hasArg()
                .argName("libsresults")
                .required()
                .build());

        options.addOption(Option.builder("o")
                .hasArg()
                .argName("output.csv")
                .build());

        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);

            File libsresults_root = new File(line.getOptionValue("d"));
            File output_csv = new File(line.getOptionValue("o", "output.csv"));
            doRecovery(libsresults_root, output_csv);


        } catch (ParseException e) {
            System.err.println(e.getMessage());

            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( "LIBZResultRecovery", options );

            System.exit(-1);
        } catch (Exception e) {
            System.err.println(e.getMessage());
            System.exit(-1);
        }

    }


}
