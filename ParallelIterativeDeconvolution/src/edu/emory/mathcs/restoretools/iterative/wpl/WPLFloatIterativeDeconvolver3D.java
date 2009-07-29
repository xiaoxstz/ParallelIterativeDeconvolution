/*  License:
 Copyright (c) 2005, OptiNav, Inc.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

 Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 Neither the name of OptiNav, Inc. nor the names of its contributors
 may be used to endorse or promote products derived from this software
 without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package edu.emory.mathcs.restoretools.iterative.wpl;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import cern.colt.matrix.tfloat.FloatMatrix3D;
import cern.colt.matrix.tfloat.impl.DenseFloatMatrix3D;
import cern.jet.math.tfloat.FloatFunctions;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.FloatCommon3D;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PaddingType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.ResizingType;
import edu.emory.mathcs.utils.ConcurrencyUtils;

/**
 * Wiener Filter Preconditioned Landweber 3D. This is a nonnegatively constrained method.
 * 
 * @author Bob Dougherty
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class WPLFloatIterativeDeconvolver3D {
    /**
     * Blurred image
     */
    protected FloatMatrix3D B;
    /**
     * Point Spread Function
     */
    protected FloatMatrix3D PSF;
    /**
     * Color model
     */
    protected java.awt.image.ColorModel cmY;

    /**
     * Number of slices in the blurred image.
     */
    protected int bSlices;

    /**
     * Number of columns in the blurred image.
     */
    protected int bColumns;

    /**
     * Number of rows in the blurred image.
     */
    protected int bRows;

    /**
     * Maximal number of iterations.
     */
    protected int maxIters;

    /**
     * If true, then the thresholding is performed.
     */
    protected boolean useThreshold;

    /**
     * The smallest nonnegative value assigned to the restored image.
     */
    protected float threshold;

    /**
     * If true, then the restored image is displayed after each iteration.
     */
    protected boolean showIteration;

    /**
     * If true, then the convergence information is displayed after each
     * iteration.
     */
    protected boolean logConvergence;

    /**
     * Type of restored image.
     */
    protected OutputType output;

    /**
     * Number of slices in the expanded blurred image.
     */
    private int slices;

    /**
     * Number of columns in the expanded blurred image.
     */
    private int columns;

    /**
     * Number of rows in the expanded blurred image.
     */
    private int rows;

    /**
     * Min value of the blurred image.
     */
    private float minB = 0;

    /**
     * Min value of the PSF.
     */
    private float minPSF = 0;

    /**
     * Sum of all elements in the PSF matrix.
     */
    private float sum;

    /**
     * Scaling factor.
     */
    private float scalePSF = 1;

    /**
     * Gaussian weights
     */
    private float[][] gweights;

    /**
     * Regularization parameter for the Wiener Filter.
     */
    protected float gamma;

    /**
     * Number of pixels in x and y directions for low-pass filter.
     */
    protected float filterXY;

    /**
     * Number of pixels in z direction for low-pass filter.
     */
    protected float filterZ;
    /**
     * If true, then PSF is normalized.
     */
    protected boolean normalize;

    /**
     * If true, then the anti-ringing step is performed.
     */
    protected boolean antiRing;

    /**
     * This parameter is used to stop the iteration if the image is not
     * changing.
     */
    protected float changeThreshPercent;

    /**
     * If true, then all the data is in decibels.
     */
    protected boolean dB;

    /**
     * If true, then the iterations are stopped when the changes appear to be
     * increasing.
     */
    protected boolean detectDivergence;

    /**
     * Crates a new instance of WPLFloatIterativeDeconvolver3D
     * 
     * @param imB
     *            blurred image
     * @param imPSF
     *            Point Spread Function
     * @param boundary
     *            type of boundary conditions
     * @param resizing
     *            type of resizing
     * @param output
     *            type of the output image
     * @param maxIters
     *            maximal number of iterations
     * @param showIteration
     *            if true, then the restored image is displayed after each
     *            iteration
     * @param options
     *            WPL options
     */
    public WPLFloatIterativeDeconvolver3D(ImagePlus imB, ImagePlus imPSF, BoundaryType boundary, ResizingType resizing, OutputType output, int maxIters, boolean showIteration, WPLOptions options) {
        IJ.showStatus("WPL initialization...");
        ImageProcessor ipB = imB.getProcessor();
        if (output == OutputType.SAME_AS_SOURCE) {
            if (ipB instanceof ByteProcessor) {
                this.output = OutputType.BYTE;
            } else if (ipB instanceof ShortProcessor) {
                this.output = OutputType.SHORT;
            } else if (ipB instanceof FloatProcessor) {
                this.output = OutputType.FLOAT;
            } else {
                throw new IllegalArgumentException("Unsupported image type.");
            }
        } else {
            this.output = output;
        }
        ImageStack isB = imB.getStack();
        cmY = ipB.getColorModel();
        bSlices = imB.getStackSize();
        bColumns = ipB.getWidth();
        bRows = ipB.getHeight();
        B = new DenseFloatMatrix3D(bSlices, bRows, bColumns);
        FloatCommon3D.assignPixelsToMatrix(isB, B);

        ImageStack isPSF = imPSF.getStack();
        ImageProcessor ipPSF = imPSF.getProcessor();
        int psfSlices = imPSF.getStackSize();
        int psfColumns = ipPSF.getWidth();
        int psfRows = ipPSF.getHeight();
        PSF = new DenseFloatMatrix3D(psfSlices, psfRows, psfColumns);
        FloatCommon3D.assignPixelsToMatrix(isPSF, PSF);
        this.maxIters = maxIters;
        this.showIteration = showIteration;
        this.gamma = (float) options.getGamma();
        this.filterXY = (float) options.getFilterXY();
        this.filterZ = (float) options.getFilterZ();
        this.normalize = options.isNormalize();
        this.antiRing = options.isAntiRing();
        this.changeThreshPercent = (float) options.getChangeThreshPercent();
        this.dB = options.isDB();
        this.detectDivergence = options.isDetectDivergence();
        this.logConvergence = options.isLogConvergence();
        if (this.dB) {
            minB = unDB(B);
            minPSF = unDB(PSF);
        }
        sum = PSF.zSum();
        if ((sum != 0) && this.normalize)
            scalePSF /= sum;

        slices = expandedSize(psfSlices, bSlices, resizing);
        columns = expandedSize(psfColumns, bColumns, resizing);
        rows = expandedSize(psfRows, bRows, resizing);
        if ((psfSlices > slices) || (psfColumns > columns) || (psfRows > rows)) {
            throw new IllegalArgumentException("PSF cannot be largest that the image.");
        }
        gweights = gaussianWeights(slices, rows, columns, this.filterXY, this.filterXY, this.filterZ);
        switch (boundary) {
        case PERIODIC:
            B = FloatCommon3D.padPeriodic(B, slices, rows, columns);
            break;
        case REFLEXIVE:
            B = FloatCommon3D.padReflexive(B, slices, rows, columns);
            break;
        case ZERO:
            B = FloatCommon3D.padZero(B, slices, rows, columns);
            break;
        }
        float[] maxLoc = PSF.getMaxLocation();
        int[] padSize = new int[3];
        padSize[0] = slices - psfSlices;
        padSize[1] = rows - psfRows;
        padSize[2] = columns - psfColumns;
        PSF = FloatCommon3D.padZero(PSF, padSize, PaddingType.POST);
        PSF = FloatCommon3D.circShift(PSF, new int[] { (int) maxLoc[1], (int) maxLoc[2], (int) maxLoc[3] });
    }

    /**
     * Performs deconvolution and returns deconvolved image.
     * 
     * @return deconvolved image
     */
    public ImagePlus deconvolve() {
        ((DenseFloatMatrix3D) PSF).dht3();
        FloatMatrix3D X;
        FloatMatrix3D AX = B.like();
        if (antiRing) {
            IJ.showStatus("WPL: performing anti-ringing step.");
            X = B.copy();
            ((DenseFloatMatrix3D) X).dht3();
            convolveFD(slices, rows, columns, PSF, X, AX);
            ((DenseFloatMatrix3D) AX).idht3(true);
            copyDataAverage(bSlices, bRows, bColumns, slices, rows, columns, sum, B, AX, B);
        }
        if (gamma > 0.0001) {
            IJ.showStatus("WPL: Wiener filter");
            float magMax = findMagMax(PSF);
            ((DenseFloatMatrix3D) B).dht3();
            X = PSF.copy();
            deconvolveFD(gamma, magMax, slices, rows, columns, X, X, PSF);
            AX = B.copy();
            deconvolveFD(gamma, magMax, slices, rows, columns, AX, X, B);
            ((DenseFloatMatrix3D) B).idht3(true);
        }

        int sOff = (slices - bSlices + 1) / 2;
        int rOff = (rows - bRows + 1) / 2;
        int cOff = (columns - bColumns + 1) / 2;

        ((DenseFloatMatrix3D) PSF).idht3(true);
        float aSum = PSF.aggregate(FloatFunctions.plus, FloatFunctions.abs);
        if (scalePSF != 1) {
            B.assign(FloatFunctions.div(scalePSF));
        }
        ((DenseFloatMatrix3D) PSF).dht3();
        X = B.copy();
        ImagePlus imX = null;
        ImageStack is = new ImageStack(bColumns, bRows);
        if (showIteration) {
            FloatCommon3D.assignPixelsToStackPadded(is, X, bSlices, bRows, bColumns, sOff, rOff, cOff, cmY);
            imX = new ImagePlus("(deblurred)", is);
            imX.show();
        }
        float oldPercentChange = Float.MAX_VALUE;
        for (int iter = 0; iter < maxIters; iter++) {
            IJ.showStatus("WPL iteration: " + (iter + 1) + "/" + maxIters);
            ((DenseFloatMatrix3D) X).dht3();
            gaussianFilter(X, gweights);
            convolveFD(slices, rows, columns, PSF, X, AX);
            ((DenseFloatMatrix3D) AX).idht3(true);
            ((DenseFloatMatrix3D) X).idht3(true);
            float meanDelta = meanDelta(B, AX, X, aSum);
            if (showIteration) {
                if (threshold == -1.0) {
                    FloatCommon3D.updatePixelsInStackPadded(is, X, bSlices, bRows, bColumns, sOff, rOff, cOff, cmY);
                } else {
                    FloatCommon3D.updatePixelsInStackPadded(is, X, bSlices, bRows, bColumns, sOff, rOff, cOff, cmY, threshold);
                }
                ImageProcessor ip1 = imX.getProcessor();
                ip1.setMinAndMax(0, 0);
                ip1.setColorModel(cmY);
                imX.updateAndDraw();
            }
            float sumPixels = energySum(X, bSlices, bRows, bColumns, sOff, rOff, cOff);
            float percentChange = 100 * meanDelta / sumPixels;
            if (logConvergence)
                IJ.write(Float.toString(percentChange));
            if ((oldPercentChange - percentChange) < changeThreshPercent) {
                if (logConvergence)
                    IJ.write("Automatically terminated after " + (iter + 1) + " iterations.");
                break;
            }
            if ((oldPercentChange < percentChange) && detectDivergence) {
                if (logConvergence)
                    IJ.write("Automatically terminated due to divergence " + (iter + 1) + " iterations.");
                break;
            }
            oldPercentChange = percentChange;
        }
        ((DenseFloatMatrix3D) X).dht3();
        gaussianFilterWithScaling(X, gweights, aSum);
        ((DenseFloatMatrix3D) X).idht3(true);
        if (dB) {
            toDB(PSF, minPSF);
            toDB(B, minB);
            toDB(X, -90);
        }
        if (!showIteration) {
            if (threshold == -1.0) {
                FloatCommon3D.assignPixelsToStackPadded(is, X, bSlices, bRows, bColumns, sOff, rOff, cOff, cmY);
            } else {
                FloatCommon3D.assignPixelsToStackPadded(is, X, bSlices, bRows, bColumns, sOff, rOff, cOff, cmY, threshold);
            }
            imX = new ImagePlus("(deblurred)", is);
        }
        FloatCommon3D.convertImage(imX, output);
        return imX;
    }

    private static void convolveFD(final int slices, final int rows, final int columns, FloatMatrix3D H1, FloatMatrix3D H2, FloatMatrix3D Result) {
        final float[] h1 = (float[]) H1.elements();
        final float[] h2 = (float[]) H2.elements();
        final float[] result = (float[]) Result.elements();
        final int sliceStride = columns * rows;
        final int rowStride = columns;

        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (slices * columns * rows >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            int k = slices / np;
            for (int j = 0; j < np; j++) {
                final int firstSlice = j * k;
                final int lastSlice = (j == np - 1) ? slices : firstSlice + k;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {
                    public void run() {
                        int sC, cC, rC, idx1, idx2;
                        float h2e, h2o;
                        for (int s = firstSlice; s < lastSlice; s++) {
                            sC = (slices - s) % slices;
                            for (int r = 0; r < rows; r++) {
                                rC = (rows - r) % rows;
                                for (int c = 0; c < columns; c += 2) {
                                    cC = (columns - c) % columns;
                                    idx1 = c + rowStride * r + sliceStride * s;
                                    idx2 = cC + rowStride * rC + sliceStride * sC;
                                    h2e = (h2[idx1] + h2[idx2]) / 2;
                                    h2o = (h2[idx1] - h2[idx2]) / 2;
                                    result[idx1] = (float) (h1[idx1] * h2e + h1[idx2] * h2o);
                                    cC = (columns - c - 1) % columns;
                                    idx1 = c + 1 + rowStride * r + sliceStride * s;
                                    idx2 = cC + rowStride * rC + sliceStride * sC;
                                    h2e = (h2[idx1] + h2[idx2]) / 2;
                                    h2o = (h2[idx1] - h2[idx2]) / 2;
                                    result[idx1] = (float) (h1[idx1] * h2e + h1[idx2] * h2o);
                                }
                            }
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int sC, cC, rC, idx1, idx2;
            float h2e, h2o;
            for (int s = 0; s < slices; s++) {
                sC = (slices - s) % slices;
                for (int r = 0; r < rows; r++) {
                    rC = (rows - r) % rows;
                    for (int c = 0; c < columns; c += 2) {
                        cC = (columns - c) % columns;
                        idx1 = c + rowStride * r + sliceStride * s;
                        idx2 = cC + rowStride * rC + sliceStride * sC;
                        h2e = (h2[idx1] + h2[idx2]) / 2;
                        h2o = (h2[idx1] - h2[idx2]) / 2;
                        result[idx1] = (float) (h1[idx1] * h2e + h1[idx2] * h2o);
                        cC = (columns - c - 1) % columns;
                        idx1 = c + 1 + rowStride * r + sliceStride * s;
                        idx2 = cC + rowStride * rC + sliceStride * sC;
                        h2e = (h2[idx1] + h2[idx2]) / 2;
                        h2o = (h2[idx1] - h2[idx2]) / 2;
                        result[idx1] = (float) (h1[idx1] * h2e + h1[idx2] * h2o);
                    }
                }
            }
        }
    }

    private static void copyDataAverage(final int slices, final int rows, final int columns, final int slicesE, final int rowsE, final int columnsE, final float sum, FloatMatrix3D DataIn, FloatMatrix3D DataOut, FloatMatrix3D Result) {
        final float[] dataIn = (float[]) DataIn.elements();
        final float[] dataOut = (float[]) DataOut.elements();
        final float[] result = (float[]) Result.elements();
        final int sOff = (slicesE - slices + 1) / 2;
        final int rOff = (rowsE - rows + 1) / 2;
        final int cOff = (columnsE - columns + 1) / 2;
        final int sliceStride = rowsE * columnsE;
        final int rowStride = columnsE;

        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (slicesE * columnsE * rowsE >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            int k = slicesE / np;
            for (int j = 0; j < np; j++) {
                final int firstSlice = -sOff + j * k;
                final int lastSlice = (j == np - 1) ? slicesE - sOff : firstSlice + k;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {

                    public void run() {
                        int sOut, cOut, rOut, idx;
                        float alphaS, alphaC, alphaR;
                        float a;
                        for (int s = firstSlice; s < lastSlice; s++) {
                            sOut = s + sOff;
                            if (s < 0) {
                                alphaS = -s / ((float) sOff);
                            } else if (s > (slices - 1)) {
                                alphaS = (s - slices) / ((float) sOff);
                            } else {
                                alphaS = 0;
                            }
                            for (int r = -rOff; r < rowsE - rOff; r++) {
                                rOut = r + rOff;
                                if (r < 0) {
                                    alphaR = -r / ((float) rOff);
                                } else if (r > (rows - 1)) {
                                    alphaR = (r - rows) / ((float) rOff);
                                } else {
                                    alphaR = 0;
                                }
                                for (int c = -cOff; c < columnsE - cOff; c++) {
                                    cOut = c + cOff;
                                    if (c < 0) {
                                        alphaC = -c / ((float) cOff);
                                    } else if (c > (columns - 1)) {
                                        alphaC = (c - columns) / ((float) cOff);
                                    } else {
                                        alphaC = 0;
                                    }
                                    a = alphaS;
                                    if (alphaR > a)
                                        a = alphaR;
                                    if (alphaC > a)
                                        a = alphaC;
                                    idx = cOut + rowStride * rOut + sliceStride * sOut;
                                    result[idx] = (1 - a) * dataIn[idx] + a * dataOut[idx] / sum;
                                }
                            }
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int sOut, cOut, rOut, idx;
            float alphaS, alphaC, alphaR;
            float a;
            for (int s = -sOff; s < slicesE - sOff; s++) {
                sOut = s + sOff;
                if (s < 0) {
                    alphaS = -s / ((float) sOff);
                } else if (s > (slices - 1)) {
                    alphaS = (s - slices) / ((float) sOff);
                } else {
                    alphaS = 0;
                }
                for (int r = -rOff; r < rowsE - rOff; r++) {
                    rOut = r + rOff;
                    if (r < 0) {
                        alphaR = -r / ((float) rOff);
                    } else if (r > (rows - 1)) {
                        alphaR = (r - rows) / ((float) rOff);
                    } else {
                        alphaR = 0;
                    }
                    for (int c = -cOff; c < columnsE - cOff; c++) {
                        cOut = c + cOff;
                        if (c < 0) {
                            alphaC = -c / ((float) cOff);
                        } else if (c > (columns - 1)) {
                            alphaC = (c - columns) / ((float) cOff);
                        } else {
                            alphaC = 0;
                        }
                        a = alphaS;
                        if (alphaR > a)
                            a = alphaR;
                        if (alphaC > a)
                            a = alphaC;
                        idx = cOut + rowStride * rOut + sliceStride * sOut;
                        result[idx] = (1 - a) * dataIn[idx] + a * dataOut[idx] / sum;
                    }
                }
            }
        }
    }

    private static void deconvolveFD(final float gamma, final float magMax, final int slices, final int rows, final int columns, FloatMatrix3D H1, FloatMatrix3D H2, FloatMatrix3D Result) {
        final float gammaScaled = gamma * magMax;
        final float[] h1 = (float[]) H1.elements();
        final float[] h2 = (float[]) H2.elements();
        final float[] result = (float[]) Result.elements();
        final int sliceStride = columns * rows;
        final int rowStride = columns;
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (slices * columns * rows >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            int k = slices / np;
            for (int j = 0; j < np; j++) {
                final int firstSlice = j * k;
                final int lastSlice = (j == np - 1) ? slices : firstSlice + k;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {
                    public void run() {
                        int sC, cC, rC, idx1, idx2;
                        float mag, h2e, h2o;
                        for (int s = firstSlice; s < lastSlice; s++) {
                            sC = (slices - s) % slices;
                            for (int r = 0; r < rows; r++) {
                                rC = (rows - r) % rows;
                                for (int c = 0; c < columns; c++) {
                                    cC = (columns - c) % columns;
                                    idx1 = c + rowStride * r + sliceStride * s;
                                    idx2 = cC + rowStride * rC + sliceStride * sC;
                                    h2e = (h2[idx1] + h2[idx2]) / 2;
                                    h2o = (h2[idx1] - h2[idx2]) / 2;
                                    mag = h2[idx1] * h2[idx1] + h2[idx2] * h2[idx2];
                                    float tmp = h1[idx1] * h2e - h1[idx2] * h2o;
                                    result[idx1] = (tmp / (mag + gammaScaled));
                                }
                            }
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int sC, cC, rC, idx1, idx2;
            float mag, h2e, h2o;
            for (int s = 0; s < slices; s++) {
                sC = (slices - s) % slices;
                for (int r = 0; r < rows; r++) {
                    rC = (rows - r) % rows;
                    for (int c = 0; c < columns; c++) {
                        cC = (columns - c) % columns;
                        idx1 = c + rowStride * r + sliceStride * s;
                        idx2 = cC + rowStride * rC + sliceStride * sC;
                        h2e = (h2[idx1] + h2[idx2]) / 2;
                        h2o = (h2[idx1] - h2[idx2]) / 2;
                        mag = h2[idx1] * h2[idx1] + h2[idx2] * h2[idx2];
                        float tmp = h1[idx1] * h2e - h1[idx2] * h2o;
                        result[idx1] = (tmp / (mag + gammaScaled));
                    }
                }
            }
        }
    }

    private static float energySum(FloatMatrix3D X, final int slices, final int rows, final int columns, final int sOff, final int rOff, final int cOff) {
        float sumPixels = 0;
        final int rowStride = X.rowStride();
        final int sliceStride = X.sliceStride();
        final float[] elemsX = (float[]) X.elements();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (slices * rows * columns >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            Float[] results = new Float[np];
            int k = slices / np;
            for (int j = 0; j < np; j++) {
                final int firstSlice = j * k;
                final int lastSlice = (j == np - 1) ? slices : firstSlice + k;
                futures[j] = ConcurrencyUtils.submit(new Callable<Float>() {
                    public Float call() throws Exception {
                        float sumPixels = 0;
                        for (int s = firstSlice; s < lastSlice; s++) {
                            for (int r = 0; r < rows; r++) {
                                for (int c = 0; c < columns; c++) {
                                    sumPixels += elemsX[c + cOff + rowStride * (r + rOff) + sliceStride * (s + sOff)];
                                }
                            }
                        }
                        return sumPixels;
                    }
                });
            }
            try {
                for (int j = 0; j < np; j++) {
                    results[j] = (Float) futures[j].get();
                }
            } catch (ExecutionException ex) {
                ex.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            for (int j = 0; j < np; j++) {
                sumPixels += results[j];
            }
        } else {
            for (int s = 0; s < slices; s++) {
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < columns; c++) {
                        sumPixels += elemsX[c + cOff + rowStride * (r + rOff) + sliceStride * (s + sOff)];
                    }
                }
            }
        }
        return sumPixels;
    }

    private static int expandedSize(int psfSize, int bSize, ResizingType resizing) {
        int result = 0;
        int minimal = psfSize + bSize;
        switch (resizing) {
        case AUTO:
            int nextPowTwo;
            if (!ConcurrencyUtils.isPowerOf2(minimal)) {
                nextPowTwo = ConcurrencyUtils.nextPow2(minimal);
            } else {
                nextPowTwo = minimal;
            }
            if (nextPowTwo >= 1.5 * minimal) {
                //use minimal padding
                result = minimal;
            } else {
                result = nextPowTwo;
            }
            break;
        case MINIMAL:
            result = minimal;
            break;
        case NEXT_POWER_OF_TWO:
            result = minimal;
            if (!ConcurrencyUtils.isPowerOf2(result)) {
                result = ConcurrencyUtils.nextPow2(result);
            }
            break;
        }
        if (result < 4) {
            result = 4;
        }
        return result;
    }

    private static float findMagMax(FloatMatrix3D H2) {
        final float[] h2 = (float[]) H2.elements();
        float magMax = 0;
        final int slices = H2.slices();
        final int rows = H2.rows();
        final int columns = H2.columns();
        final int sliceStride = rows * columns;
        final int rowStride = columns;
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (slices * rows * columns >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            Float[] results = new Float[np];
            int k = slices / np;
            for (int j = 0; j < np; j++) {
                final int firstSlice = j * k;
                final int lastSlice = (j == np - 1) ? slices : firstSlice + k;
                futures[j] = ConcurrencyUtils.submit(new Callable<Float>() {
                    public Float call() throws Exception {
                        int sC, cC, rC, idx1, idx2;
                        float magMax = 0;
                        float mag;
                        for (int s = firstSlice; s < lastSlice; s++) {
                            sC = (slices - s) % slices;
                            for (int r = 0; r < rows; r++) {
                                rC = (rows - r) % rows;
                                for (int c = 0; c < columns; c++) {
                                    cC = (columns - c) % columns;
                                    idx1 = c + rowStride * r + sliceStride * s;
                                    idx2 = cC + rowStride * rC + sliceStride * sC;
                                    mag = h2[idx1] * h2[idx1] + h2[idx2] * h2[idx2];
                                    if (mag > magMax)
                                        magMax = mag;
                                }
                            }
                        }
                        return magMax;
                    }
                });
            }
            try {
                for (int j = 0; j < np; j++) {
                    results[j] = (Float) futures[j].get();
                }
            } catch (ExecutionException ex) {
                ex.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            magMax = results[0];
            for (int j = 1; j < np; j++) {
                if (results[j] > magMax)
                    magMax = results[j];
            }
        } else {
            int sC, cC, rC, idx1, idx2;
            float mag;
            for (int s = 0; s < slices; s++) {
                sC = (slices - s) % slices;
                for (int r = 0; r < rows; r++) {
                    rC = (rows - r) % rows;
                    for (int c = 0; c < columns; c++) {
                        cC = (columns - c) % columns;
                        idx1 = c + rowStride * r + sliceStride * s;
                        idx2 = cC + rowStride * rC + sliceStride * sC;
                        mag = h2[idx1] * h2[idx1] + h2[idx2] * h2[idx2];
                        if (mag > magMax)
                            magMax = mag;
                    }
                }
            }
        }
        return magMax;
    }

    private static void gaussianFilter(FloatMatrix3D X, final float[][] weights) {
        final float[] elems = (float[]) X.elements();
        final int sliceStride = X.sliceStride();
        final int rowStride = X.rowStride();
        final int columnStride = X.columnStride();
        final int slices = X.slices();
        final int rows = X.rows();
        final int columns = X.columns();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (slices * columns * rows >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            int k = slices / np;
            for (int j = 0; j < np; j++) {
                final int firstSlice = j * k;
                final int lastSlice = (j == np - 1) ? slices : firstSlice + k;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {
                    public void run() {
                        int idx;
                        for (int s = firstSlice; s < lastSlice; s++) {
                            for (int r = 0; r < rows; r++) {
                                idx = s * sliceStride + r * rowStride;
                                for (int c = 0; c < columns; c++) {
                                    elems[idx] *= weights[2][s] * weights[1][r] * weights[0][c];
                                    idx += columnStride;
                                }
                            }
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx;
            for (int s = 0; s < slices; s++) {
                for (int r = 0; r < rows; r++) {
                    idx = s * sliceStride + r * rowStride;
                    for (int c = 0; c < columns; c++) {
                        elems[idx] *= weights[2][s] * weights[1][r] * weights[0][c];
                        idx += columnStride;
                    }
                }
            }
        }
    }

    private static void gaussianFilterWithScaling(FloatMatrix3D X, final float[][] weights, final float scale) {
        final float[] elems = (float[]) X.elements();
        final int sliceStride = X.sliceStride();
        final int rowStride = X.rowStride();
        final int columnStride = X.columnStride();
        final int slices = X.slices();
        final int rows = X.rows();
        final int columns = X.columns();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (slices * columns * rows >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            int k = slices / np;
            for (int j = 0; j < np; j++) {
                final int firstSlice = j * k;
                final int lastSlice = (j == np - 1) ? slices : firstSlice + k;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {
                    public void run() {
                        int idx;
                        for (int s = firstSlice; s < lastSlice; s++) {
                            for (int r = 0; r < rows; r++) {
                                idx = s * sliceStride + r * rowStride;
                                for (int c = 0; c < columns; c++) {
                                    elems[idx] *= weights[2][s] * weights[1][r] * weights[0][c] / scale;
                                    idx += columnStride;
                                }
                            }
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx;
            for (int s = 0; s < slices; s++) {
                for (int r = 0; r < rows; r++) {
                    idx = s * sliceStride + r * rowStride;
                    for (int c = 0; c < columns; c++) {
                        elems[idx] *= weights[2][s] * weights[1][r] * weights[0][c] / scale;
                        idx += columnStride;
                    }
                }
            }
        }
    }

    private static float[][] gaussianWeights(final int slices, final int rows, final int columns, final float filterX, final float filterY, final float filterZ) {
        final float[][] weights = new float[3][];
        weights[0] = new float[columns];
        weights[1] = new float[rows];
        weights[2] = new float[slices];
        final float cc = (float) (columns / (filterX + 0.000001));
        final float rc = (float) (rows / (filterY + 0.000001));
        final float sc = (float) (slices / (filterZ + 0.000001));
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (Math.max(slices, Math.max(columns, rows)) >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            Future<?>[] futures = new Future[np];
            int kcol = columns / np;
            int krow = rows / np;
            int ksls = slices / np;
            for (int j = 0; j < np; j++) {
                final int firstCol = j * kcol;
                final int lastCol = (j == np - 1) ? columns : firstCol + kcol;
                final int firstRow = j * krow;
                final int lastRow = (j == np - 1) ? rows : firstRow + krow;
                final int firstSlice = j * ksls;
                final int lastSlice = (j == np - 1) ? slices : firstSlice + ksls;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {

                    public void run() {
                        for (int c = firstCol; c < lastCol; c++) {
                            int cShifted = c;
                            if (cShifted > columns / 2)
                                cShifted = columns - cShifted;
                            float tmp = (cShifted / cc);
                            weights[0][c] = (float) Math.exp(-tmp * tmp);
                        }
                        for (int r = firstRow; r < lastRow; r++) {
                            int rShifted = r;
                            if (rShifted > rows / 2)
                                rShifted = rows - rShifted;
                            float tmp = (rShifted / rc);
                            weights[1][r] = (float) Math.exp(-tmp * tmp);
                        }
                        for (int s = firstSlice; s < lastSlice; s++) {
                            int sShifted = s;
                            if (sShifted > slices / 2)
                                sShifted = slices - sShifted;
                            float tmp = (sShifted / sc);
                            weights[2][s] = (float) Math.exp(-tmp * tmp);
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            for (int c = 0; c < columns; c++) {
                int cShifted = c;
                if (cShifted > columns / 2)
                    cShifted = columns - cShifted;
                float tmp = (cShifted / cc);
                weights[0][c] = (float) Math.exp(-tmp * tmp);
            }
            for (int r = 0; r < rows; r++) {
                int rShifted = r;
                if (rShifted > rows / 2)
                    rShifted = rows - rShifted;
                float tmp = (rShifted / rc);
                weights[1][r] = (float) Math.exp(-tmp * tmp);
            }
            for (int s = 0; s < slices; s++) {
                int sShifted = s;
                if (sShifted > slices / 2)
                    sShifted = slices - sShifted;
                float tmp = (sShifted / sc);
                weights[2][s] = (float) Math.exp(-tmp * tmp);
            }
        }
        return weights;
    }

    private static float meanDelta(FloatMatrix3D B, FloatMatrix3D AX, FloatMatrix3D X, final float aSum) {
        float meanDelta = 0;
        final float[] elemsB = (float[]) B.elements();
        final float[] elemsAX = (float[]) AX.elements();
        final float[] elemsX = (float[]) X.elements();
        final int size = (int)B.size();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            Float[] results = new Float[np];
            int k = size / np;
            for (int j = 0; j < np; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == np - 1) ? size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(new Callable<Float>() {
                    public Float call() throws Exception {
                        float meanDelta = 0;
                        float delta;
                        for (int i = firstIdx; i < lastIdx; i++) {
                            delta = (elemsB[i] - elemsAX[i] / aSum);
                            elemsX[i] += delta;
                            if (elemsX[i] < 0) {
                                elemsX[i] = 0;
                            } else {
                                meanDelta += Math.abs(delta);
                            }
                        }
                        return meanDelta;
                    }
                });
            }
            try {
                for (int j = 0; j < np; j++) {
                    results[j] = (Float) futures[j].get();
                }
            } catch (ExecutionException ex) {
                ex.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            for (int j = 0; j < np; j++) {
                meanDelta += results[j];
            }
        } else {
            float delta;
            for (int i = 0; i < size; i++) {
                delta = (elemsB[i] - elemsAX[i] / aSum);
                elemsX[i] += delta;
                if (elemsX[i] < 0) {
                    elemsX[i] = 0;
                } else {
                    meanDelta += Math.abs(delta);
                }
            }
        }
        return meanDelta;
    }

    private static void toDB(FloatMatrix3D X, final float minDB) {
        final float[] x = (float[]) X.elements();
        final float SCALE = (float) (10 / Math.log(10));
        final float minVal = (float) Math.exp(minDB / SCALE);
        int size = (int)X.size();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            int k = size / np;
            for (int j = 0; j < np; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == np - 1) ? size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {
                    public void run() {
                        for (int i = firstIdx; i < lastIdx; i++) {
                            if (x[i] > minVal)
                                x[i] = (float) (SCALE * Math.log(x[i]));
                            else
                                x[i] = minDB;
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            for (int i = 0; i < size; i++) {
                if (x[i] > minVal)
                    x[i] = (float) (SCALE * Math.log(x[i]));
                else
                    x[i] = minDB;
            }
        }
    }

    private static float unDB(FloatMatrix3D X) {
        final float[] x = (float[]) X.elements();
        final float SCALE = (float) (10 / Math.log(10));
        final int size = (int)X.size();
        float min = Float.MAX_VALUE;
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            Float[] results = new Float[np];
            int k = size / np;
            for (int j = 0; j < np; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == np - 1) ? size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(new Callable<Float>() {
                    public Float call() throws Exception {
                        float min = Float.MAX_VALUE;
                        for (int i = firstIdx; i < lastIdx; i++) {
                            if (x[i] < min)
                                min = x[i];
                            x[i] = (float) Math.exp(x[i] / SCALE);
                        }
                        return min;
                    }
                });
            }
            try {
                for (int j = 0; j < np; j++) {
                    results[j] = (Float) futures[j].get();
                }
            } catch (ExecutionException ex) {
                ex.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            min = results[0];
            for (int j = 1; j < np; j++) {
                if (results[j] < min)
                    min = results[j];
            }
        } else {
            for (int i = 0; i < size; i++) {
                if (x[i] < min)
                    min = x[i];
                x[i] = (float) Math.exp(x[i] / SCALE);
            }
        }
        return min;
    }

}
