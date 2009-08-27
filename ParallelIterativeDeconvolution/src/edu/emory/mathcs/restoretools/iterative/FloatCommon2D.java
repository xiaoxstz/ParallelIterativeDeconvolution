/*
 *  Copyright (C) 2008-2009 Piotr Wendykier
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.emory.mathcs.restoretools.iterative;

import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;

import java.util.concurrent.Future;

import cern.colt.matrix.tfloat.FloatMatrix1D;
import cern.colt.matrix.tfloat.FloatMatrix2D;
import cern.colt.matrix.tfloat.impl.DenseFloatMatrix2D;
import cern.jet.math.tfloat.FloatFunctions;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PaddingType;
import edu.emory.mathcs.utils.ConcurrencyUtils;

/**
 * Common methods for 2D iterative deblurring. Some code is from Bob Dougherty's Iterative
 * Deconvolve 3D.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 * 
 */
public class FloatCommon2D {

    private FloatCommon2D() {
    }

    /**
     * Tolerance for optimization.Fmin.fmin.
     */
    public static final float FMIN_TOL = 1.0e-4f;

    /**
     * Relative accuracy for float.
     */
    public static final float eps = (float) Math.pow(2, -23);

    /**
     * Square root of the relative accuracy for float.
     */
    public static final float sqrteps = (float) Math.sqrt(eps);

    /**
     * Copies pixel values from image processor <code>ip</code> to matrix
     * <code>X</code>.
     * 
     * @param ip
     *            image processor
     * @return matrix
     * 
     */
    public static FloatMatrix2D assignPixelsToMatrix(final ImageProcessor ip) {
        FloatMatrix2D X;
        if (ip instanceof FloatProcessor) {
            X = new DenseFloatMatrix2D(ip.getHeight(), ip.getWidth(), ((float[]) ip.getPixels()).clone(), 0, 0, ip.getWidth(), 1, false);
        } else {
            X = new DenseFloatMatrix2D(ip.getHeight(), ip.getWidth(), (float[]) ip.convertToFloat().getPixels(), 0, 0, ip.getWidth(), 1, false);
        }
        return X;
    }

    /**
     * Copies pixel values from image processor <code>ip</code> to matrix
     * <code>X</code>.
     * 
     * @param X
     *            matrix
     * @param ip
     *            image processor
     */
    public static void assignPixelsToMatrix(final FloatMatrix2D X, final ImageProcessor ip) {
        if (ip instanceof FloatProcessor) {
            X.assign((float[]) ip.getPixels());
        } else {
            X.assign((float[]) ip.convertToFloat().getPixels());
        }
    }

    /**
     * Copies pixel values from matrix <code>x</code> to image processor
     * <code>ip</code>.
     * 
     * @param ip
     *            image processor
     * @param x
     *            matrix
     * @param cmY
     *            color model
     */
    public static void assignPixelsToProcessor(final FloatProcessor ip, final FloatMatrix1D x, final java.awt.image.ColorModel cmY) {
        if (x.isView()) {
            ip.setPixels((float[]) x.copy().elements());
        } else {
            ip.setPixels((float[]) x.elements());
        }
        ip.setMinAndMax(0, 0);
        ip.setColorModel(cmY);
    }

    /**
     * Copies pixel values from matrix <code>x</code> to image processor
     * <code>ip</code>.
     * 
     * @param ip
     *            image processor
     * @param x
     *            matrix
     * @param cmY
     *            color model
     * @param threshold
     *            the smallest positive value assigned to the image processor,
     *            all the values less than the threshold are set to zero
     * 
     */
    public static void assignPixelsToProcessor(final FloatProcessor ip, final FloatMatrix1D x, final java.awt.image.ColorModel cmY, final float threshold) {
        final int rows = ip.getHeight();
        final int cols = ip.getWidth();
        final int size = rows * cols;
        final float[] px = (float[]) ip.getPixels();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if (x.isView()) {
            if ((np > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = size / np;
                for (int j = 0; j < np; j++) {
                    final int firstIdx = j * k;
                    final int lastIdx = (j == np - 1) ? size : firstIdx + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            float elem;
                            for (int i = firstIdx; i < lastIdx; i++) {
                                elem = (float) x.getQuick(i);
                                if (elem >= threshold) {
                                    px[i] = elem;
                                } else {
                                    px[i] = 0;
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                float elem;
                for (int i = 0; i < size; i++) {
                    elem = (float) x.getQuick(i);
                    if (elem >= threshold) {
                        px[i] = elem;
                    } else {
                        px[i] = 0;
                    }
                }
            }
        } else {
            final float[] elems = (float[]) x.elements();
            if ((np > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = size / np;
                for (int j = 0; j < np; j++) {
                    final int firstIdx = j * k;
                    final int lastIdx = (j == np - 1) ? size : firstIdx + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            float elem;
                            for (int i = firstIdx; i < lastIdx; i++) {
                                elem = (float) elems[i];
                                if (elem >= threshold) {
                                    px[i] = (float) elem;
                                } else {
                                    px[i] = 0;
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                float elem;
                for (int i = 0; i < size; i++) {
                    elem = (float) elems[i];
                    if (elem >= threshold) {
                        px[i] = (float) elem;
                    } else {
                        px[i] = 0;
                    }
                }
            }
        }
        ip.setMinAndMax(0, 0);
        ip.setColorModel(cmY);
    }

    /**
     * Copies pixel values from matrix <code>X</code> to image processor
     * <code>ip</code>.
     * 
     * @param ip
     *            image processor
     * @param X
     *            matrix
     * @param cmY
     *            color model
     * 
     */
    public static void assignPixelsToProcessor(final FloatProcessor ip, final FloatMatrix2D X, final java.awt.image.ColorModel cmY) {
        final int rows = X.rows();
        final int cols = X.columns();
        final float[] px = (float[]) ip.getPixels();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if (X.isView()) {
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rows / np;
                for (int j = 0; j < np; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == np - 1) ? rows : firstRow + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int c = 0; c < cols; c++) {
                                    px[c + cols * r] = (float) X.getQuick(r, c);
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < cols; c++) {
                        px[c + cols * r] = (float) X.getQuick(r, c);
                    }
                }
            }
        } else {
            final float[] elems = (float[]) X.elements();
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rows / np;
                for (int j = 0; j < np; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == np - 1) ? rows : firstRow + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            int idx = firstRow * cols;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int c = 0; c < cols; c++) {
                                    px[idx] = (float) elems[idx];
                                    idx++;
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int idx = 0;
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < cols; c++) {
                        px[idx] = (float) elems[idx];
                        idx++;
                    }
                }
            }
        }
        ip.setMinAndMax(0, 0);
        ip.setColorModel(cmY);
    }

    /**
     * Copies pixel values from complex matrix <code>X</code> to image processor
     * <code>ip</code>
     * 
     * @param ip
     *            image processor
     * @param X
     *            padded matrix
     * @param rows
     *            original number of rows
     * @param cols
     *            original number of columns
     * @param rOff
     *            row offset
     * @param cOff
     *            column offset
     * @param cmY
     *            color model
     * 
     */
    public static void assignPixelsToProcessorPadded(final FloatProcessor ip, final FloatMatrix2D X, final int rows, final int cols, final int rOff, final int cOff, final java.awt.image.ColorModel cmY) {
        final float[] px = (float[]) ip.getPixels();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if (X.isView()) {
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rows / np;
                for (int j = 0; j < np; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == np - 1) ? rows : firstRow + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int c = 0; c < cols; c++) {
                                    px[c + cols * r] = (float) X.getQuick(r + rOff, c + cOff);
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < cols; c++) {
                        px[c + cols * r] = (float) X.getQuick(r + rOff, c + cOff);
                    }
                }
            }
        } else {
            final float[] elems = (float[]) X.elements();
            final int rowStride = X.columns();
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rows / np;
                for (int j = 0; j < np; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == np - 1) ? rows : firstRow + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            int idx;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int c = 0; c < cols; c++) {
                                    idx = (r + rOff) * rowStride + (c + cOff);
                                    px[r * cols + c] = (float) elems[idx];
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int idx;
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < cols; c++) {
                        idx = (r + rOff) * rowStride + (c + cOff);
                        px[r * cols + c] = (float) elems[idx];
                        idx++;
                    }
                }
            }
        }
        ip.setMinAndMax(0, 0);
        ip.setColorModel(cmY);
    }

    /**
     * Copies pixel values from complex matrix <code>X</code> to image processor
     * <code>ip</code>
     * 
     * @param ip
     *            image processor
     * @param X
     *            padded matrix
     * @param rows
     *            original number of rows
     * @param cols
     *            original number of columns
     * @param rOff
     *            row offset
     * @param cOff
     *            column offset
     * @param cmY
     *            color model
     * @param threshold
     *            the smallest positive value assigned to the image processor,
     *            all the values less than the threshold are set to zero
     */
    public static void assignPixelsToProcessorPadded(final FloatProcessor ip, final FloatMatrix2D X, final int rows, final int cols, final int rOff, final int cOff, final java.awt.image.ColorModel cmY, final float threshold) {
        final float[] px = (float[]) ip.getPixels();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if (X.isView()) {
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rows / np;
                for (int j = 0; j < np; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == np - 1) ? rows : firstRow + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            float elem;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int c = 0; c < cols; c++) {
                                    elem = (float) X.getQuick(r + rOff, c + cOff);
                                    if (elem >= threshold) {
                                        px[c + cols * r] = elem;
                                    } else {
                                        px[c + cols * r] = 0;
                                    }
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                float elem;
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < cols; c++) {
                        elem = (float) X.getQuick(r + rOff, c + cOff);
                        if (elem >= threshold) {
                            px[c + cols * r] = elem;
                        } else {
                            px[c + cols * r] = 0;
                        }
                    }
                }
            }
        } else {
            final float[] elems = (float[]) X.elements();
            final int rowStride = X.columns();
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rows / np;
                for (int j = 0; j < np; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == np - 1) ? rows : firstRow + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            int idx;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int c = 0; c < cols; c++) {
                                    idx = (r + rOff) * rowStride + (c + cOff);
                                    if ((float) elems[idx] >= threshold) {
                                        px[r * cols + c] = (float) elems[idx];
                                    } else {
                                        px[r * cols + c] = 0;
                                    }
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int idx;
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < cols; c++) {
                        idx = (r + rOff) * rowStride + (c + cOff);
                        if ((float) elems[idx] >= threshold) {
                            px[r * cols + c] = (float) elems[idx];
                        } else {
                            px[r * cols + c] = 0;
                        }
                    }
                }
            }
        }
        ip.setMinAndMax(0, 0);
        ip.setColorModel(cmY);
    }

    /**
     * Copies pixel values from matrix <code>X</code> to image processor
     * <code>ip</code>.
     * 
     * @param ip
     *            image processor
     * @param X
     *            matrix
     * @param cmY
     *            color model
     * @param threshold
     *            the smallest positive value assigned to the image processor,
     *            all the values less than the threshold are set to zero
     * 
     */
    public static void assignPixelsToProcessor(final FloatProcessor ip, final FloatMatrix2D X, final java.awt.image.ColorModel cmY, final float threshold) {
        final int rows = X.rows();
        final int cols = X.columns();
        final float[] px = (float[]) ip.getPixels();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if (X.isView()) {
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rows / np;
                for (int j = 0; j < np; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == np - 1) ? rows : firstRow + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            float elem;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int c = 0; c < cols; c++) {
                                    elem = (float) X.getQuick(r, c);
                                    if (elem >= threshold) {
                                        px[c + cols * r] = elem;
                                    } else {
                                        px[c + cols * r] = 0;
                                    }
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                float elem;
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < cols; c++) {
                        elem = (float) X.getQuick(r, c);
                        if (elem >= threshold) {
                            px[c + cols * r] = elem;
                        } else {
                            px[c + cols * r] = 0;
                        }
                    }
                }
            }
        } else {
            final float[] elems = (float[]) X.elements();
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rows / np;
                for (int j = 0; j < np; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == np - 1) ? rows : firstRow + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            int idx = firstRow * cols;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int c = 0; c < cols; c++) {
                                    if ((float) elems[idx] >= threshold) {
                                        px[idx] = (float) elems[idx];
                                    } else {
                                        px[idx] = 0;
                                    }
                                    idx++;
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int idx = 0;
                for (int r = 0; r < rows; r++) {
                    for (int c = 0; c < cols; c++) {
                        if ((float) elems[idx] >= threshold) {
                            px[idx] = (float) elems[idx];
                        } else {
                            px[idx] = 0;
                        }
                        idx++;
                    }
                }
            }
        }
        ip.setMinAndMax(0, 0);
        ip.setColorModel(cmY);
    }

    /**
     * Computes the circular shift of <code>PSF</code> matrix. This method
     * computes a matrix containing first column of a blurring matrix when
     * implementing periodic boundary conditions.
     * 
     * @param PSF
     *            real matrix containing the point spread function.
     * @param center
     *            indices of center of <code>PSF</code>.
     * @return real matrix containing first column of a blurring matrix
     */
    public static FloatMatrix2D circShift(FloatMatrix2D PSF, int[] center) {

        int rows = PSF.rows();
        int cols = PSF.columns();
        int cr = center[0];
        int cc = center[1];
        FloatMatrix2D P1 = new DenseFloatMatrix2D(rows, cols);
        P1.viewPart(0, 0, rows - cr, cols - cc).assign(PSF.viewPart(cr, cc, rows - cr, cols - cc));
        P1.viewPart(0, cols - cc, rows - cr, cc).assign(PSF.viewPart(cr, 0, rows - cr, cc));
        P1.viewPart(rows - cr, 0, cr, cols - cc).assign(PSF.viewPart(0, cc, cr, cols - cc));
        P1.viewPart(rows - cr, cols - cc, cr, cc).assign(PSF.viewPart(0, 0, cr, cc));
        return P1;
    }

    /**
     * Converts an image into a given output type.
     * 
     * @param image
     *            image
     * @param output
     *            output type
     */
    public static void convertImage(ImagePlus image, OutputType output) {
        switch (output) {
        case BYTE:
            new ImageConverter(image).convertToGray8();
            break;
        case SHORT:
            new ImageConverter(image).convertToGray16();
            break;
        case FLOAT:
            //image is always in 32-bit precision
            break;
        }
    }

    /**
     * Pads matrix <code>X</code> with periodic boundary conditions.
     * 
     * @param X
     *            matrix to be padded
     * @param rowsPad
     *            number of rows in padded matrix
     * @param colsPad
     *            number of columns in padded matrix
     * @return padded matrix
     */
    public static FloatMatrix2D padPeriodic(final FloatMatrix2D X, final int rowsPad, final int colsPad) {
        final int rows = X.rows();
        final int cols = X.columns();
        if ((rows == rowsPad) && (cols == colsPad)) {
            return X;
        }
        final FloatMatrix2D Xpad = new DenseFloatMatrix2D(rowsPad, colsPad);
        final int rOff = (rowsPad - rows + 1) / 2;
        final int cOff = (colsPad - cols + 1) / 2;
        final float[] elemsXpad = (float[]) Xpad.elements();
        final int rowStrideXpad = Xpad.columns();
        int np = ConcurrencyUtils.getNumberOfThreads();

        if (X.isView()) {
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rowsPad / np;
                for (int j = 0; j < np; j++) {
                    final int firstIdx = -rOff + j * k;
                    final int lastIdx = (j == np - 1) ? rowsPad - rOff : firstIdx + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            int cIn, rIn, cOut, rOut;
                            int idxXpad;
                            for (int r = firstIdx; r < lastIdx; r++) {
                                rOut = r + rOff;
                                rIn = periodic(r, rows);
                                for (int c = -cOff; c < colsPad - cOff; c++) {
                                    cOut = c + cOff;
                                    cIn = periodic(c, cols);
                                    idxXpad = rOut * rowStrideXpad + cOut;
                                    elemsXpad[idxXpad] = X.getQuick(rIn, cIn);
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int cIn, rIn, cOut, rOut;
                int idxXpad;
                for (int r = -rOff; r < rowsPad - rOff; r++) {
                    rOut = r + rOff;
                    rIn = periodic(r, rows);
                    for (int c = -cOff; c < colsPad - cOff; c++) {
                        cOut = c + cOff;
                        cIn = periodic(c, cols);
                        idxXpad = rOut * rowStrideXpad + cOut;
                        elemsXpad[idxXpad] = X.getQuick(rIn, cIn);
                    }
                }
            }
        } else {
            final float[] elemsX = (float[]) X.elements();
            final int rowStrideX = X.columns();
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rowsPad / np;
                for (int j = 0; j < np; j++) {
                    final int firstIdx = -rOff + j * k;
                    final int lastIdx = (j == np - 1) ? rowsPad - rOff : firstIdx + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            int cIn, rIn, cOut, rOut;
                            int idxX;
                            int idxXpad;
                            for (int r = firstIdx; r < lastIdx; r++) {
                                rOut = r + rOff;
                                rIn = periodic(r, rows);
                                for (int c = -cOff; c < colsPad - cOff; c++) {
                                    cOut = c + cOff;
                                    cIn = periodic(c, cols);
                                    idxX = rIn * rowStrideX + cIn;
                                    idxXpad = rOut * rowStrideXpad + cOut;
                                    elemsXpad[idxXpad] = elemsX[idxX];
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int cIn, rIn, cOut, rOut;
                int idxXpad;
                int idxX;
                for (int r = -rOff; r < rowsPad - rOff; r++) {
                    rOut = r + rOff;
                    rIn = periodic(r, rows);
                    for (int c = -cOff; c < colsPad - cOff; c++) {
                        cOut = c + cOff;
                        cIn = periodic(c, cols);
                        idxX = rIn * rowStrideX + cIn;
                        idxXpad = rOut * rowStrideXpad + cOut;
                        elemsXpad[idxXpad] = elemsX[idxX];
                    }
                }
            }
        }
        return Xpad;
    }

    /**
     * Pads matrix <code>X</code> with reflexive boundary conditions.
     * 
     * @param X
     *            matrix to be padded
     * @param rowsPad
     *            number of rows in padded matrix
     * @param colsPad
     *            number of columns in padded matrix
     * @return padded matrix
     */
    public static FloatMatrix2D padReflexive(final FloatMatrix2D X, final int rowsPad, final int colsPad) {
        final int rows = X.rows();
        final int cols = X.columns();
        if ((rows == rowsPad) && (cols == colsPad)) {
            return X;
        }
        final FloatMatrix2D Xpad = new DenseFloatMatrix2D(rowsPad, colsPad);
        final int rOff = (rowsPad - rows + 1) / 2;
        final int cOff = (colsPad - cols + 1) / 2;
        final float[] elemsXpad = (float[]) Xpad.elements();
        final int rowStrideXpad = Xpad.columns();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if (X.isView()) {
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rowsPad / np;
                for (int j = 0; j < np; j++) {
                    final int firstIdx = -rOff + j * k;
                    final int lastIdx = (j == np - 1) ? rowsPad - rOff : firstIdx + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            int cIn, rIn, cOut, rOut;
                            int idxXpad;
                            for (int r = firstIdx; r < lastIdx; r++) {
                                rOut = r + rOff;
                                rIn = mirror(r, rows);
                                for (int c = -cOff; c < colsPad - cOff; c++) {
                                    cOut = c + cOff;
                                    cIn = mirror(c, cols);
                                    idxXpad = rOut * rowStrideXpad + cOut;
                                    elemsXpad[idxXpad] = X.getQuick(rIn, cIn);
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int cIn, rIn, cOut, rOut;
                int idxXpad;
                for (int r = -rOff; r < rowsPad - rOff; r++) {
                    rOut = r + rOff;
                    rIn = mirror(r, rows);
                    for (int c = -cOff; c < colsPad - cOff; c++) {
                        cOut = c + cOff;
                        cIn = mirror(c, cols);
                        idxXpad = rOut * rowStrideXpad + cOut;
                        elemsXpad[idxXpad] = X.getQuick(rIn, cIn);
                    }
                }
            }
        } else {
            final float[] elemsX = (float[]) X.elements();
            final int rowStrideX = X.columns();
            if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                Future<?>[] futures = new Future[np];
                int k = rowsPad / np;
                for (int j = 0; j < np; j++) {
                    final int firstIdx = -rOff + j * k;
                    final int lastIdx = (j == np - 1) ? rowsPad - rOff : firstIdx + k;
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {
                        public void run() {
                            int cIn, rIn, cOut, rOut;
                            int idxX;
                            int idxXpad;
                            for (int r = firstIdx; r < lastIdx; r++) {
                                rOut = r + rOff;
                                rIn = mirror(r, rows);
                                for (int c = -cOff; c < colsPad - cOff; c++) {
                                    cOut = c + cOff;
                                    cIn = mirror(c, cols);
                                    idxX = rIn * rowStrideX + cIn;
                                    idxXpad = rOut * rowStrideXpad + cOut;
                                    elemsXpad[idxXpad] = elemsX[idxX];
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int cIn, rIn, cOut, rOut;
                int idxX;
                int idxXpad;
                for (int r = -rOff; r < rowsPad - rOff; r++) {
                    rOut = r + rOff;
                    rIn = mirror(r, rows);
                    for (int c = -cOff; c < colsPad - cOff; c++) {
                        cOut = c + cOff;
                        cIn = mirror(c, cols);
                        idxX = rIn * rowStrideX + cIn;
                        idxXpad = rOut * rowStrideXpad + cOut;
                        elemsXpad[idxXpad] = elemsX[idxX];
                    }
                }
            }

        }
        return Xpad;
    }

    /**
     * Pads matrix <code>X</code> with zero boundary conditions.
     * 
     * @param X
     *            matrix to be padded
     * @param rowsPad
     *            number of rows in padded matrix
     * @param colsPad
     *            number of columns in padded matrix
     * @return padded matrix
     */
    public static FloatMatrix2D padZero(FloatMatrix2D X, int rowsPad, int colsPad) {
        final int rows = X.rows();
        final int cols = X.columns();
        if ((rows == rowsPad) && (cols == colsPad)) {
            return X;
        }
        FloatMatrix2D Xpad = new DenseFloatMatrix2D(rowsPad, colsPad);
        final int rOff = (rowsPad - rows + 1) / 2;
        final int cOff = (colsPad - cols + 1) / 2;
        Xpad.viewPart(rOff, cOff, rows, cols).assign(X);
        return Xpad;
    }

    /**
     * Pads matrix <code>X</code> with zero boundary conditions.
     * 
     * @param X
     *            matrix to be padded
     * @param padSize
     *            padding size
     * @param padding
     *            type of padding
     * @return padded matrix
     */
    public static FloatMatrix2D padZero(FloatMatrix2D X, int[] padSize, PaddingType padding) {
        if ((padSize[0] == 0) && (padSize[1] == 0)) {
            return X;
        }
        FloatMatrix2D Xpad = null;
        switch (padding) {
        case BOTH:
            Xpad = new DenseFloatMatrix2D(X.rows() + 2 * padSize[0], X.columns() + 2 * padSize[1]);
            Xpad.viewPart(padSize[0], padSize[1], X.rows(), X.columns()).assign(X);
            break;
        case POST:
            Xpad = new DenseFloatMatrix2D(X.rows() + padSize[0], X.columns() + padSize[1]);
            Xpad.viewPart(0, 0, X.rows(), X.columns()).assign(X);
            break;
        case PRE:
            Xpad = new DenseFloatMatrix2D(X.rows() + padSize[0], X.columns() + padSize[1]);
            Xpad.viewPart(padSize[0], padSize[1], X.rows(), X.columns()).assign(X);
            break;
        }
        return Xpad;
    }

    public static FloatMatrix2D dctShift(FloatMatrix2D PSF, int[] center) {
        int rows = PSF.rows();
        int cols = PSF.columns();
        int cr = center[0];
        int cc = center[1];
        int k = Math.min(Math.min(Math.min(cr, rows - cr - 1), cc), cols - cc - 1);
        int frow = cr - k;
        int lrow = cr + k;
        int rowSize = lrow - frow + 1;
        int fcol = cc - k;
        int lcol = cc + k;
        int colSize = lcol - fcol + 1;

        FloatMatrix2D PP = new DenseFloatMatrix2D(rowSize, colSize);
        FloatMatrix2D P1 = new DenseFloatMatrix2D(rowSize, colSize);
        FloatMatrix2D P2 = new DenseFloatMatrix2D(rowSize, colSize);
        FloatMatrix2D P3 = new DenseFloatMatrix2D(rowSize, colSize);
        FloatMatrix2D P4 = new DenseFloatMatrix2D(rowSize, colSize);
        FloatMatrix2D Ps = new DenseFloatMatrix2D(rows, cols);
        PP.assign(PSF.viewPart(frow, fcol, rowSize, colSize));

        P1.viewPart(0, 0, rowSize - cr + frow, colSize - cc + fcol).assign(PP.viewPart(cr - frow, cc - fcol, rowSize - cr + frow, colSize - cc + fcol));
        P2.viewPart(0, 0, rowSize - cr + frow, colSize - cc + fcol - 1).assign(PP.viewPart(cr - frow, cc - fcol + 1, rowSize - cr + frow, colSize - cc + fcol - 1));
        P3.viewPart(0, 0, rowSize - cr + frow - 1, colSize - cc + fcol).assign(PP.viewPart(cr - frow + 1, cc - fcol, rowSize - cr + frow - 1, colSize - cc + fcol));
        P4.viewPart(0, 0, rowSize - cr + frow - 1, colSize - cc + fcol - 1).assign(PP.viewPart(cr - frow + 1, cc - fcol + 1, rowSize - cr + frow - 1, colSize - cc + fcol - 1));
        P1.assign(P2, FloatFunctions.plus);
        P1.assign(P3, FloatFunctions.plus);
        P1.assign(P4, FloatFunctions.plus);
        Ps.viewPart(0, 0, 2 * k + 1, 2 * k + 1).assign(P1);
        return Ps.copy();
    }

    public static void swapQuadrants(final int rows, final int columns, FloatMatrix2D X) {
        final float[] x = (float[]) X.elements();
        final int cHalf = columns / 2;
        final int rHalf = rows / 2;
        int size = columns * rows / 2;
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            Future<?>[] futures = new Future[np];
            int krow = rHalf / np;
            for (int j = 0; j < np; j++) {
                final int firstRow = j * krow;
                final int lastRow = (j == np - 1) ? rHalf : firstRow + krow;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {

                    public void run() {
                        int rP, idx1, idx2;
                        float temp;
                        for (int r = firstRow; r < lastRow; r++) {
                            rP = r + rHalf;
                            for (int c = 0; c < columns; c++) {
                                idx1 = c + columns * r;
                                idx2 = c + columns * rP;
                                temp = x[idx1];
                                x[idx1] = x[idx2];
                                x[idx2] = temp;
                            }
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int rP, idx1, idx2;
            float temp;
            for (int r = 0; r < rHalf; r++) {
                rP = r + rHalf;
                for (int c = 0; c < columns; c++) {
                    idx1 = c + columns * r;
                    idx2 = c + columns * rP;
                    temp = x[idx1];
                    x[idx1] = x[idx2];
                    x[idx2] = temp;
                }
            }
        }
        if ((np > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            Future<?>[] futures = new Future[np];
            int kcol = cHalf / np;
            for (int j = 0; j < np; j++) {
                final int firstCol = j * kcol;
                final int lastCol = (j == np - 1) ? cHalf : firstCol + kcol;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {
                    public void run() {
                        int cP, idx1, idx2;
                        float temp;
                        for (int c = firstCol; c < lastCol; c++) {
                            cP = c + cHalf;
                            for (int r = 0; r < rows; r++) {
                                idx1 = c + columns * r;
                                idx2 = cP + columns * r;
                                temp = x[idx1];
                                x[idx1] = x[idx2];
                                x[idx2] = temp;
                            }
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int cP, idx1, idx2;
            float temp;
            for (int c = 0; c < cHalf; c++) {
                cP = c + cHalf;
                for (int r = 0; r < rows; r++) {
                    idx1 = c + columns * r;
                    idx2 = cP + columns * r;
                    temp = x[idx1];
                    x[idx1] = x[idx2];
                    x[idx2] = temp;
                }
            }
        }
    }

    public static void convolveTransposeFD(FloatMatrix2D H1, FloatMatrix2D H2, FloatMatrix2D Result) {
        final int rows = H1.rows();
        final int columns = H1.columns();
        final float[] h1 = (float[]) H1.elements();
        final float[] h2 = (float[]) H2.elements();
        final float[] result = (float[]) Result.elements();
        int cC, rC, idx1, idx2;
        float h2e, h2o;
        for (int r = 0; r < rows; r++) {
            rC = (rows - r) % rows;
            for (int c = 0; c < columns; c++) {
                cC = (columns - c) % columns;
                idx1 = c + columns * r;
                idx2 = cC + columns * rC;
                h2e = (h2[idx1] + h2[idx2]) / 2;
                h2o = (h2[idx1] - h2[idx2]) / 2;
                result[idx2] = (float) (h1[idx1] * h2e - h1[idx2] * h2o);
            }
        }
    }

    public static void convolveFD(FloatMatrix2D H1, FloatMatrix2D H2, FloatMatrix2D Result) {
        final int rows = H1.rows();
        final int columns = H1.columns();
        final float[] h1 = (float[]) H1.elements();
        final float[] h2 = (float[]) H2.elements();
        final float[] result = (float[]) Result.elements();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (columns * rows >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            Future<?>[] futures = new Future[np];
            int k = rows / np;
            for (int j = 0; j < np; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == np - 1) ? rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(new Runnable() {
                    public void run() {
                        int cC, rC, idx1, idx2;
                        float h2e, h2o;
                        for (int r = firstRow; r < lastRow; r++) {
                            rC = (rows - r) % rows;
                            for (int c = 0; c < columns; c++) {
                                cC = (columns - c) % columns;
                                idx1 = c + columns * r;
                                idx2 = cC + columns * rC;
                                h2e = (h2[idx1] + h2[idx2]) / 2;
                                h2o = (h2[idx1] - h2[idx2]) / 2;
                                result[idx1] = (float) (h1[idx1] * h2e + h1[idx2] * h2o);
                            }
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int cC, rC, idx1, idx2;
            float h2e, h2o;
            for (int r = 0; r < rows; r++) {
                rC = (rows - r) % rows;
                for (int c = 0; c < columns; c++) {
                    cC = (columns - c) % columns;
                    idx1 = c + columns * r;
                    idx2 = cC + columns * rC;
                    h2e = (h2[idx1] + h2[idx2]) / 2;
                    h2o = (h2[idx1] - h2[idx2]) / 2;
                    result[idx1] = (float) (h1[idx1] * h2e + h1[idx2] * h2o);
                }
            }
        }
    }

    private static int mirror(int i, int n) {
        int ip = mod(i, 2 * n);
        if (ip < n) {
            return ip;
        } else {
            return n - (ip % n) - 1;
        }
    }

    private static int mod(int i, int n) {
        return ((i % n) + n) % n;
    }

    private static int periodic(int i, int n) {
        int ip = mod(i, 2 * n);
        if (ip < n) {
            return ip;
        } else {
            return (ip % n);
        }
    }

}
