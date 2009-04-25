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
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.StackConverter;

import java.util.concurrent.Future;

import cern.colt.matrix.tfloat.FloatMatrix1D;
import cern.colt.matrix.tfloat.FloatMatrix3D;
import cern.colt.matrix.tfloat.impl.DenseFloatMatrix3D;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PaddingType;
import edu.emory.mathcs.utils.ConcurrencyUtils;

/**
 * Common methods for 3D iterative deblurring. Some code is from Bob Dougherty's Iterative
 * Deconvolve 3D.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class FloatCommon3D {

    private FloatCommon3D() {
    }

    /**
     * Copies pixel values from image stack <code>stack</code> to matrix
     * <code>X</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            matrix
     */
    public static void assignPixelsToMatrix(final ImageStack stack, final FloatMatrix3D X) {
        int slices = X.slices();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (X.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
            Future<?>[] futures = new Future[np];
            int k = slices / np;
            for (int j = 0; j < np; j++) {
                final int startslice = j * k;
                final int stopslice;
                if (j == np - 1) {
                    stopslice = slices;
                } else {
                    stopslice = startslice + k;
                }
                futures[j] = ConcurrencyUtils.submit(new Runnable() {
                    public void run() {
                        for (int s = startslice; s < stopslice; s++) {
                            ImageProcessor ip = stack.getProcessor(s + 1);
                            FloatCommon2D.assignPixelsToMatrix(X.viewSlice(s), ip);
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            for (int s = 0; s < slices; s++) {
                ImageProcessor ip = stack.getProcessor(s + 1);
                FloatCommon2D.assignPixelsToMatrix(X.viewSlice(s), ip);
            }
        }
    }

    /**
     * Copies pixel values from real matrix <code>X</code> to image stack
     * <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param cmY
     *            color model
     */
    public static void assignPixelsToStack(final ImageStack stack, final FloatMatrix3D X, final java.awt.image.ColorModel cmY) {
        final int slices = X.slices();
        final int rows = X.rows();
        final int cols = X.columns();
        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = new FloatProcessor(cols, rows);
            FloatCommon2D.assignPixelsToProcessor(ip, X.viewSlice(s), cmY);
            stack.addSlice(null, ip, s);
        }
    }

    /**
     * Copies pixel values from real matrix <code>X</code> to image stack
     * <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param cmY
     *            color model
     * @param threshold
     *            the smallest positive value assigned to the image processor,
     *            all the values less than the threshold are set to zero
     * 
     */
    public static void assignPixelsToStack(final ImageStack stack, final FloatMatrix3D X, final java.awt.image.ColorModel cmY, float threshold) {
        final int slices = X.slices();
        final int rows = X.rows();
        final int cols = X.columns();
        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = new FloatProcessor(cols, rows);
            FloatCommon2D.assignPixelsToProcessor(ip, X.viewSlice(s), cmY, threshold);
            stack.addSlice(null, ip, s);
        }
    }

    /**
     * Copies pixel values from real padded matrix <code>X</code> to image stack
     * <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            padded real matrix
     * @param slices
     *            original number of slices
     * @param rows
     *            original number of rows
     * @param cols
     *            original number of columns
     * @param sOff
     *            slice offset
     * @param rOff
     *            row offset
     * @param cOff
     *            column offset
     * @param cmY
     *            color model
     */
    public static void assignPixelsToStackPadded(final ImageStack stack, final FloatMatrix3D X, final int slices, final int rows, final int cols, final int sOff, final int rOff, final int cOff, final java.awt.image.ColorModel cmY) {
        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = new FloatProcessor(cols, rows);
            FloatCommon2D.assignPixelsToProcessorPadded(ip, X.viewSlice(s + sOff), rows, cols, rOff, cOff, cmY);
            stack.addSlice(null, ip, s);
        }
    }

    /**
     * Copies pixel values from real padded matrix <code>X</code> to image stack
     * <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            padded real matrix
     * @param slices
     *            original number of slices
     * @param rows
     *            original number of rows
     * @param cols
     *            original number of columns
     * @param sOff
     *            slice offset
     * @param rOff
     *            row offset
     * @param cOff
     *            column offset
     * @param cmY
     *            color model
     * @param threshold
     *            the smallest positive value assigned to the stack, all the
     *            values less than the threshold are set to zero
     * 
     */
    public static void assignPixelsToStackPadded(final ImageStack stack, final FloatMatrix3D X, final int slices, final int rows, final int cols, final int sOff, final int rOff, final int cOff, final java.awt.image.ColorModel cmY, final float threshold) {
        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = new FloatProcessor(cols, rows);
            FloatCommon2D.assignPixelsToProcessorPadded(ip, X.viewSlice(s + sOff), rows, cols, rOff, cOff, cmY, threshold);
            stack.addSlice(null, ip, s);
        }
    }

    /**
     * Updates pixel values in image stack <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param cmY
     *            color model
     * 
     */
    public static void updatePixelsInStack(final ImageStack stack, final FloatMatrix3D X, final java.awt.image.ColorModel cmY) {
        final int slices = X.slices();
        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = (FloatProcessor) stack.getProcessor(s + 1);
            FloatCommon2D.assignPixelsToProcessor(ip, X.viewSlice(s), cmY);
        }
    }

    /**
     * Updates pixel values in image stack <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param slices
     *            original number of slices
     * @param rows
     *            original number of rows
     * @param cols
     *            original number of columns
     * @param sOff
     *            slice offset
     * @param rOff
     *            row offset
     * @param cOff
     *            column offset
     * @param cmY
     *            color model
     * @param threshold
     *            the smallest positive value assigned to the image processor,
     *            all the values less than the threshold are set to zero
     * 
     */
    public static void updatePixelsInStackPadded(final ImageStack stack, final FloatMatrix3D X, final int slices, final int rows, final int cols, final int sOff, final int rOff, final int cOff, final java.awt.image.ColorModel cmY, float threshold) {
        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = (FloatProcessor) stack.getProcessor(s + 1);
            FloatCommon2D.assignPixelsToProcessorPadded(ip, X.viewSlice(s + sOff), rows, cols, rOff, cOff, cmY, threshold);
        }
    }

    /**
     * Updates pixel values in image stack <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param slices
     *            original number of slices
     * @param rows
     *            original number of rows
     * @param cols
     *            original number of columns
     * @param sOff
     *            slice offset
     * @param rOff
     *            row offset
     * @param cOff
     *            column offset
     * @param cmY
     *            color model
     */
    public static void updatePixelsInStackPadded(final ImageStack stack, final FloatMatrix3D X, final int slices, final int rows, final int cols, final int sOff, final int rOff, final int cOff, final java.awt.image.ColorModel cmY) {
        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = (FloatProcessor) stack.getProcessor(s + 1);
            FloatCommon2D.assignPixelsToProcessorPadded(ip, X.viewSlice(s + sOff), rows, cols, rOff, cOff, cmY);
        }
    }

    /**
     * Updates pixel values in image stack <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param cmY
     *            color model
     * @param threshold
     *            the smallest positive value assigned to the image processor,
     *            all the values less than the threshold are set to zero
     * 
     */
    public static void updatePixelsInStack(final ImageStack stack, final FloatMatrix3D X, final java.awt.image.ColorModel cmY, float threshold) {
        final int slices = X.slices();
        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = (FloatProcessor) stack.getProcessor(s + 1);
            FloatCommon2D.assignPixelsToProcessor(ip, X.viewSlice(s), cmY, threshold);
        }
    }

    /**
     * Copies pixel values from real matrix <code>X</code> to image stack
     * <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param xsize3D
     *            size of X as 3D matrix
     * @param cmY
     *            color model
     * 
     */
    public static void assignPixelsToStack(final ImageStack stack, final FloatMatrix1D X, int[] xsize3D, final java.awt.image.ColorModel cmY) {
        final int slices = xsize3D[0];
        final int rows = xsize3D[1];
        final int cols = xsize3D[2];
        final int sliceStride = rows * cols;

        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = new FloatProcessor(cols, rows);
            FloatCommon2D.assignPixelsToProcessor(ip, X.viewPart(s * sliceStride, sliceStride), cmY);
            ip.setColorModel(cmY);
            ip.setMinAndMax(0, 0);
            stack.addSlice(null, ip, s);
        }
    }

    /**
     * Copies pixel values from real matrix <code>X</code> to image stack
     * <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param xsize3D
     *            size of X as 3D matrix
     * @param cmY
     *            color model
     * @param threshold
     *            the smallest positive value assigned to the image processor,
     *            all the values less than the threshold are set to zero
     * 
     */
    public static void assignPixelsToStack(final ImageStack stack, final FloatMatrix1D X, int[] xsize3D, final java.awt.image.ColorModel cmY, float threshold) {
        final int slices = xsize3D[0];
        final int rows = xsize3D[1];
        final int cols = xsize3D[2];
        final int sliceStride = rows * cols;

        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = new FloatProcessor(cols, rows);
            FloatCommon2D.assignPixelsToProcessor(ip, X.viewPart(s * sliceStride, sliceStride), cmY, threshold);
            ip.setColorModel(cmY);
            ip.setMinAndMax(0, 0);
            stack.addSlice(null, ip, s);
        }
    }

    /**
     * Updates pixel values in image stack <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param xsize3D
     *            size of X as 3D matrix
     * @param cmY
     *            color model
     */
    public static void updatePixelsInStack(final ImageStack stack, final FloatMatrix1D X, int[] xsize3D, final java.awt.image.ColorModel cmY) {
        final int slices = xsize3D[0];
        final int rows = xsize3D[1];
        final int cols = xsize3D[2];
        final int sliceStride = rows * cols;

        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = (FloatProcessor) stack.getProcessor(s + 1);
            FloatCommon2D.assignPixelsToProcessor(ip, X.viewPart(s * sliceStride, sliceStride), cmY);
        }
    }

    /**
     * Updates pixel values in image stack <code>stack</code>
     * 
     * @param stack
     *            image stack
     * @param X
     *            real matrix
     * @param xsize3D
     *            size of X as 3D matrix
     * @param cmY
     *            color model
     * @param threshold
     *            the smallest positive value assigned to the image processor,
     *            all the values less than the threshold are set to zero
     * 
     */
    public static void updatePixelsInStack(final ImageStack stack, final FloatMatrix1D X, int[] xsize3D, final java.awt.image.ColorModel cmY, float threshold) {
        final int slices = xsize3D[0];
        final int rows = xsize3D[1];
        final int cols = xsize3D[2];
        final int sliceStride = rows * cols;

        for (int s = 0; s < slices; s++) {
            FloatProcessor ip = (FloatProcessor) stack.getProcessor(s + 1);
            FloatCommon2D.assignPixelsToProcessor(ip, X.viewPart(s * sliceStride, sliceStride), cmY, threshold);
        }
    }

    /**
     * Computes the circular shift of <code>PSF</code> matrix.
     * 
     * @param PSF
     *            real matrix containing the point spread function.
     * @param center
     *            indices of center of <code>PSF</code>.
     * @return shifted matrix
     */
    public static FloatMatrix3D circShift(FloatMatrix3D PSF, int[] center) {
        int slices = PSF.slices();
        int rows = PSF.rows();
        int cols = PSF.columns();
        int cs = center[0];
        int cr = center[1];
        int cc = center[2];
        FloatMatrix3D P1 = new DenseFloatMatrix3D(slices, rows, cols);

        P1.viewPart(0, 0, 0, slices - cs, rows - cr, cols - cc).assign(PSF.viewPart(cs, cr, cc, slices - cs, rows - cr, cols - cc));
        P1.viewPart(0, rows - cr, 0, slices - cs, cr, cols - cc).assign(PSF.viewPart(cs, 0, cc, slices - cs, cr, cols - cc));
        P1.viewPart(0, 0, cols - cc, slices - cs, rows - cr, cc).assign(PSF.viewPart(cs, cr, 0, slices - cs, rows - cr, cc));
        P1.viewPart(0, rows - cr, cols - cc, slices - cs, cr, cc).assign(PSF.viewPart(cs, 0, 0, slices - cs, cr, cc));

        P1.viewPart(slices - cs, 0, 0, cs, rows - cr, cols - cc).assign(PSF.viewPart(0, cr, cc, cs, rows - cr, cols - cc));
        P1.viewPart(slices - cs, 0, cols - cc, cs, rows - cr, cc).assign(PSF.viewPart(0, cr, 0, cs, rows - cr, cc));
        P1.viewPart(slices - cs, rows - cr, 0, cs, cr, cols - cc).assign(PSF.viewPart(0, 0, cc, cs, cr, cols - cc));
        P1.viewPart(slices - cs, rows - cr, cols - cc, cs, cr, cc).assign(PSF.viewPart(0, 0, 0, cs, cr, cc));
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
            new StackConverter(image).convertToGray8();
            break;
        case SHORT:
            new StackConverter(image).convertToGray16();
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
     * @param slicesPad
     *            number of slices in padded matrix
     * @param rowsPad
     *            number of rows in padded matrix
     * @param colsPad
     *            number of columns in padded matrix
     * 
     * @return padded matrix
     */
    public static FloatMatrix3D padPeriodic(final FloatMatrix3D X, final int slicesPad, final int rowsPad, final int colsPad) {
        final int slices = X.slices();
        final int rows = X.rows();
        final int cols = X.columns();
        if ((slices == slicesPad) && (rows == rowsPad) && (cols == colsPad)) {
            return X;
        }
        final FloatMatrix3D Xpad = new DenseFloatMatrix3D(slicesPad, rowsPad, colsPad);
        final int sOff = (slicesPad - slices + 1) / 2;
        final int rOff = (rowsPad - rows + 1) / 2;
        final int cOff = (colsPad - cols + 1) / 2;
        final float[] elemsXpad = (float[]) Xpad.elements();
        final int slicesStrideXpad = Xpad.rows() * Xpad.columns();
        final int rowsStrideXpad = Xpad.columns();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if (X.isView()) {
            if ((np > 1) && (X.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
                Future<?>[] futures = new Future[np];
                int k = slicesPad / np;
                for (int j = 0; j < np; j++) {
                    final int startslice = -sOff + j * k;
                    final int stopslice;
                    if (j == np - 1) {
                        stopslice = slicesPad - sOff;
                    } else {
                        stopslice = startslice + k;
                    }
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {

                        public void run() {
                            int cIn, rIn, sIn, cOut, rOut, sOut;
                            int idxXpad;
                            for (int s = startslice; s < stopslice; s++) {
                                sOut = s + sOff;
                                sIn = periodic(s, slices);
                                for (int r = -rOff; r < rowsPad - rOff; r++) {
                                    rOut = r + rOff;
                                    rIn = periodic(r, rows);
                                    for (int c = -cOff; c < colsPad - cOff; c++) {
                                        cOut = c + cOff;
                                        cIn = periodic(c, cols);
                                        idxXpad = sOut * slicesStrideXpad + rOut * rowsStrideXpad + cOut;
                                        elemsXpad[idxXpad] = X.getQuick(sIn, rIn, cIn);
                                    }
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int cIn, rIn, sIn, cOut, rOut, sOut;
                int idxXpad;
                for (int s = -sOff; s < slicesPad - sOff; s++) {
                    sOut = s + sOff;
                    sIn = periodic(s, slices);
                    for (int r = -rOff; r < rowsPad - rOff; r++) {
                        rOut = r + rOff;
                        rIn = periodic(r, rows);
                        for (int c = -cOff; c < colsPad - cOff; c++) {
                            cOut = c + cOff;
                            cIn = periodic(c, cols);
                            idxXpad = sOut * slicesStrideXpad + rOut * rowsStrideXpad + cOut;
                            elemsXpad[idxXpad] = X.getQuick(sIn, rIn, cIn);
                        }
                    }
                }
            }
        } else {
            final float[] elemsX = (float[]) X.elements();
            final int slicesStrideX = X.rows() * X.columns();
            final int rowsStrideX = X.columns();
            if ((np > 1) && (X.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
                Future<?>[] futures = new Future[np];
                int k = slicesPad / np;
                for (int j = 0; j < np; j++) {
                    final int startslice = -sOff + j * k;
                    final int stopslice;
                    if (j == np - 1) {
                        stopslice = slicesPad - sOff;
                    } else {
                        stopslice = startslice + k;
                    }
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {

                        public void run() {
                            int cIn, rIn, sIn, cOut, rOut, sOut;
                            int idxXpad;
                            int idxX;
                            for (int s = startslice; s < stopslice; s++) {
                                sOut = s + sOff;
                                sIn = periodic(s, slices);
                                for (int r = -rOff; r < rowsPad - rOff; r++) {
                                    rOut = r + rOff;
                                    rIn = periodic(r, rows);
                                    for (int c = -cOff; c < colsPad - cOff; c++) {
                                        cOut = c + cOff;
                                        cIn = periodic(c, cols);
                                        idxX = sIn * slicesStrideX + rIn * rowsStrideX + cIn;
                                        idxXpad = sOut * slicesStrideXpad + rOut * rowsStrideXpad + cOut;
                                        elemsXpad[idxXpad] = elemsX[idxX];
                                    }
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int cIn, rIn, sIn, cOut, rOut, sOut;
                int idxX;
                int idxXpad;
                for (int s = -sOff; s < slicesPad - sOff; s++) {
                    sOut = s + sOff;
                    sIn = periodic(s, slices);
                    for (int r = -rOff; r < rowsPad - rOff; r++) {
                        rOut = r + rOff;
                        rIn = periodic(r, rows);
                        for (int c = -cOff; c < colsPad - cOff; c++) {
                            cOut = c + cOff;
                            cIn = periodic(c, cols);
                            idxX = sIn * slicesStrideX + rIn * rowsStrideX + cIn;
                            idxXpad = sOut * slicesStrideXpad + rOut * rowsStrideXpad + cOut;
                            elemsXpad[idxXpad] = elemsX[idxX];
                        }
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
     * @param slicesPad
     *            number of slices in padded matrix
     * @param rowsPad
     *            number of rows in padded matrix
     * @param colsPad
     *            number of columns in padded matrix
     * 
     * @return padded matrix
     */
    public static FloatMatrix3D padReflexive(final FloatMatrix3D X, final int slicesPad, final int rowsPad, final int colsPad) {
        final int slices = X.slices();
        final int rows = X.rows();
        final int cols = X.columns();
        if ((slices == slicesPad) && (rows == rowsPad) && (cols == colsPad)) {
            return X;
        }
        final FloatMatrix3D Xpad = new DenseFloatMatrix3D(slicesPad, rowsPad, colsPad);
        final int sOff = (slicesPad - slices + 1) / 2;
        final int rOff = (rowsPad - rows + 1) / 2;
        final int cOff = (colsPad - cols + 1) / 2;
        final float[] elemsXpad = (float[]) Xpad.elements();
        final int slicesStrideXpad = Xpad.rows() * Xpad.columns();
        final int rowsStrideXpad = Xpad.columns();
        int np = ConcurrencyUtils.getNumberOfThreads();
        if (X.isView()) {
            if ((np > 1) && (X.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
                Future<?>[] futures = new Future[np];
                int k = slicesPad / np;
                for (int j = 0; j < np; j++) {
                    final int startslice = -sOff + j * k;
                    final int stopslice;
                    if (j == np - 1) {
                        stopslice = slicesPad - sOff;
                    } else {
                        stopslice = startslice + k;
                    }
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {

                        public void run() {
                            int cIn, rIn, sIn, cOut, rOut, sOut;
                            int idxXpad;
                            for (int s = startslice; s < stopslice; s++) {
                                sOut = s + sOff;
                                sIn = mirror(s, slices);
                                for (int r = -rOff; r < rowsPad - rOff; r++) {
                                    rOut = r + rOff;
                                    rIn = mirror(r, rows);
                                    for (int c = -cOff; c < colsPad - cOff; c++) {
                                        cOut = c + cOff;
                                        cIn = mirror(c, cols);
                                        idxXpad = sOut * slicesStrideXpad + rOut * rowsStrideXpad + cOut;
                                        elemsXpad[idxXpad] = X.getQuick(sIn, rIn, cIn);
                                    }
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int cIn, rIn, sIn, cOut, rOut, sOut;
                int idxXpad;
                for (int s = -sOff; s < slicesPad - sOff; s++) {
                    sOut = s + sOff;
                    sIn = mirror(s, slices);
                    for (int r = -rOff; r < rowsPad - rOff; r++) {
                        rOut = r + rOff;
                        rIn = mirror(r, rows);
                        for (int c = -cOff; c < colsPad - cOff; c++) {
                            cOut = c + cOff;
                            cIn = mirror(c, cols);
                            idxXpad = sOut * slicesStrideXpad + rOut * rowsStrideXpad + cOut;
                            elemsXpad[idxXpad] = X.getQuick(sIn, rIn, cIn);
                        }
                    }
                }
            }
        } else {
            final float[] elemsX = (float[]) X.elements();
            final int slicesStrideX = X.rows() * X.columns();
            final int rowsStrideX = X.columns();
            if ((np > 1) && (X.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
                Future<?>[] futures = new Future[np];
                int k = slicesPad / np;
                for (int j = 0; j < np; j++) {
                    final int startslice = -sOff + j * k;
                    final int stopslice;
                    if (j == np - 1) {
                        stopslice = slicesPad - sOff;
                    } else {
                        stopslice = startslice + k;
                    }
                    futures[j] = ConcurrencyUtils.submit(new Runnable() {

                        public void run() {
                            int cIn, rIn, sIn, cOut, rOut, sOut;
                            int idxX;
                            int idxXpad;
                            for (int s = startslice; s < stopslice; s++) {
                                sOut = s + sOff;
                                sIn = mirror(s, slices);
                                for (int r = -rOff; r < rowsPad - rOff; r++) {
                                    rOut = r + rOff;
                                    rIn = mirror(r, rows);
                                    for (int c = -cOff; c < colsPad - cOff; c++) {
                                        cOut = c + cOff;
                                        cIn = mirror(c, cols);
                                        idxX = sIn * slicesStrideX + rIn * rowsStrideX + cIn;
                                        idxXpad = sOut * slicesStrideXpad + rOut * rowsStrideXpad + cOut;
                                        elemsXpad[idxXpad] = elemsX[idxX];
                                    }
                                }
                            }
                        }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int cIn, rIn, sIn, cOut, rOut, sOut;
                int idxX;
                int idxXpad;
                for (int s = -sOff; s < slicesPad - sOff; s++) {
                    sOut = s + sOff;
                    sIn = mirror(s, slices);
                    for (int r = -rOff; r < rowsPad - rOff; r++) {
                        rOut = r + rOff;
                        rIn = mirror(r, rows);
                        for (int c = -cOff; c < colsPad - cOff; c++) {
                            cOut = c + cOff;
                            cIn = mirror(c, cols);
                            idxX = sIn * slicesStrideX + rIn * rowsStrideX + cIn;
                            idxXpad = sOut * slicesStrideXpad + rOut * rowsStrideXpad + cOut;
                            elemsXpad[idxXpad] = elemsX[idxX];
                        }
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
     * @param slicesPad
     *            number of slices in padded matrix
     * @param rowsPad
     *            number of rows in padded matrix
     * @param colsPad
     *            number of columns in padded matrix
     * 
     * @return padded matrix
     */
    public static FloatMatrix3D padZero(FloatMatrix3D X, int slicesPad, int rowsPad, int colsPad) {
        int slices = X.slices();
        int rows = X.rows();
        int cols = X.columns();
        if ((slices == slicesPad) && (rows == rowsPad) && (cols == colsPad)) {
            return X;
        }
        int kOff = (slicesPad - slices + 1) / 2;
        int jOff = (rowsPad - rows + 1) / 2;
        int iOff = (colsPad - cols + 1) / 2;
        FloatMatrix3D Xpad = new DenseFloatMatrix3D(slicesPad, rowsPad, colsPad);
        Xpad.viewPart(kOff, jOff, iOff, slices, rows, cols).assign(X);
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
    public static FloatMatrix3D padZero(FloatMatrix3D X, int[] padSize, PaddingType padding) {
        if ((padSize[0] == 0) && (padSize[1] == 0) & (padSize[2] == 0)) {
            return X;
        }
        FloatMatrix3D Xpad = null;
        switch (padding) {
        case BOTH:
            Xpad = new DenseFloatMatrix3D(X.slices() + 2 * padSize[0], X.rows() + 2 * padSize[1], X.columns() + 2 * padSize[2]);
            Xpad.viewPart(padSize[0], padSize[1], padSize[2], X.slices(), X.rows(), X.columns()).assign(X);
            break;
        case POST:
            Xpad = new DenseFloatMatrix3D(X.slices() + padSize[0], X.rows() + padSize[1], X.columns() + padSize[2]);
            Xpad.viewPart(0, 0, 0, X.slices(), X.rows(), X.columns()).assign(X);
            break;
        case PRE:
            Xpad = new DenseFloatMatrix3D(X.slices() + padSize[0], X.rows() + padSize[1], X.columns() + padSize[2]);
            Xpad.viewPart(padSize[0], padSize[1], padSize[2], X.slices(), X.rows(), X.columns()).assign(X);
            break;
        }
        return Xpad;
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
