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

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;
import cern.colt.matrix.tdouble.DoubleMatrix3D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix3D;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PreconditionerType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.ResizingType;
import edu.emory.mathcs.restoretools.iterative.preconditioner.DoublePreconditioner3D;
import edu.emory.mathcs.restoretools.iterative.preconditioner.FFTDoublePreconditioner3D;
import edu.emory.mathcs.restoretools.iterative.psf.DoublePSFMatrix3D;

/**
 * Abstract iterative deconvolver 3D.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public abstract class AbstractDoubleIterativeDeconvolver3D implements IterativeDeconvolver {

    /**
     * Algebra.
     */
    protected static final DenseDoubleAlgebra alg = DenseDoubleAlgebra.DEFAULT;

    /**
     * Blurred image.
     */
    protected DoubleMatrix3D B;

    /**
     * Point Spread Function.
     */
    protected DoublePSFMatrix3D A;

    /**
     * Preconditioner.
     */
    protected DoublePreconditioner3D P;

    /**
     * Color model.
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
    protected double threshold;

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
     * The name of a deconvolution algorithm.
     */
    protected String name;

    /**
     * Creates a new instance of AbstractDoubleIterativeDeconvolver3D
     * 
     * @param name
     *            name of a deconvolution algorithm
     * @param imB
     *            blurred image
     * @param imPSF
     *            Point Spread Function
     * @param preconditioner
     *            type of a preconditioner
     * @param preconditionerTol
     *            tolerance for the preconditioner
     * @param boundary
     *            type of boundary conditions
     * @param resizing
     *            type of resizing
     * @param output
     *            type of the output image
     * @param useThreshold
     *            if true, then the thresholding is performed
     * @param threshold
     *            the smallest nonnegative value assigned to the restored image
     * @param maxIters
     *            maximal number of iterations
     * @param showIteration
     *            if true, then the restored image is displayed after each
     *            iteration
     * @param logConvergence
     *            if true, then the convergence information is displayed after
     *            each iteration
     */
    protected AbstractDoubleIterativeDeconvolver3D(String name, ImagePlus imB, ImagePlus[][][] imPSF, PreconditionerType preconditioner, double preconditionerTol, BoundaryType boundary, ResizingType resizing, OutputType output, boolean useThreshold, double threshold, int maxIters,
            boolean showIteration, boolean logConvergence) {
        this.name = name;
        IJ.showStatus(name + ": initialization...");
        ImageStack isB = imB.getStack();
        ImageProcessor ipB = imB.getProcessor();
        this.cmY = ipB.getColorModel();
        this.bSlices = isB.getSize();
        this.bColumns = ipB.getWidth();
        this.bRows = ipB.getHeight();
        B = new DenseDoubleMatrix3D(this.bSlices, this.bRows, this.bColumns);
        A = new DoublePSFMatrix3D(imPSF, boundary, resizing, new int[] { this.bSlices, this.bRows, this.bColumns });
        DoubleCommon3D.assignPixelsToMatrix(isB, B);
        switch (preconditioner) {
        case FFT:
            this.P = new FFTDoublePreconditioner3D(A, B, preconditionerTol);
            break;
        case NONE:
            this.P = null;
            break;
        default:
            throw new IllegalArgumentException("Unsupported preconditioner type.");
        }
        IJ.showStatus(name + ": initialization...");
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
        this.useThreshold = useThreshold;
        this.threshold = threshold;
        this.maxIters = maxIters;
        this.showIteration = showIteration;
        this.logConvergence = logConvergence;

    }

    /**
     * Returns 3D preconditioner.
     * 
     * @return preconditioner
     */
    public DoublePreconditioner3D getPreconditioner() {
        return P;
    }
}
