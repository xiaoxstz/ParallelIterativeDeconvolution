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
package edu.emory.mathcs.restoretools.iterative.mrnsd;

import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.AbstractDoubleIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.DoubleCommon2D;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PreconditionerType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.ResizingType;

/**
 * Modified Residual Norm Steepest Descent 2D. This is a nonnegatively
 * constrained steepest descent method.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class MRNSDDoubleIterativeDeconvolver2D extends AbstractDoubleIterativeDeconvolver2D {

    /**
     * If true, then the stopping tolerance is computed automatically.
     */
    protected boolean autoStoppingTol;

    /**
     * Stopping tolerance.
     */
    protected double stoppingTol;

    /**
     * Creates a new instance of MRNSDDoubleIterativeDeconvolver2D
     * 
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
     * @param maxIters
     *            maximal number of iterations
     * @param showIteration
     *            if true, then the restored image is displayed after each
     *            iteration
     * @param options
     *            MRNSD options
     */
    public MRNSDDoubleIterativeDeconvolver2D(ImagePlus imB, ImagePlus[][] imPSF, PreconditionerType preconditioner, double preconditionerTol, BoundaryType boundary, ResizingType resizing, OutputType output, int maxIters, boolean showIteration, MRNSDOptions options) {
        super("MRNSD", imB, imPSF, preconditioner, preconditionerTol, boundary, resizing, output, options.getUseThreshold(), options.getThreshold(), maxIters, showIteration, options.getLogConvergence());
        this.autoStoppingTol = options.getAutoStoppingTol();
        this.stoppingTol = options.getStoppingTol();
    }

    public ImagePlus deconvolve() {
        double alpha;
        double gamma;
        double theta;
        double rnrm;
        DoubleMatrix2D r, s, u, w;
        IntArrayList rowList, columnList;
        DoubleArrayList valueList;
        double tau = DoubleCommon2D.sqrteps;
        double sigsq = tau;
        double[] minAndLoc = B.getMinLocation();
        double minB = minAndLoc[0];
        if (minB < 0) {
            B.assign(DoubleFunctions.plus(-minB + sigsq));
        }

        if (autoStoppingTol) {
            stoppingTol = DoubleCommon2D.sqrteps * alg.vectorNorm2(A.times(B, true));
        }
        r = A.times(B, false);
        r.assign(B, DoubleFunctions.plusMultFirst(-1));
        if (P != null) {
            r = P.solve(r, false);
            r = P.solve(r, true);
            r = A.times(r, true);
            r.assign(DoubleFunctions.neg);
            gamma = B.aggregate(r, DoubleFunctions.plus, DoubleFunctions.multSquare);
            rnrm = alg.vectorNorm2(r);
        } else {
            r = A.times(r, true);
            r.assign(DoubleFunctions.neg);
            gamma = B.aggregate(r, DoubleFunctions.plus, DoubleFunctions.multSquare);
            rnrm = Math.sqrt(gamma);
        }
        ImagePlus imX = null;
        FloatProcessor ip = new FloatProcessor(bColumns, bRows);
        if (showIteration == true) {
            DoubleCommon2D.assignPixelsToProcessor(ip, B, cmY);
            imX = new ImagePlus("(deblurred)", ip);
            imX.show();
        }
        int k;
        rowList = new IntArrayList(B.size() / 2);
        columnList = new IntArrayList(B.size() / 2);
        valueList = new DoubleArrayList(B.size() / 2);
        for (k = 0; k < maxIters; k++) {
            if (rnrm <= stoppingTol) {
                IJ.log("MRNSD converged after " + k + "iterations.");
                break;
            }
            IJ.showStatus(name + " iteration: " + (k + 1) + "/" + maxIters);
            s = B.copy();
            s.assign(r, DoubleFunctions.multNeg);
            u = A.times(s, false);
            if (P != null) {
                u = P.solve(u, false);
            }
            theta = gamma / u.aggregate(DoubleFunctions.plus, DoubleFunctions.square);
            s.getNegativeValues(rowList, columnList, valueList);
            w = B.copy();
            w.assign(s, DoubleFunctions.divNeg, rowList, columnList);
            alpha = Math.min(theta, w.aggregate(DoubleFunctions.min, DoubleFunctions.identity, rowList, columnList));
            B.assign(s, DoubleFunctions.plusMultSecond(alpha));
            if (P != null) {
                w = P.solve(u, true);
                w = A.times(w, true);
                r.assign(w, DoubleFunctions.plusMultSecond(alpha));
                gamma = B.aggregate(r, DoubleFunctions.plus, DoubleFunctions.multSquare);
                rnrm = alg.vectorNorm2(r);
            } else {
                w = A.times(u, true);
                r.assign(w, DoubleFunctions.plusMultSecond(alpha));
                gamma = B.aggregate(r, DoubleFunctions.plus, DoubleFunctions.multSquare);
                rnrm = Math.sqrt(gamma);
            }
            if (logConvergence) {
                IJ.log((k + 1) + ".  Norm of the residual = " + rnrm);
            }
            if (showIteration == true) {
                if (useThreshold) {
                    DoubleCommon2D.assignPixelsToProcessor(ip, B, cmY, threshold);
                } else {
                    DoubleCommon2D.assignPixelsToProcessor(ip, B, cmY);
                }
                imX.updateAndDraw();
            }
        }
        if (logConvergence && k == maxIters) {
            IJ.log("MRNSD didn't converge. Reason: maximum number of iterations performed.");
        }
        if (showIteration == false) {
            if (useThreshold) {
                DoubleCommon2D.assignPixelsToProcessor(ip, B, cmY, threshold);
            } else {
                DoubleCommon2D.assignPixelsToProcessor(ip, B, cmY);
            }
            imX = new ImagePlus("(deblurred)", ip);
        }
        DoubleCommon2D.convertImage(imX, output);
        return imX;
    }
}
