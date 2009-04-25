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
import cern.colt.list.tfloat.FloatArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tfloat.FloatMatrix2D;
import cern.jet.math.tfloat.FloatFunctions;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.AbstractFloatIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.FloatCommon2D;
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
public class MRNSDFloatIterativeDeconvolver2D extends AbstractFloatIterativeDeconvolver2D {

    /**
     * If true, then the stopping tolerance is computed automatically.
     */
    protected boolean autoStoppingTol;

    /**
     * Stopping tolerance.
     */
    protected float stoppingTol;

    /**
     * Creates a new instance of MRNSDFloatIterativeDeconvolver2D
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
    public MRNSDFloatIterativeDeconvolver2D(ImagePlus imB, ImagePlus[][] imPSF, PreconditionerType preconditioner, float preconditionerTol, BoundaryType boundary, ResizingType resizing, OutputType output, int maxIters, boolean showIteration, MRNSDOptions options) {
        super("MRNSD", imB, imPSF, preconditioner, preconditionerTol, boundary, resizing, output, options.getUseThreshold(), (float) options.getThreshold(), maxIters, showIteration, options.getLogConvergence());
        this.autoStoppingTol = options.getAutoStoppingTol();
        this.stoppingTol = (float) options.getStoppingTol();
    }

    public ImagePlus deconvolve() {
        float alpha;
        float gamma;
        float theta;
        float rnrm;
        FloatMatrix2D r, s, u, w;
        IntArrayList rowList, columnList;
        FloatArrayList valueList;
        float tau = FloatCommon2D.sqrteps;
        float sigsq = tau;
        float[] minAndLoc = B.getMinLocation();
        float minB = minAndLoc[0];
        if (minB < 0) {
            B.assign(FloatFunctions.plus(-minB + sigsq));
        }
        if (autoStoppingTol) {
            stoppingTol = FloatCommon2D.sqrteps * alg.vectorNorm2(A.times(B, true));
        }
        r = A.times(B, false);
        r.assign(B, FloatFunctions.plusMultFirst(-1));
        if (P != null) {
            r = P.solve(r, false);
            r = P.solve(r, true);
            r = A.times(r, true);
            r.assign(FloatFunctions.neg);
            gamma = B.aggregate(r, FloatFunctions.plus, FloatFunctions.multSquare);
            rnrm = alg.vectorNorm2(r);
        } else {
            r = A.times(r, true);
            r.assign(FloatFunctions.neg);
            gamma = B.aggregate(r, FloatFunctions.plus, FloatFunctions.multSquare);
            rnrm = (float) Math.sqrt(gamma);
        }
        ImagePlus imX = null;
        FloatProcessor ip = new FloatProcessor(bColumns, bRows);
        if (showIteration == true) {
            FloatCommon2D.assignPixelsToProcessor(ip, B, cmY);
            imX = new ImagePlus("(deblurred)", ip);
            imX.show();
        }
        int k;
        rowList = new IntArrayList(B.size() / 2);
        columnList = new IntArrayList(B.size() / 2);
        valueList = new FloatArrayList(B.size() / 2);
        for (k = 0; k < maxIters; k++) {
            if (rnrm <= stoppingTol) {
                IJ.log("MRNSD converged after " + k + "iterations.");
                break;
            }
            IJ.showStatus(name + " iteration: " + (k + 1) + "/" + maxIters);
            s = B.copy();
            s.assign(r, FloatFunctions.multNeg);
            u = A.times(s, false);
            if (P != null) {
                u = P.solve(u, false);
            }
            theta = gamma / u.aggregate(FloatFunctions.plus, FloatFunctions.square);
            s.getNegativeValues(rowList, columnList, valueList);
            w = B.copy();
            w.assign(s, FloatFunctions.divNeg, rowList, columnList);
            alpha = Math.min(theta, w.aggregate(FloatFunctions.min, FloatFunctions.identity, rowList, columnList));
            B.assign(s, FloatFunctions.plusMultSecond(alpha));
            if (P != null) {
                w = P.solve(u, true);
                w = A.times(w, true);
                r.assign(w, FloatFunctions.plusMultSecond(alpha));
                gamma = B.aggregate(r, FloatFunctions.plus, FloatFunctions.multSquare);
                rnrm = alg.vectorNorm2(r);
            } else {
                w = A.times(u, true);
                r.assign(w, FloatFunctions.plusMultSecond(alpha));
                gamma = B.aggregate(r, FloatFunctions.plus, FloatFunctions.multSquare);
                rnrm = (float) Math.sqrt(gamma);
            }
            if (logConvergence) {
                IJ.log((k + 1) + ".  Norm of the residual = " + rnrm);
            }
            if (showIteration == true) {
                if (useThreshold) {
                    FloatCommon2D.assignPixelsToProcessor(ip, B, cmY, threshold);
                } else {
                    FloatCommon2D.assignPixelsToProcessor(ip, B, cmY);
                }
                imX.updateAndDraw();
            }
        }
        if (logConvergence && k == maxIters) {
            IJ.log("MRNSD didn't converge. Reason: maximum number of iterations performed.");
        }
        if (showIteration == false) {
            if (useThreshold) {
                FloatCommon2D.assignPixelsToProcessor(ip, B, cmY, threshold);
            } else {
                FloatCommon2D.assignPixelsToProcessor(ip, B, cmY);
            }
            imX = new ImagePlus("(deblurred)", ip);
        }
        FloatCommon2D.convertImage(imX, output);
        return imX;
    }
}
