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
package edu.emory.mathcs.restoretools.iterative.cgls;

import ij.IJ;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import cern.colt.matrix.tfloat.FloatMatrix2D;
import cern.jet.math.tfloat.FloatFunctions;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.AbstractFloatIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.FloatCommon2D;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PreconditionerType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.ResizingType;

/**
 * Conjugate Gradient for Least Squares 2D.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class CGLSFloatIterativeDeconvolver2D extends AbstractFloatIterativeDeconvolver2D {

    /**
     * If true, then the stopping tolerance is computed automatically.
     */
    protected boolean autoStoppingTol;

    /**
     * Stopping tolerance.
     */
    protected float stoppingTol;

    /**
     * Creates a new instance of CGLSFloatIterativeDeconvolver2D
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
     *            CGLS options
     */
    public CGLSFloatIterativeDeconvolver2D(ImagePlus imB, ImagePlus[][] imPSF, PreconditionerType preconditioner, float preconditionerTol, BoundaryType boundary, ResizingType resizing, OutputType output, int maxIters, boolean showIteration, CGLSOptions options) {
        super("CGLS", imB, imPSF, preconditioner, preconditionerTol, boundary, resizing, output, options.getUseThreshold(), (float) options.getThreshold(), maxIters, showIteration, options.getLogConvergence());
        this.autoStoppingTol = options.getAutoStoppingTol();
        this.stoppingTol = (float) options.getStoppingTol();
    }

    public ImagePlus deconvolve() {
        FloatMatrix2D p, q, r, s;
        float alpha;
        float beta;
        float gamma;
        float oldgamma = 0;
        float nq;
        float rnrm;

        if (autoStoppingTol) {
            stoppingTol = FloatCommon2D.sqrteps * alg.vectorNorm2(A.times(B, true));
        }
        s = A.times(B, false);
        s.assign(B, FloatFunctions.plusMultFirst(-1));
        r = A.times(s, true);
        rnrm = alg.vectorNorm2(r);
        if (P != null) {
            r = P.solve(r, true);
            gamma = alg.vectorNorm2(r);
            gamma *= gamma;
        } else {
            gamma = rnrm;
            gamma *= gamma;
        }
        ImagePlus imX = null;
        FloatProcessor ip = new FloatProcessor(bColumns, bRows);
        if (showIteration == true) {
            FloatCommon2D.assignPixelsToProcessor(ip, B, cmY);
            imX = new ImagePlus("(deblurred)", ip);
            imX.show();
        }
        p = r.copy();
        int k;
        for (k = 0; k < maxIters; k++) {
            if (rnrm <= stoppingTol) {
                IJ.log("CGLS converged after " + k + "iterations.");
                break;
            }
            IJ.showStatus(name + "iteration: " + (k + 1) + "/" + maxIters);
            if (k >= 1) {
                beta = gamma / oldgamma;
                p.assign(r, FloatFunctions.plusMultFirst(beta));
            }
            if (P != null) {
                r = P.solve(p, false);
                q = A.times(r, false);
            } else {
                q = A.times(p, false);
            }
            nq = alg.vectorNorm2(q);
            nq = nq * nq;
            alpha = gamma / nq;
            if (P != null) {
                B.assign(r, FloatFunctions.plusMultSecond(alpha));
            } else {
                B.assign(p, FloatFunctions.plusMultSecond(alpha));
            }
            s.assign(q, FloatFunctions.plusMultSecond(-alpha));
            r = A.times(s, true);
            rnrm = alg.vectorNorm2(r);
            if (P != null) {
                r = P.solve(r, true);
                oldgamma = gamma;
                gamma = alg.vectorNorm2(r);
                gamma *= gamma;
            } else {
                oldgamma = gamma;
                gamma = rnrm;
                gamma *= gamma;
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
            IJ.log("CGLS didn't converge. Reason: maximum number of iterations performed.");
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
