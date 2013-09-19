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
package edu.emory.mathcs.restoretools.iterative.hybr;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import optimization.FloatFmin;
import optimization.FloatFmin_methods;
import cern.colt.list.tfloat.FloatArrayList;
import cern.colt.matrix.tdouble.algo.solver.HyBRInnerSolver;
import cern.colt.matrix.tdouble.algo.solver.HyBRRegularizationMethod;
import cern.colt.matrix.tfloat.FloatFactory2D;
import cern.colt.matrix.tfloat.FloatMatrix1D;
import cern.colt.matrix.tfloat.FloatMatrix2D;
import cern.colt.matrix.tfloat.algo.decomposition.DenseFloatSingularValueDecomposition;
import cern.colt.matrix.tfloat.impl.DenseColumnFloatMatrix2D;
import cern.colt.matrix.tfloat.impl.DenseFloatMatrix1D;
import cern.colt.matrix.tfloat.impl.DenseFloatMatrix2D;
import cern.jet.math.tfloat.FloatFunctions;
import cern.jet.stat.tfloat.FloatDescriptive;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.AbstractFloatIterativeDeconvolver3D;
import edu.emory.mathcs.restoretools.iterative.FloatCommon2D;
import edu.emory.mathcs.restoretools.iterative.FloatCommon3D;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PreconditionerType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.ResizingType;
import edu.emory.mathcs.restoretools.iterative.preconditioner.FloatPreconditioner3D;
import edu.emory.mathcs.restoretools.iterative.psf.FloatPSFMatrix3D;

/**
 * Hybrid Bidiagonalization Regularization 3D.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class HyBRFloatIterativeDeconvolver3D extends AbstractFloatIterativeDeconvolver3D {

    /**
     * Inner solver.
     */
    private HyBRInnerSolver innerSolver;

    /**
     * Regularization method.
     */
    private HyBRRegularizationMethod regMethod;

    /**
     * Regularization parameter.
     */
    private float regParam;

    /**
     * Omega parameter for weighted GCV.
     */
    private float omega;

    /**
     * If true then the reorthogonalization of Lanczos subspaces is performed.
     */
    private boolean reorth;

    /**
     * Begin regularization after this iteration.
     */
    private int begReg;

    /**
     * Tolerance for detecting flatness in the GCV curve as a % stopping
     * criteria.
     */
    private float flatTol;

    /**
     * The size of blurred image.
     */
    private int[] bsize;

    /**
     * Creates a new instance of HyBRFloatIterativeDeconvolver3D
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
     *            HyBR options
     */
    public HyBRFloatIterativeDeconvolver3D(ImagePlus imB, ImagePlus[][][] imPSF, PreconditionerType preconditioner, float preconditionerTol, BoundaryType boundary, ResizingType resizing, OutputType output, int maxIters, boolean showIteration, HyBROptions options) {
        super("HyBR", imB, imPSF, preconditioner, (float)preconditionerTol, boundary, resizing, output, options.getUseThreshold(), (float)options.getThreshold(), maxIters, showIteration, options.getLogConvergence());
        this.innerSolver = options.getInnerSolver();
        this.regMethod = options.getRegMethod();
        this.regParam = (float)options.getRegParameter();
        this.omega = (float)options.getOmega();
        this.reorth = options.getReorthogonalize();
        this.begReg = options.getBeginReg();
        this.flatTol = (float)options.getFlatTolerance();
        bsize = new int[] { bSlices, bRows, bColumns };
    }

    
    public ImagePlus deconvolve() {
        int k;
        int columns = A.getSize()[1];
        boolean bump = false;
        boolean warning = false;
        float rnrm = -1.0f;
        int iterationsSave = 0;
        float alpha, beta;
        HyBRInnerSolver inSolver = HyBRInnerSolver.NONE;
        FloatLBD lbd;
        FloatMatrix1D v;
        FloatMatrix1D work;
        FloatMatrix2D Ub, Vb;
        FloatMatrix1D f = null;
        FloatMatrix1D x = null;
        FloatMatrix1D xSave = null;
        float[] sv;
        FloatArrayList omegaList = new FloatArrayList(new float[begReg - 2]);
        FloatArrayList GCV = new FloatArrayList(new float[begReg - 2]);
        FloatMatrix1D b = new DenseFloatMatrix1D((int)B.size(), (float[]) B.elements(), 0, 1, false);
        FloatMatrix2D U = new DenseFloatMatrix2D(1, (int)B.size());
        FloatMatrix2D C = null;
        FloatMatrix2D V = null;
        DenseFloatSingularValueDecomposition svd;
        if (P == null) {
            beta = alg.norm2(b);
            U.viewRow(0).assign(b, FloatFunctions.multSecond(1.0f / beta));
            lbd = new FloatSimpleLBD(A, U, reorth);
        } else {
            work = P.solve(b, false);
            beta = alg.norm2(work);
            U.viewRow(0).assign(work, FloatFunctions.multSecond(1.0f / beta));
            lbd = new FloatPLBD(P, A, U, reorth);
        }
        ImagePlus imX = null;
        ImageStack is = new ImageStack(bColumns, bRows);
        if (showIteration) {
            FloatCommon3D.assignPixelsToStack(is, B, cmY);
            imX = new ImagePlus("(deblurred)", is);
            imX.show();
            IJ.showStatus("HyBR initialization...");
        }
        for (k = 0; k <= maxIters; k++) {
            lbd.apply();
            U = lbd.getU();
            C = lbd.getC();
            V = lbd.getV();
            v = new DenseFloatMatrix1D(C.columns() + 1);
            v.setQuick(0, beta);
            if (k >= 1) {
                IJ.showStatus("HyBR iteration: " + k + "/" + maxIters);
                if (k >= begReg - 1) {
                    inSolver = innerSolver;
                }
                switch (inSolver) {
                case TIKHONOV:
                    svd = alg.svd(C);
                    Ub = svd.getU();
                    sv = svd.getSingularValues();
                    Vb = svd.getV();
                    if (regMethod == HyBRRegularizationMethod.ADAPTWGCV) {
                        work = new DenseFloatMatrix1D(Ub.rows());
                        Ub.zMult(v, work, 1, 0, true);
                        omegaList.add(Math.min(1, findOmega(work, sv)));
                        omega = FloatDescriptive.mean(omegaList);
                    }
                    f = new DenseFloatMatrix1D(Vb.rows());
                    alpha = tikhonovSolver(Ub, sv, Vb, v, f);
                    GCV.add(GCVstopfun(alpha, Ub.viewRow(0), sv, beta, columns));
                    if (k > 1) {
                        if (Math.abs((GCV.getQuick(k - 1) - GCV.getQuick(k - 2))) / GCV.get(begReg - 2) < flatTol) {
                            x = V.zMult(f, null);
                            if (logConvergence) {
                                work = A.times(x, false);
                                rnrm = alg.norm2(work.assign(b, FloatFunctions.minus));
                                IJ.log(k + ".  Norm of the residual = " + rnrm);
                                IJ.log("HyBR stopped after " + k + " iterations.");
                                IJ.log("Reason for stopping: flat GCV curve.");
                            }
                            if (showIteration == false) {
                                if (useThreshold) {
                                    FloatCommon3D.assignPixelsToStack(is, x, bsize, cmY, threshold);
                                } else {
                                    FloatCommon3D.assignPixelsToStack(is, x, bsize, cmY);
                                }
                                imX = new ImagePlus("(deblurred)", is);
                            } else {
                                if (useThreshold) {
                                    FloatCommon3D.updatePixelsInStack(is, x, bsize, cmY, threshold);
                                } else {
                                    FloatCommon3D.updatePixelsInStack(is, x, bsize, cmY);
                                }
                                imX.updateAndDraw();
                            }
                            FloatCommon3D.convertImage(imX, output);
                            return imX;
                        } else if ((warning == true) && (GCV.size() > iterationsSave + 3)) {
                            for (int j = iterationsSave; j < GCV.size(); j++) {
                                if (GCV.getQuick(iterationsSave - 1) > GCV.get(j)) {
                                    bump = true;
                                }
                            }
                            if (bump == false) {
                                x.assign(xSave);
                                if (logConvergence) {
                                    work = A.times(x, false);
                                    rnrm = alg.norm2(work.assign(b, FloatFunctions.minus));
                                    IJ.log("HyBR stopped after " + iterationsSave + " iterations.");
                                    IJ.log("Reason for stopping: min of GCV curve within window of 4 iterations.");
                                }
                                if (showIteration == false) {
                                    if (useThreshold) {
                                        FloatCommon3D.assignPixelsToStack(is, x, bsize, cmY, threshold);
                                    } else {
                                        FloatCommon3D.assignPixelsToStack(is, x, bsize, cmY);
                                    }
                                    imX = new ImagePlus("(deblurred)", is);
                                } else {
                                    if (useThreshold) {
                                        FloatCommon3D.updatePixelsInStack(is, x, bsize, cmY, threshold);
                                    } else {
                                        FloatCommon3D.updatePixelsInStack(is, x, bsize, cmY);
                                    }
                                    imX.updateAndDraw();
                                }
                                FloatCommon3D.convertImage(imX, output);
                                return imX;
                            } else {
                                bump = false;
                                warning = false;
                                iterationsSave = maxIters;
                            }
                        } else if (warning == false) {
                            if (GCV.get(k - 2) < GCV.get(k - 1)) {
                                warning = true;
                                xSave = V.zMult(f, null);
                                iterationsSave = k;
                            }
                        }
                    }
                    break;
                case NONE:
                    f = alg.solve(C, v);
                    break;
                }
                x = V.zMult(f, null);
                if (logConvergence) {
                    work = A.times(x, false);
                    rnrm = alg.norm2(work.assign(b, FloatFunctions.minus));
                    IJ.log(k + ".  Norm of the residual = " + rnrm);
                }
                if (showIteration == true) {
                    if (useThreshold) {
                        FloatCommon3D.updatePixelsInStack(is, x, bsize, cmY, threshold);
                    } else {
                        FloatCommon3D.updatePixelsInStack(is, x, bsize, cmY);
                    }
                    imX.updateAndDraw();
                }
            }
        }
        if (logConvergence && k == (maxIters + 1)) {
            IJ.log("HyBR didn't converge. Reason: maximum number of iterations performed.");
        }
        if (showIteration == false) {
            if (useThreshold) {
                FloatCommon3D.assignPixelsToStack(is, x, bsize, cmY, threshold);
            } else {
                FloatCommon3D.assignPixelsToStack(is, x, bsize, cmY);
            }
            imX = new ImagePlus("(deblurred)", is);
        }
        FloatCommon3D.convertImage(imX, output);
        return imX;

    }

    private float findOmega(FloatMatrix1D bhat, float[] s) {
        int m = (int)bhat.size();
        int n = s.length;
        float alpha = s[n - 1];
        float t0 = bhat.viewPart(n, m - n).aggregate(FloatFunctions.plus, FloatFunctions.square);
        FloatMatrix1D s2 = new DenseFloatMatrix1D(s);
        s2.assign(FloatFunctions.square);
        float alpha2 = alpha * alpha;
        FloatMatrix1D tt = s2.copy();
        tt.assign(FloatFunctions.plus(alpha2));
        tt.assign(FloatFunctions.inv);
        float t1 = s2.aggregate(tt, FloatFunctions.plus, FloatFunctions.mult);
        s2 = new DenseFloatMatrix1D(s);
        s2.assign(FloatFunctions.mult(alpha));
        s2.assign(bhat.viewPart(0, n), FloatFunctions.mult);
        s2.assign(FloatFunctions.square);
        FloatMatrix1D work = tt.copy();
        work.assign(FloatFunctions.pow(3));
        work.assign(FloatFunctions.abs);
        float t3 = work.aggregate(s2, FloatFunctions.plus, FloatFunctions.mult);
        work = new DenseFloatMatrix1D(s);
        work.assign(tt, FloatFunctions.mult);
        float t4 = work.aggregate(FloatFunctions.plus, FloatFunctions.square);
        work = tt.copy();
        work.assign(bhat.viewPart(0, n), FloatFunctions.mult);
        work.assign(FloatFunctions.mult(alpha2));
        float t5 = work.aggregate(FloatFunctions.plus, FloatFunctions.square);
        s2 = new DenseFloatMatrix1D(s);
        s2.assign(bhat.viewPart(0, n), FloatFunctions.mult);
        s2.assign(FloatFunctions.square);
        tt.assign(FloatFunctions.pow(3));
        tt.assign(FloatFunctions.abs);
        float v2 = tt.aggregate(s2, FloatFunctions.plus, FloatFunctions.mult);
        return (m * alpha2 * v2) / (t1 * t3 + t4 * (t5 + t0));
    }

    private float tikhonovSolver(FloatMatrix2D U, float[] s, FloatMatrix2D V, FloatMatrix1D b, FloatMatrix1D x) {
        TikFmin fmin;
        FloatMatrix1D bhat = new DenseFloatMatrix1D(U.rows());
        U.zMult(b, bhat, 1, 0, true);
        float alpha = 0;
        switch (regMethod) {
        case GCV:
            fmin = new TikFmin(bhat, s, 1);
            alpha = FloatFmin.fmin(0, 1, fmin, FloatCommon2D.FMIN_TOL);
            break;
        case WGCV:
            fmin = new TikFmin(bhat, s, omega);
            alpha = FloatFmin.fmin(0, 1, fmin, FloatCommon2D.FMIN_TOL);
            break;
        case ADAPTWGCV:
            fmin = new TikFmin(bhat, s, omega);
            alpha = FloatFmin.fmin(0, 1, fmin, FloatCommon2D.FMIN_TOL);
            break;
        case NONE: // regularization parameter is given
            alpha = regParam;
            break;
        }
        FloatMatrix1D d = new DenseFloatMatrix1D(s);
        d.assign(FloatFunctions.square);
        d.assign(FloatFunctions.plus(alpha * alpha));
        bhat = bhat.viewPart(0, s.length);
        FloatMatrix1D S = new DenseFloatMatrix1D(s);
        bhat.assign(S, FloatFunctions.mult);
        bhat.assign(d, FloatFunctions.div);
        V.zMult(bhat, x);
        return alpha;
    }

    private static class TikFmin implements FloatFmin_methods {
        FloatMatrix1D bhat;

        float[] s;

        float omega;

        public TikFmin(FloatMatrix1D bhat, float[] s, float omega) {
            this.bhat = bhat;
            this.s = s;
            this.omega = omega;
        }

        public float f_to_minimize(float alpha) {
            int m = (int)bhat.size();
            int n = s.length;
            float t0 = bhat.viewPart(n, m - n).aggregate(FloatFunctions.plus, FloatFunctions.square);
            FloatMatrix1D s2 = new DenseFloatMatrix1D(s);
            s2.assign(FloatFunctions.square);
            float alpha2 = alpha * alpha;
            FloatMatrix1D work = s2.copy();
            work.assign(FloatFunctions.plus(alpha2));
            work.assign(FloatFunctions.inv);
            FloatMatrix1D t1 = work.copy();
            t1.assign(FloatFunctions.mult(alpha2));
            FloatMatrix1D t2 = t1.copy();
            t2.assign(bhat.viewPart(0, n), FloatFunctions.mult);
            FloatMatrix1D t3 = work.copy();
            t3.assign(s2, FloatFunctions.mult);
            t3.assign(FloatFunctions.mult(1 - omega));
            float denom = t3.aggregate(t1, FloatFunctions.plus, FloatFunctions.plus) + m - n;
            return n * (t2.aggregate(FloatFunctions.plus, FloatFunctions.square) + t0) / (denom * denom);
        }

    }

    private static float GCVstopfun(float alpha, FloatMatrix1D u, float[] s, float beta, int n) {
        int k = s.length;
        float beta2 = beta * beta;
        FloatMatrix1D s2 = new DenseFloatMatrix1D(s);
        s2.assign(FloatFunctions.square);
        float alpha2 = alpha * alpha;
        FloatMatrix1D t1 = s2.copy();
        t1.assign(FloatFunctions.plus(alpha2));
        t1.assign(FloatFunctions.inv);
        FloatMatrix1D t2 = t1.copy();
        t2.assign(u.viewPart(0, k), FloatFunctions.mult);
        t2.assign(FloatFunctions.mult(alpha2));
        float num = (float)(beta2 * (t2.aggregate(FloatFunctions.plus, FloatFunctions.square) + Math.pow(Math.abs(u.getQuick(k)), 2)) / (float) n);
        float den = (n - t1.aggregate(s2, FloatFunctions.plus, FloatFunctions.mult)) / (float) n;
        den = den * den;
        return num / den;
    }

    private interface FloatLBD {
        public void apply();

        public FloatMatrix2D getC();

        public FloatMatrix2D getU();

        public FloatMatrix2D getV();
    }

    private class FloatSimpleLBD implements FloatLBD {
        private final FloatFactory2D factory = FloatFactory2D.dense;

        private final FloatMatrix2D alphaBeta = new DenseFloatMatrix2D(2, 1);

        private final FloatPSFMatrix3D A;

        private FloatMatrix2D C;

        private FloatMatrix2D U;

        private FloatMatrix2D V;

        private boolean reorth;

        private int counter = 1;

        public FloatSimpleLBD(FloatPSFMatrix3D A, FloatMatrix2D U, boolean reorth) {
            this.A = A;
            this.reorth = reorth;
            this.U = U;
            this.V = null;
            this.C = null;
        }

        public void apply() {
            if (reorth) {
                int k = U.rows();
                FloatMatrix1D u = null;
                FloatMatrix1D v = null;
                FloatMatrix1D column = null;
                if (k == 1) {
                    v = A.times(U.viewRow(k - 1), true);
                } else {
                    v = A.times(U.viewRow(k - 1), true);
                    column = V.viewColumn(k - 2);
                    v.assign(column, FloatFunctions.plusMultSecond(-C.getQuick(k - 1, k - 2)));
                    for (int j = 0; j < k - 1; j++) {
                        column = V.viewColumn(j);
                        v.assign(column, FloatFunctions.plusMultSecond(-column.zDotProduct(v)));
                    }
                }
                float alpha = alg.norm2(v);
                v.assign(FloatFunctions.div(alpha));
                u = A.times(v, false);
                column = U.viewRow(k - 1);
                u.assign(column, FloatFunctions.plusMultSecond(-alpha));
                for (int j = 0; j < k; j++) {
                    column = U.viewRow(j);
                    u.assign(column, FloatFunctions.plusMultSecond(-column.zDotProduct(u)));
                }
                float beta = alg.norm2(u);
                alphaBeta.setQuick(0, 0, alpha);
                alphaBeta.setQuick(1, 0, beta);
                u.assign(FloatFunctions.div(beta));
                U = factory.appendRow(U, u);
                if (V == null) {
                    V = new DenseColumnFloatMatrix2D((int) v.size(), 1);
                    V.assign((float[]) v.elements());
                } else {
                    V = factory.appendColumn(V, v);
                }
                if (C == null) {
                    C = new DenseFloatMatrix2D(2, 1);
                    C.assign(alphaBeta);
                } else {
                    C = factory.composeBidiagonal(C, alphaBeta);
                }
            } else {
                FloatMatrix1D u = null;
                FloatMatrix1D v = null;
                FloatMatrix1D column = null;
                if (counter == 1) {
                    v = A.times(U.viewRow(0), true);
                } else {
                    v = A.times(U.viewRow(0), true);
                    column = V.viewColumn(counter - 2);
                    v.assign(column, FloatFunctions.plusMultSecond(-C.getQuick(counter - 1, counter - 2)));
                }
                float alpha = alg.norm2(v);
                v.assign(FloatFunctions.div(alpha));
                u = A.times(v, false);
                column = U.viewRow(0);
                u.assign(column, FloatFunctions.plusMultSecond(-alpha));
                float beta = alg.norm2(u);
                alphaBeta.setQuick(0, 0, alpha);
                alphaBeta.setQuick(1, 0, beta);
                u.assign(FloatFunctions.div(beta));
                U.viewRow(0).assign(u);
                if (V == null) {
                    V = new DenseColumnFloatMatrix2D((int) v.size(), 1);
                    V.assign((float[]) v.elements());
                } else {
                    V = factory.appendColumn(V, v);
                }
                if (C == null) {
                    C = new DenseFloatMatrix2D(2, 1);
                    C.assign(alphaBeta);
                } else {
                    C = factory.composeBidiagonal(C, alphaBeta);
                }
                counter++;
            }
        }

        public FloatMatrix2D getC() {
            return C;
        }

        public FloatMatrix2D getU() {
            return U;
        }

        public FloatMatrix2D getV() {
            return V;
        }
    }

    private class FloatPLBD implements FloatLBD {

        private final FloatFactory2D factory = FloatFactory2D.dense;

        private final FloatMatrix2D alphaBeta = new DenseFloatMatrix2D(2, 1);

        private final FloatPreconditioner3D P;

        private final FloatPSFMatrix3D A;

        private FloatMatrix2D C;

        private FloatMatrix2D U;

        private FloatMatrix2D V;

        private boolean reorth;
        
        private int counter = 1;

        public FloatPLBD(FloatPreconditioner3D P, FloatPSFMatrix3D A, FloatMatrix2D U, boolean reorth) {
            this.P = P;
            this.A = A;
            this.reorth = reorth;
            this.U = U;
            this.V = null;
            this.C = null;
        }

        public void apply() {
            if (reorth) {
                int k = U.rows();
                FloatMatrix1D u = null;
                FloatMatrix1D v = null;
                FloatMatrix1D row = null;
                if (k == 1) {
                    row = U.viewRow(k - 1).copy();
                    row = P.solve(row, true);
                    v = A.times(row, true);
                } else {
                    row = U.viewRow(k - 1).copy();
                    row = P.solve(row, true);
                    v = A.times(row, true);
                    row = V.viewColumn(k - 2);
                    v.assign(row, FloatFunctions.plusMultSecond(-C.getQuick(k - 1, k - 2)));
                    for (int j = 0; j < k - 1; j++) {
                        row = V.viewColumn(j);
                        v.assign(row, FloatFunctions.plusMultSecond(-row.zDotProduct(v)));
                    }
                }
                float alpha = alg.norm2(v);
                v.assign(FloatFunctions.div(alpha));
                row = A.times(v, false);
                u = P.solve(row, false);
                row = U.viewRow(k - 1);
                u.assign(row, FloatFunctions.plusMultSecond(-alpha));
                for (int j = 0; j < k; j++) {
                    row = U.viewRow(j);
                    u.assign(row, FloatFunctions.plusMultSecond(-row.zDotProduct(u)));
                }
                float beta = alg.norm2(u);
                alphaBeta.setQuick(0, 0, alpha);
                alphaBeta.setQuick(1, 0, beta);
                u.assign(FloatFunctions.div(beta));
                U = factory.appendRow(U, u);
                if (V == null) {
                    V = new DenseColumnFloatMatrix2D((int) v.size(), 1);
                    V.assign((float[]) v.elements());
                } else {
                    V = factory.appendColumn(V, v);
                }
                if (C == null) {
                    C = new DenseFloatMatrix2D(2, 1);
                    C.assign(alphaBeta);
                } else {
                    C = factory.composeBidiagonal(C, alphaBeta);
                }
            } else {
                FloatMatrix1D u = null;
                FloatMatrix1D v = null;
                FloatMatrix1D row = null;
                if (counter == 1) {
                    row = U.viewRow(0).copy();
                    row = P.solve(row, true);
                    v = A.times(row, true);
                } else {
                    row = U.viewRow(0).copy();
                    row = P.solve(row, true);
                    v = A.times(row, true);
                    row = V.viewColumn(counter - 2);
                    v.assign(row, FloatFunctions.plusMultSecond(-C.getQuick(counter - 1, counter - 2)));
                }
                float alpha = alg.norm2(v);
                v.assign(FloatFunctions.div(alpha));
                row = A.times(v, false);
                u = P.solve(row, false);
                row = U.viewRow(0);
                u.assign(row, FloatFunctions.plusMultSecond(-alpha));
                float beta = alg.norm2(u);
                alphaBeta.setQuick(0, 0, alpha);
                alphaBeta.setQuick(1, 0, beta);
                u.assign(FloatFunctions.div(beta));
                U.viewRow(0).assign(u);
                if (V == null) {
                    V = new DenseColumnFloatMatrix2D((int) v.size(), 1);
                    V.assign((float[]) v.elements());
                } else {
                    V = factory.appendColumn(V, v);
                }
                if (C == null) {
                    C = new DenseFloatMatrix2D(2, 1);
                    C.assign(alphaBeta);
                } else {
                    C = factory.composeBidiagonal(C, alphaBeta);
                }
                counter++;
            }
        }
        public FloatMatrix2D getC() {
            return C;
        }

        public FloatMatrix2D getU() {
            return U;
        }

        public FloatMatrix2D getV() {
            return V;
        }
    }
}
