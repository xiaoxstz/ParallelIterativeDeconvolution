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

package edu.emory.mathcs.restoretools.iterative.preconditioner;

import ij.IJ;

import java.util.concurrent.Future;

import cern.colt.function.tint.IntComparator;
import cern.colt.matrix.AbstractMatrix3D;
import cern.colt.matrix.tfcomplex.FComplexMatrix3D;
import cern.colt.matrix.tfcomplex.impl.DenseFComplexMatrix3D;
import cern.colt.matrix.tfloat.FloatMatrix1D;
import cern.colt.matrix.tfloat.FloatMatrix3D;
import cern.colt.matrix.tfloat.algo.FloatSorting;
import cern.colt.matrix.tfloat.impl.DenseFloatMatrix1D;
import cern.colt.matrix.tfloat.impl.DenseFloatMatrix3D;
import cern.jet.math.tfcomplex.FComplexFunctions;
import cern.jet.math.tfloat.FloatFunctions;
import edu.emory.mathcs.restoretools.iterative.FloatCommon3D;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PSFType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PaddingType;
import edu.emory.mathcs.restoretools.iterative.psf.FloatPSFMatrix3D;
import edu.emory.mathcs.utils.ConcurrencyUtils;

/**
 * 3D preconditioner based on the Fast Fourier Transform.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class FFTFloatPreconditioner3D implements FloatPreconditioner3D {
    private AbstractMatrix3D matdata;

    private float tol;

    private BoundaryType boundary;

    private int[] imSize;

    private int[] psfSize;

    private int[] padSize;

    /**
     * Creates a new instance of FloatFFTPreconditioner3D.
     * 
     * @param PSFMatrix
     *            PSF matrix
     * @param B
     *            blurred image
     * @param tol
     *            tolerance
     */
    public FFTFloatPreconditioner3D(FloatPSFMatrix3D PSFMatrix, FloatMatrix3D B, float tol) {
        this.tol = tol;
        this.boundary = PSFMatrix.getBoundary();
        this.imSize = new int[3];
        imSize[0] = B.slices();
        imSize[1] = B.rows();
        imSize[2] = B.columns();
        if (PSFMatrix.getType() == PSFType.INVARIANT) {
            this.psfSize = PSFMatrix.getInvPsfSize();
            this.padSize = PSFMatrix.getInvPadSize();
        } else {
            this.psfSize = PSFMatrix.getPSF().getSize();
            int[] minimal = new int[3];
            minimal[0] = psfSize[0] + imSize[0];
            minimal[1] = psfSize[1] + imSize[1];
            minimal[2] = psfSize[2] + imSize[2];
            switch (PSFMatrix.getResizing()) {
            case AUTO:
                int[] nextPowTwo = new int[3];
                if (!ConcurrencyUtils.isPowerOf2(minimal[0])) {
                    nextPowTwo[0] = ConcurrencyUtils.nextPow2(minimal[0]);
                } else {
                    nextPowTwo[0] = minimal[0];
                }
                if (!ConcurrencyUtils.isPowerOf2(minimal[1])) {
                    nextPowTwo[1] = ConcurrencyUtils.nextPow2(minimal[1]);
                } else {
                    nextPowTwo[1] = minimal[1];
                }
                if (!ConcurrencyUtils.isPowerOf2(minimal[2])) {
                    nextPowTwo[2] = ConcurrencyUtils.nextPow2(minimal[2]);
                } else {
                    nextPowTwo[2] = minimal[2];
                }
                if ((nextPowTwo[0] >= 1.5 * minimal[0]) || (nextPowTwo[1] >= 1.5 * minimal[1]) || (nextPowTwo[2] >= 1.5 * minimal[2])) {
                    //use minimal padding
                    psfSize[0] = minimal[0];
                    psfSize[1] = minimal[1];
                    psfSize[2] = minimal[2];
                } else {
                    psfSize[0] = nextPowTwo[0];
                    psfSize[1] = nextPowTwo[1];
                    psfSize[2] = nextPowTwo[2];
                }
                break;
            case MINIMAL:
                psfSize[0] = minimal[0];
                psfSize[1] = minimal[1];
                psfSize[2] = minimal[2];
                break;
            case NEXT_POWER_OF_TWO:
                psfSize[0] = minimal[0];
                psfSize[1] = minimal[1];
                psfSize[2] = minimal[2];
                if (!ConcurrencyUtils.isPowerOf2(psfSize[0])) {
                    psfSize[0] = ConcurrencyUtils.nextPow2(psfSize[0]);
                }
                if (!ConcurrencyUtils.isPowerOf2(psfSize[1])) {
                    psfSize[1] = ConcurrencyUtils.nextPow2(psfSize[1]);
                }
                if (!ConcurrencyUtils.isPowerOf2(psfSize[2])) {
                    psfSize[2] = ConcurrencyUtils.nextPow2(psfSize[2]);
                }
                break;
            }
            padSize = new int[3];
            if (imSize[0] < psfSize[0]) {
                padSize[0] = (psfSize[0] - imSize[0] + 1) / 2;
            }
            if (imSize[1] < psfSize[1]) {
                padSize[1] = (psfSize[1] - imSize[1] + 1) / 2;
            }
            if (imSize[2] < psfSize[2]) {
                padSize[2] = (psfSize[2] - imSize[2] + 1) / 2;
            }
        }
        constructMatrix(PSFMatrix.getPSF().getImage(), B, PSFMatrix.getPSF().getCenter());
    }

    public float getTolerance() {
        return tol;
    }

    public FloatMatrix1D solve(FloatMatrix1D b, boolean transpose) {
        FloatMatrix3D B = null;
        if (b.isView()) {
            B = new DenseFloatMatrix3D(imSize[0], imSize[1], imSize[2], (float[]) b.copy().elements(), 0, 0, 0, imSize[1] * imSize[2], imSize[2], 1, false);
        } else {
            B = new DenseFloatMatrix3D(imSize[0], imSize[1], imSize[2], (float[]) b.elements(), 0, 0, 0, imSize[1] * imSize[2], imSize[2], 1, false);
        }
        B = solve(B, transpose);
        return new DenseFloatMatrix1D(B.size(), (float[]) B.elements(), 0, 1, false);
    }

    public FloatMatrix3D solve(AbstractMatrix3D B, boolean transpose) {
        switch (boundary) {
        case ZERO:
            B = FloatCommon3D.padZero((FloatMatrix3D) B, psfSize[0], psfSize[1], psfSize[2]);
            break;
        case PERIODIC:
            B = FloatCommon3D.padPeriodic((FloatMatrix3D) B, psfSize[0], psfSize[1], psfSize[2]);
            break;
        case REFLEXIVE:
            B = FloatCommon3D.padReflexive((FloatMatrix3D) B, psfSize[0], psfSize[1], psfSize[2]);
            break;
        }
        B = ((DenseFloatMatrix3D) B).getFft3();
        if (transpose) {
            ((FComplexMatrix3D) B).assign((FComplexMatrix3D) matdata, FComplexFunctions.multConjSecond);
        } else {
            ((FComplexMatrix3D) B).assign((FComplexMatrix3D) matdata, FComplexFunctions.mult);
        }
        ((DenseFComplexMatrix3D) B).ifft3(true);
        return ((FComplexMatrix3D) B).viewPart(padSize[0], padSize[1], padSize[2], imSize[0], imSize[1], imSize[2]).getRealPart();
    }

    private void constructMatrix(FloatMatrix3D[][][] PSFs, FloatMatrix3D B, int[][][][] center) {
        matdata = PSFs[0][0][0].like();
        int[] center1 = center[0][0][0];
        int slices = PSFs.length;
        int rows = PSFs[0].length;
        int columns = PSFs[0][0].length;
        int size = slices * rows * columns;
        for (int s = 0; s < slices; s++) {
            for (int r = 0; r < rows; r++) {
                for (int c = 0; c < columns; c++) {
                    ((FloatMatrix3D) matdata).assign(PSFs[s][r][c], FloatFunctions.plus);
                }
            }
        }
        if (size != 1) {
            ((FloatMatrix3D) matdata).assign(FloatFunctions.div(size));
        }
        switch (boundary) {
        case ZERO:
            B = FloatCommon3D.padZero(B, psfSize[0], psfSize[1], psfSize[2]);
            break;
        case PERIODIC:
            B = FloatCommon3D.padPeriodic(B, psfSize[0], psfSize[1], psfSize[2]);
            break;
        case REFLEXIVE:
            B = FloatCommon3D.padReflexive(B, psfSize[0], psfSize[1], psfSize[2]);
            break;
        }
        precMatrixOnePsf(center1, B);
    }

    private void precMatrixOnePsf(int[] center, FloatMatrix3D Bpad) {
        int[] padSize = new int[3];
        padSize[0] = Bpad.slices() - matdata.slices();
        padSize[1] = Bpad.rows() - matdata.rows();
        padSize[2] = Bpad.columns() - matdata.columns();
        if ((padSize[0] > 0) || (padSize[1] > 0) || (padSize[2] > 0)) {
            matdata = FloatCommon3D.padZero((FloatMatrix3D) matdata, padSize, PaddingType.POST);
        }
        matdata = FloatCommon3D.circShift((FloatMatrix3D) matdata, center);
        matdata = ((DenseFloatMatrix3D) matdata).getFft3();
        AbstractMatrix3D E = ((FComplexMatrix3D) matdata).copy();
        ((FComplexMatrix3D) E).assign(FComplexFunctions.abs);
        E = ((FComplexMatrix3D) E).getRealPart();
        float[] maxAndLoc = ((FloatMatrix3D) E).getMaxLocation();
        final float maxE = maxAndLoc[0];

        if (tol == -1) {
            IJ.showStatus("Computing tolerance for preconditioner...");
            float[] minAndLoc = ((FloatMatrix3D) E).getMinLocation();
            float minE = minAndLoc[0];
            if (maxE / minE < 100) {
                tol = 0;
            } else {
                tol = defaultTol2(((FloatMatrix3D) E), Bpad);
            }
            IJ.showStatus("Computing tolerance for preconditioner...done.");
        }

        final float[] one = new float[] { 1, 0 };
        if (maxE != 1.0) {
            ((FComplexMatrix3D) matdata).assign(FComplexFunctions.div(new float[] { maxE, 0 }));
        }
        final int slices = E.slices();
        final int rows = E.rows();
        final int cols = E.columns();
        final float[] elementsE = (float[]) ((FloatMatrix3D) E).elements();
        final int zeroE = (int) ((FloatMatrix3D) E).index(0, 0, 0);
        final int sliceStrideE = ((FloatMatrix3D) E).sliceStride();
        final int rowStrideE = ((FloatMatrix3D) E).rowStride();
        final int columnStrideE = ((FloatMatrix3D) E).columnStride();
        final float[] elementsM = (float[]) ((FComplexMatrix3D) matdata).elements();
        final int zeroM = (int) ((FComplexMatrix3D) matdata).index(0, 0, 0);
        final int sliceStrideM = ((FComplexMatrix3D) matdata).sliceStride();
        final int rowStrideM = ((FComplexMatrix3D) matdata).rowStride();
        final int columnStrideM = ((FComplexMatrix3D) matdata).columnStride();

        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (slices * rows * cols >= ConcurrencyUtils.getThreadsBeginN_3D())) {
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
                        float[] elem = new float[2];
                        if (maxE != 1.0) {
                            for (int s = startslice; s < stopslice; s++) {
                                for (int r = 0; r < rows; r++) {
                                    for (int c = 0; c < cols; c++) {
                                        int idxE = zeroE + s * sliceStrideE + r * rowStrideE + c * columnStrideE;
                                        int idxM = zeroM + s * sliceStrideM + r * rowStrideM + c * columnStrideM;
                                        elem[0] = elementsM[idxM];
                                        elem[1] = elementsM[idxM + 1];
                                        if (elementsE[idxE] >= tol) {
                                            if (elem[1] != 0.0) {
                                                float scalar;
                                                if (Math.abs(elem[0]) >= Math.abs(elem[1])) {
                                                    scalar = (float) (1.0 / (elem[0] + elem[1] * (elem[1] / elem[0])));
                                                    elem[0] = scalar;
                                                    elem[1] = scalar * (-elem[1] / elem[0]);
                                                } else {
                                                    scalar = (float) (1.0 / (elem[0] * (elem[0] / elem[1]) + elem[1]));
                                                    elem[0] = scalar * (elem[0] / elem[1]);
                                                    elem[1] = -scalar;
                                                }
                                            } else {
                                                elem[0] = 1 / elem[0];
                                                elem[1] = 0;
                                            }
                                            elem[0] *= maxE;
                                            elem[1] *= maxE;
                                            elementsM[idxM] = elem[0];
                                            elementsM[idxM + 1] = elem[1];
                                        } else {
                                            elementsM[idxM] = one[0];
                                            elementsM[idxM + 1] = one[1];
                                        }
                                    }
                                }
                            }
                        } else {
                            for (int s = startslice; s < stopslice; s++) {
                                for (int r = 0; r < rows; r++) {
                                    for (int c = 0; c < cols; c++) {
                                        int idxE = zeroE + s * sliceStrideE + r * rowStrideE + c * columnStrideE;
                                        int idxM = zeroM + s * sliceStrideM + r * rowStrideM + c * columnStrideM;
                                        elem[0] = elementsM[idxM];
                                        elem[1] = elementsM[idxM + 1];

                                        if (elementsE[idxE] >= tol) {
                                            if (elem[1] != 0.0) {
                                                float scalar;
                                                if (Math.abs(elem[0]) >= Math.abs(elem[1])) {
                                                    scalar = (float) (1.0 / (elem[0] + elem[1] * (elem[1] / elem[0])));
                                                    elem[0] = scalar;
                                                    elem[1] = scalar * (-elem[1] / elem[0]);
                                                } else {
                                                    scalar = (float) (1.0 / (elem[0] * (elem[0] / elem[1]) + elem[1]));
                                                    elem[0] = scalar * (elem[0] / elem[1]);
                                                    elem[1] = -scalar;
                                                }
                                            } else {
                                                elem[0] = 1 / elem[0];
                                                elem[1] = 0;
                                            }
                                            elementsM[idxM] = elem[0];
                                            elementsM[idxM + 1] = elem[1];
                                        } else {
                                            elementsM[idxM] = one[0];
                                            elementsM[idxM + 1] = one[1];
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            float[] elem = new float[2];
            if (maxE != 1.0) {
                for (int s = 0; s < slices; s++) {
                    for (int r = 0; r < rows; r++) {
                        for (int c = 0; c < cols; c++) {
                            int idxE = zeroE + s * sliceStrideE + r * rowStrideE + c * columnStrideE;
                            int idxM = zeroM + s * sliceStrideM + r * rowStrideM + c * columnStrideM;
                            elem[0] = elementsM[idxM];
                            elem[1] = elementsM[idxM + 1];
                            if (elementsE[idxE] >= tol) {
                                if (elem[1] != 0.0) {
                                    float scalar;
                                    if (Math.abs(elem[0]) >= Math.abs(elem[1])) {
                                        scalar = (float) (1.0 / (elem[0] + elem[1] * (elem[1] / elem[0])));
                                        elem[0] = scalar;
                                        elem[1] = scalar * (-elem[1] / elem[0]);
                                    } else {
                                        scalar = (float) (1.0 / (elem[0] * (elem[0] / elem[1]) + elem[1]));
                                        elem[0] = scalar * (elem[0] / elem[1]);
                                        elem[1] = -scalar;
                                    }
                                } else {
                                    elem[0] = 1 / elem[0];
                                    elem[1] = 0;
                                }
                                elem[0] *= maxE;
                                elem[1] *= maxE;
                                elementsM[idxM] = elem[0];
                                elementsM[idxM + 1] = elem[1];
                            } else {
                                elementsM[idxM] = one[0];
                                elementsM[idxM + 1] = one[1];
                            }
                        }
                    }
                }
            } else {
                for (int s = 0; s < slices; s++) {
                    for (int r = 0; r < rows; r++) {
                        for (int c = 0; c < cols; c++) {
                            int idxE = zeroE + s * sliceStrideE + r * rowStrideE + c * columnStrideE;
                            int idxM = zeroM + s * sliceStrideM + r * rowStrideM + c * columnStrideM;
                            elem[0] = elementsM[idxM];
                            elem[1] = elementsM[idxM + 1];

                            if (elementsE[idxE] >= tol) {
                                if (elem[1] != 0.0) {
                                    float scalar;
                                    if (Math.abs(elem[0]) >= Math.abs(elem[1])) {
                                        scalar = (float) (1.0 / (elem[0] + elem[1] * (elem[1] / elem[0])));
                                        elem[0] = scalar;
                                        elem[1] = scalar * (-elem[1] / elem[0]);
                                    } else {
                                        scalar = (float) (1.0 / (elem[0] * (elem[0] / elem[1]) + elem[1]));
                                        elem[0] = scalar * (elem[0] / elem[1]);
                                        elem[1] = -scalar;
                                    }
                                } else {
                                    elem[0] = 1 / elem[0];
                                    elem[1] = 0;
                                }
                                elementsM[idxM] = elem[0];
                                elementsM[idxM + 1] = elem[1];
                            } else {
                                elementsM[idxM] = one[0];
                                elementsM[idxM + 1] = one[1];
                            }
                        }
                    }
                }
            }
        }
    }

    private float defaultTol2(FloatMatrix3D E, FloatMatrix3D B) {
        FloatMatrix1D s = new DenseFloatMatrix1D(E.size());
        System.arraycopy((float[]) E.elements(), 0, (float[]) s.elements(), 0, s.size());
        final float[] evalues = (float[]) s.elements();
        IntComparator compDec = new IntComparator() {
            public int compare(int a, int b) {
                if (evalues[a] != evalues[a] || evalues[b] != evalues[b])
                    return compareNaN(evalues[a], evalues[b]); // swap NaNs to
                // the end
                return evalues[a] < evalues[b] ? 1 : (evalues[a] == evalues[b] ? 0 : -1);
            }
        };
        int[] indices = FloatSorting.quickSort.sortIndex(s, compDec);
        s = s.viewSelection(indices);
        AbstractMatrix3D Bhat = ((DenseFloatMatrix3D) B).getFft3();
        ((FComplexMatrix3D) Bhat).assign(FComplexFunctions.abs);
        Bhat = ((FComplexMatrix3D) Bhat).getRealPart();
        FloatMatrix1D bhat = new DenseFloatMatrix1D(Bhat.size(), (float[]) ((FloatMatrix3D) Bhat).elements(), 0, 1, false);
        bhat = bhat.viewSelection(indices);
        bhat.assign(FloatFunctions.div((float) Math.sqrt(B.size())));
        int n = s.size();
        float[] rho = new float[n - 1];
        rho[n - 2] = bhat.getQuick(n - 1) * bhat.getQuick(n - 1);
        FloatMatrix1D G = new DenseFloatMatrix1D(n - 1);
        float[] elemsG = (float[]) G.elements();
        elemsG[n - 2] = rho[n - 2];
        float bhatel, temp1;
        for (int k = n - 2; k > 0; k--) {
            bhatel = bhat.getQuick(k);
            rho[k - 1] = rho[k] + bhatel * bhatel;
            temp1 = n - k;
            temp1 = temp1 * temp1;
            elemsG[k - 1] = rho[k - 1] / temp1;
        }
        for (int k = 0; k < n - 3; k++) {
            if (s.getQuick(k) == s.getQuick(k + 1)) {
                elemsG[k] = Float.POSITIVE_INFINITY;
            }
        }
        return s.getQuick((int) G.getMinLocation()[1]);
    }

    private final int compareNaN(float a, float b) {
        if (a != a) {
            if (b != b)
                return 0; // NaN equals NaN
            else
                return 1; // e.g. NaN > 5
        }
        return -1; // e.g. 5 < NaN
    }
}
