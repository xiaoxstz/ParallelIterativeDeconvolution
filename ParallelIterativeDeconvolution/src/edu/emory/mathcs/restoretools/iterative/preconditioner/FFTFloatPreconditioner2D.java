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

import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import cern.colt.function.tint.IntComparator;
import cern.colt.matrix.AbstractMatrix2D;
import cern.colt.matrix.tfcomplex.FComplexMatrix2D;
import cern.colt.matrix.tfcomplex.impl.DenseFComplexMatrix2D;
import cern.colt.matrix.tfloat.FloatMatrix1D;
import cern.colt.matrix.tfloat.FloatMatrix2D;
import cern.colt.matrix.tfloat.algo.FloatSorting;
import cern.colt.matrix.tfloat.impl.DenseFloatMatrix1D;
import cern.colt.matrix.tfloat.impl.DenseFloatMatrix2D;
import cern.jet.math.tfcomplex.FComplex;
import cern.jet.math.tfcomplex.FComplexFunctions;
import cern.jet.math.tfloat.FloatFunctions;
import edu.emory.mathcs.restoretools.iterative.FloatCommon2D;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PSFType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PaddingType;
import edu.emory.mathcs.restoretools.iterative.psf.FloatPSFMatrix2D;
import edu.emory.mathcs.utils.ConcurrencyUtils;

/**
 * 2D preconditioner based on the Fast Fourier Transform.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class FFTFloatPreconditioner2D implements FloatPreconditioner2D {
    private AbstractMatrix2D matdata;

    private float tol;

    private BoundaryType boundary;

    private int[] imSize;

    private int[] psfSize;

    private int[] padSize;

    /**
     * Creates a new instance of FloatFFTPreconditioner2D.
     * 
     * @param PSFMatrix
     *            PSF matrix
     * @param B
     *            blurred image
     * @param tol
     *            tolerance
     */
    public FFTFloatPreconditioner2D(FloatPSFMatrix2D PSFMatrix, FloatMatrix2D B, float tol) {
        this.tol = tol;
        this.boundary = PSFMatrix.getBoundary();
        this.imSize = new int[2];
        imSize[0] = B.rows();
        imSize[1] = B.columns();
        if (PSFMatrix.getType() == PSFType.INVARIANT) {
            this.psfSize = PSFMatrix.getInvPsfSize();
            this.padSize = PSFMatrix.getInvPadSize();
        } else {
            this.psfSize = PSFMatrix.getPSF().getSize();
            int[] minimal = new int[2];
            minimal[0] = psfSize[0] + imSize[0];
            minimal[1] = psfSize[1] + imSize[1];
            switch (PSFMatrix.getResizing()) {
            case AUTO:
                int[] nextPowTwo = new int[2];
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
                if ((nextPowTwo[0] >= 1.5 * minimal[0]) || (nextPowTwo[1] >= 1.5 * minimal[1])) {
                    //use minimal padding
                    psfSize[0] = minimal[0];
                    psfSize[1] = minimal[1];
                } else {
                    psfSize[0] = nextPowTwo[0];
                    psfSize[1] = nextPowTwo[1];
                }
                break;
            case MINIMAL:
                psfSize[0] = minimal[0];
                psfSize[1] = minimal[1];
                break;
            case NEXT_POWER_OF_TWO:
                psfSize[0] = minimal[0];
                psfSize[1] = minimal[1];
                if (!ConcurrencyUtils.isPowerOf2(psfSize[0])) {
                    psfSize[0] = ConcurrencyUtils.nextPow2(psfSize[0]);
                }
                if (!ConcurrencyUtils.isPowerOf2(psfSize[1])) {
                    psfSize[1] = ConcurrencyUtils.nextPow2(psfSize[1]);
                }
                break;
            }
            padSize = new int[2];
            if (imSize[0] < psfSize[0]) {
                padSize[0] = (psfSize[0] - imSize[0] + 1) / 2;
            }
            if (imSize[1] < psfSize[1]) {
                padSize[1] = (psfSize[1] - imSize[1] + 1) / 2;
            }

        }
        constructMatrix(PSFMatrix.getPSF().getImage(), B, PSFMatrix.getPSF().getCenter());
    }

    public float getTolerance() {
        return tol;
    }

    public FloatMatrix1D solve(FloatMatrix1D b, boolean transpose) {
        FloatMatrix2D B = null;
        if (b.isView()) {
            B = new DenseFloatMatrix2D(imSize[0], imSize[1], (float[]) b.copy().elements(), 0, 0, imSize[1], 1, false);
        } else {
            B = new DenseFloatMatrix2D(imSize[0], imSize[1], (float[]) b.elements(), 0, 0, imSize[1], 1, false);
        }
        B = solve(B, transpose);
        return new DenseFloatMatrix1D((int)B.size(), (float[]) B.elements(), 0, 1, false);
    }

    public FloatMatrix2D solve(AbstractMatrix2D B, boolean transpose) {
        switch (boundary) {
        case ZERO:
            B = FloatCommon2D.padZero((FloatMatrix2D) B, psfSize[0], psfSize[1]);
            break;
        case PERIODIC:
            B = FloatCommon2D.padPeriodic((FloatMatrix2D) B, psfSize[0], psfSize[1]);
            break;
        case REFLEXIVE:
            B = FloatCommon2D.padReflexive((FloatMatrix2D) B, psfSize[0], psfSize[1]);
            break;
        }
        B = ((DenseFloatMatrix2D) B).getFft2();
        if (transpose) {
            ((FComplexMatrix2D) B).assign((FComplexMatrix2D) matdata, FComplexFunctions.multConjSecond);
        } else {
            ((FComplexMatrix2D) B).assign((FComplexMatrix2D) matdata, FComplexFunctions.mult);
        }
        ((DenseFComplexMatrix2D) B).ifft2(true);
        return ((FComplexMatrix2D) B).viewPart(padSize[0], padSize[1], imSize[0], imSize[1]).getRealPart();
    }

    private void constructMatrix(FloatMatrix2D[][] PSFs, FloatMatrix2D B, int[][][] center) {
        matdata = PSFs[0][0].like();
        int[] center1 = center[0][0];
        int rows = PSFs.length;
        int columns = PSFs[0].length;
        int size = rows * columns;
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < columns; c++) {
                ((FloatMatrix2D) matdata).assign(PSFs[r][c], FloatFunctions.plus);
            }
        }
        if (size != 1) {
            ((FloatMatrix2D) matdata).assign(FloatFunctions.div(size));
        }
        switch (boundary) {
        case ZERO:
            B = FloatCommon2D.padZero(B, psfSize[0], psfSize[1]);
            break;
        case PERIODIC:
            B = FloatCommon2D.padPeriodic(B, psfSize[0], psfSize[1]);
            break;
        case REFLEXIVE:
            B = FloatCommon2D.padReflexive(B, psfSize[0], psfSize[1]);
            break;
        }
        precMatrixOnePsf(center1, B);
    }

    private void precMatrixOnePsf(int[] center, FloatMatrix2D Bpad) {
        int[] padSize = new int[2];
        padSize[0] = Bpad.rows() - matdata.rows();
        padSize[1] = Bpad.columns() - matdata.columns();
        if ((padSize[0] > 0) || (padSize[1] > 0)) {
            matdata = FloatCommon2D.padZero((FloatMatrix2D) matdata, padSize, PaddingType.POST);
        }
        matdata = FloatCommon2D.circShift((FloatMatrix2D) matdata, center);
        matdata = ((DenseFloatMatrix2D) matdata).getFft2();
        AbstractMatrix2D E = ((FComplexMatrix2D) matdata).copy();
        ((FComplexMatrix2D) E).assign(FComplexFunctions.abs);
        E = ((FComplexMatrix2D) E).getRealPart();
        float[] maxAndLoc = ((FloatMatrix2D) E).getMaxLocation();
        final float maxE = maxAndLoc[0];

        if (tol == -1) {
            IJ.showStatus("Computing tolerance for preconditioner...");
            float[] minAndLoc = ((FloatMatrix2D) E).getMinLocation();
            float minE = minAndLoc[0];
            if (maxE / minE < 100) {
                tol = 0;
            } else {
                tol = defaultTol2(((FloatMatrix2D) E), Bpad);
            }
            IJ.showStatus("Computing tolerance for preconditioner...done.");
        }

        final float[] one = new float[] { 1, 0 };
        if (maxE != 1.0) {
            ((FComplexMatrix2D) matdata).assign(FComplexFunctions.div(new float[] { maxE, 0 }));
        }
        final int rows = E.rows();
        final int cols = E.columns();
        final float[] elementsE = (float[]) ((FloatMatrix2D) E).elements();
        final int zeroE = (int) ((FloatMatrix2D) E).index(0, 0);
        final int rowStrideE = ((FloatMatrix2D) E).rowStride();
        final int columnStrideE = ((FloatMatrix2D) E).columnStride();
        final float[] elementsM = (float[]) ((FComplexMatrix2D) matdata).elements();
        final int zeroM = (int) ((FComplexMatrix2D) matdata).index(0, 0);
        final int rowStrideM = ((FComplexMatrix2D) matdata).rowStride();
        final int columnStrideM = ((FComplexMatrix2D) matdata).columnStride();

        int np = ConcurrencyUtils.getNumberOfThreads();
        if ((np > 1) && (rows * cols >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            Future<?>[] futures = new Future[np];
            int k = rows / np;
            for (int j = 0; j < np; j++) {
                final int startrow = j * k;
                final int stoprow;
                if (j == np - 1) {
                    stoprow = rows;
                } else {
                    stoprow = startrow + k;
                }
                futures[j] = ConcurrencyUtils.submit(new Runnable() {
                    public void run() {
                        float[] elem = new float[2];
                        int idxE = zeroE + startrow * rowStrideE;
                        int idxM = zeroM + startrow * rowStrideM;
                        if (maxE != 1.0) {
                            for (int r = startrow; r < stoprow; r++) {
                                for (int iE = idxE, iM = idxM, c = 0; c < cols; c++) {
                                    elem[0] = elementsM[iM];
                                    elem[1] = elementsM[iM + 1];
                                    if (elementsE[iE] >= tol) {
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
                                        elementsM[iM] = elem[0];
                                        elementsM[iM + 1] = elem[1];
                                    } else {
                                        elementsM[iM] = one[0];
                                        elementsM[iM + 1] = one[1];
                                    }
                                    iE += columnStrideE;
                                    iM += columnStrideM;
                                }
                                idxE += rowStrideE;
                                idxM += rowStrideM;
                            }
                        } else {
                            for (int r = startrow; r < stoprow; r++) {
                                for (int iE = idxE, iM = idxM, c = 0; c < cols; c++) {
                                    elem[0] = elementsM[iM];
                                    elem[1] = elementsM[iM + 1];
                                    if (elementsE[iE] >= tol) {
                                        elem = FComplex.inv(elem);
                                        elementsM[iM] = elem[0];
                                        elementsM[iM + 1] = elem[1];
                                    } else {
                                        elementsM[iM] = one[0];
                                        elementsM[iM + 1] = one[1];
                                    }
                                    iE += columnStrideE;
                                    iM += columnStrideM;
                                }
                                idxE += rowStrideE;
                                idxM += rowStrideM;
                            }
                        }
                    }
                });
            }
            try {
                for (int j = 0; j < np; j++) {
                    futures[j].get();
                }
            } catch (ExecutionException ex) {
                ex.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        } else {
            float[] elem = new float[2];
            int idxE = zeroE;
            int idxM = zeroM;
            if (maxE != 1.0) {
                for (int r = 0; r < rows; r++) {
                    for (int iE = idxE, iM = idxM, c = 0; c < cols; c++) {
                        elem[0] = elementsM[iM];
                        elem[1] = elementsM[iM + 1];
                        if (elementsE[iE] >= tol) {
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
                            elementsM[iM] = elem[0];
                            elementsM[iM + 1] = elem[1];
                        } else {
                            elementsM[iM] = one[0];
                            elementsM[iM + 1] = one[1];
                        }
                        iE += columnStrideE;
                        iM += columnStrideM;
                    }
                    idxE += rowStrideE;
                    idxM += rowStrideM;
                }
            } else {
                for (int r = 0; r < rows; r++) {
                    for (int iE = idxE, iM = idxM, c = 0; c < cols; c++) {
                        elem[0] = elementsM[iM];
                        elem[1] = elementsM[iM + 1];
                        if (elementsE[iE] >= tol) {
                            elem = FComplex.inv(elem);
                            elementsM[iM] = elem[0];
                            elementsM[iM + 1] = elem[1];
                        } else {
                            elementsM[iM] = one[0];
                            elementsM[iM + 1] = one[1];
                        }
                        iE += columnStrideE;
                        iM += columnStrideM;
                    }
                    idxE += rowStrideE;
                    idxM += rowStrideM;
                }
            }
        }
    }

    private float defaultTol2(FloatMatrix2D E, FloatMatrix2D B) {
        FloatMatrix1D s = new DenseFloatMatrix1D((int)E.size());
        System.arraycopy((float[]) E.elements(), 0, (float[]) s.elements(), 0, (int)s.size());
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
        AbstractMatrix2D Bhat = ((DenseFloatMatrix2D) B).getFft2();
        ((FComplexMatrix2D) Bhat).assign(FComplexFunctions.abs);
        Bhat = ((FComplexMatrix2D) Bhat).getRealPart();
        FloatMatrix1D bhat = new DenseFloatMatrix1D((int)Bhat.size(), (float[]) ((FloatMatrix2D) Bhat).elements(), 0, 1, false);
        bhat = bhat.viewSelection(indices);
        bhat.assign(FloatFunctions.div((float) Math.sqrt(B.size())));
        int n = (int)s.size();
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
