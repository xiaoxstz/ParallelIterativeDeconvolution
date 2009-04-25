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
import ij.ImageJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.WindowManager;
import ij.io.OpenDialog;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.PrintWriter;
import java.io.StringWriter;

import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.ToolTipManager;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;

import cern.colt.Arrays;
import cern.colt.Timer;
import cern.colt.matrix.tdouble.algo.solver.HyBRInnerSolver;
import cern.colt.matrix.tdouble.algo.solver.HyBRRegularizationMethod;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.Enums.PrecisionType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.MethodType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PreconditionerType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.ResizingType;
import edu.emory.mathcs.restoretools.iterative.cgls.CGLSDoubleIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.cgls.CGLSFloatIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.cgls.CGLSOptions;
import edu.emory.mathcs.restoretools.iterative.hybr.HyBRDoubleIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.hybr.HyBRFloatIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.hybr.HyBROptions;
import edu.emory.mathcs.restoretools.iterative.mrnsd.MRNSDDoubleIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.mrnsd.MRNSDFloatIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.mrnsd.MRNSDOptions;
import edu.emory.mathcs.restoretools.iterative.wpl.WPLDoubleIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.wpl.WPLFloatIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.wpl.WPLOptions;
import edu.emory.mathcs.utils.ConcurrencyUtils;

/**
 * Parallel Iterative Deconvolution 2D GUI
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class ParallelIterativeDeconvolution2D implements PlugIn, ImageListener {

    /**
     * Method used to deconvolve with WPL from a macro.
     * 
     * @param pathToBlurredImage
     *            path to a blurred image file
     * @param pathToPsf
     *            path to a PSF image file
     * @param pathToDeblurredImage
     *            path to a restored image file
     * @param boundaryStr
     *            type of boundary conditions
     * @param resizingStr
     *            type of resizing
     * @param outputStr
     *            type of the output image
     * @param precisionStr
     *            type of precision
     * @param thresholdStr
     *            threshold, set thresholdStr=-1 to disable the thresholding
     * @param maxItersStr
     *            maximal number of iterations
     * @param nOfThreadsStr
     *            maximal number of threads
     * @param showIterationsStr
     *            if true, then the restored image is displayed after each
     *            iteration
     * @param gammaStr
     *            regularization parameter for the Wiener Filter, set gammaStr=0
     *            to disable Wiener Filter
     * @param filterXYStr
     *            number of pixels in x and y directions for low-pass filter
     * @param normalizeStr
     *            if true, then PSF is normalized
     * @param logMeanStr
     *            if true, then the convergence information is displayed after
     *            each iteration
     * @param antiRingStr
     *            if true, then the anti-ringing step is performed
     * @param changeThreshPercentStr
     *            this parameter is used to stop the iteration if the image is
     *            not changing
     * @param dbStr
     *            if true, then all the data is in decibels
     * @param detectDivergenceStr
     *            if true, then the iterations are stopped when the changes
     *            appear to be increasing
     * @return path to deblurred image or error message
     */
    public static String deconvolveWPL(String pathToBlurredImage, String pathToPsf, String pathToDeblurredImage, String boundaryStr, String resizingStr, String outputStr, String precisionStr, String thresholdStr, String maxItersStr, String nOfThreadsStr, String showIterationsStr, String gammaStr,
            String filterXYStr, String normalizeStr, String logMeanStr, String antiRingStr, String changeThreshPercentStr, String dbStr, String detectDivergenceStr) {
        boolean showIterations, normalize, logMean, antiRing, db, detectDivergence;
        double threshold, gamma, filterXY, changeThreshPercent;
        int maxIters;
        int nOfThreads;
        BoundaryType boundary = null;
        ResizingType resizing = null;
        OutputType output = null;
        PrecisionType precision = null;
        ImagePlus imX = null;
        ImagePlus imB = IJ.openImage(pathToBlurredImage);
        if (imB == null) {
            return "Cannot open image " + pathToBlurredImage;
        }
        ImagePlus imPSF = IJ.openImage(pathToPsf);
        if (imPSF == null) {
            return "Cannot open image " + pathToPsf;
        }
        ImageProcessor ipB = imB.getProcessor();
        if (ipB instanceof ColorProcessor) {
            return "RGB images are not currently supported";
        }
        if (imB.getStackSize() > 1) {
            return "For 3D images use Parallel Iterative Deconvolution 3D";
        }
        ImageProcessor ipPSF = imPSF.getProcessor();
        if (ipPSF instanceof ColorProcessor) {
            return "RGB images are not currently supported";
        }
        if (imPSF.getStackSize() > 1) {
            return "For 3D images use Parallel Iterative Deconvolution 3D";
        }
        try {
            maxIters = Integer.parseInt(maxItersStr);
        } catch (Exception ex) {
            return "maxIters must be a positive integer";
        }
        if (maxIters < 1) {
            return "maxIters must be a positive integer";
        }
        for (BoundaryType elem : BoundaryType.values()) {
            if (elem.toString().equals(boundaryStr)) {
                boundary = elem;
                break;
            }
        }
        if (boundary == null) {
            return "boundary must be in " + Arrays.toString(BoundaryType.values());
        }
        for (ResizingType elem : ResizingType.values()) {
            if (elem.toString().equals(resizingStr)) {
                resizing = elem;
                break;
            }
        }
        if (resizing == null) {
            return "resizing must be in " + Arrays.toString(ResizingType.values());
        }
        for (OutputType elem : OutputType.values()) {
            if (elem.toString().equals(outputStr)) {
                output = elem;
                break;
            }
        }
        if (output == null) {
            return "output must be in " + Arrays.toString(OutputType.values());
        }
        for (PrecisionType elem : PrecisionType.values()) {
            if (elem.toString().equals(precisionStr)) {
                precision = elem;
                break;
            }
        }
        if (precision == null) {
            return "precision must be in " + Arrays.toString(PrecisionType.values());
        }
        try {
            threshold = Double.parseDouble(thresholdStr);
        } catch (Exception ex) {
            return "threshold must be a nonnegative number or -1 to disable";
        }
        if ((threshold != -1) && (threshold < 0)) {
            return "threshold must be a nonnegative number or -1 to disable";
        }
        try {
            nOfThreads = Integer.parseInt(nOfThreadsStr);
        } catch (Exception ex) {
            return "nOfThreads must be power of 2";
        }
        if (nOfThreads < 1) {
            return "nOfThreads must be power of 2";
        }
        if (!ConcurrencyUtils.isPowerOf2(nOfThreads)) {
            return "nOfThreads must be power of 2";
        }
        try {
            showIterations = Boolean.parseBoolean(showIterationsStr);
        } catch (Exception ex) {
            return "showItrations must be a boolean value (true or false)";
        }
        try {
            gamma = Double.parseDouble(gammaStr);
        } catch (Exception ex) {
            return "gamma must be a nonnegative value";
        }
        if (gamma < 0.0) {
            return "gamma must be a nonnegative value";
        }

        try {
            filterXY = Double.parseDouble(filterXYStr);
        } catch (Exception ex) {
            return "filterXY must be a nonnegative value";
        }
        if (filterXY < 0.0) {
            return "filterXY must be a nonnegative value";
        }
        try {
            normalize = Boolean.parseBoolean(normalizeStr);
        } catch (Exception ex) {
            return "normalize must be a boolean value (true or false)";
        }
        try {
            logMean = Boolean.parseBoolean(logMeanStr);
        } catch (Exception ex) {
            return "logMean must be a boolean value (true or false)";
        }
        try {
            antiRing = Boolean.parseBoolean(antiRingStr);
        } catch (Exception ex) {
            return "antiRing must be a boolean value (true or false)";
        }
        try {
            db = Boolean.parseBoolean(dbStr);
        } catch (Exception ex) {
            return "db must be a boolean value (true or false)";
        }
        try {
            detectDivergence = Boolean.parseBoolean(detectDivergenceStr);
        } catch (Exception ex) {
            return "detectDivergence must be a boolean value (true or false)";
        }
        try {
            changeThreshPercent = Double.parseDouble(changeThreshPercentStr);
        } catch (Exception ex) {
            return "changeThreshPercent must be a nonnegative value";
        }
        if (changeThreshPercent < 0.0) {
            IJ.error("changeThreshPercent must be a nonnegative value");
        }
        ConcurrencyUtils.setNumberOfThreads(nOfThreads);
        WPLOptions options = new WPLOptions(gamma, filterXY, 0, normalize, logMean, antiRing, changeThreshPercent, db, detectDivergence, (threshold == -1) ? false : true, threshold);
        switch (precision) {
        case DOUBLE:
            WPLDoubleIterativeDeconvolver2D dwpl = new WPLDoubleIterativeDeconvolver2D(imB, imPSF, boundary, resizing, output, maxIters, showIterations, options);
            imX = dwpl.deconvolve();
            break;
        case SINGLE:
            WPLFloatIterativeDeconvolver2D fwpl = new WPLFloatIterativeDeconvolver2D(imB, imPSF, boundary, resizing, output, maxIters, showIterations, options);
            imX = fwpl.deconvolve();
            break;
        }
        IJ.save(imX, pathToDeblurredImage);
        return pathToDeblurredImage;
    }

    /**
     * Method used to deconvolve with MRNSD from a macro.
     * 
     * @param pathToBlurredImage
     *            path to a blurred image file
     * @param pathToPsf
     *            path to a PSF image file
     * @param pathToDeblurredImage
     *            path to a restored image file
     * @param preconditionerStr
     *            type of preconditioner
     * @param preconditionerTolStr
     *            tolerance for the preconditioner
     * @param boundaryStr
     *            type of boundary conditions
     * @param resizingStr
     *            type of resizing
     * @param outputStr
     *            type of the output image
     * @param precisionStr
     *            type of precision
     * @param stoppingTolStr
     *            stopping tolerance, if stoppingTolStr==-1, then the tolerance
     *            is computed automatically
     * @param thresholdStr
     *            the smallest nonnegative value assigned to the restored image,
     *            if thresholdStr==-1, then the thresholding is disabled
     * @param logConvergenceStr
     *            if true, then the convergence information is displayed after
     *            each iteration
     * @param maxItersStr
     *            maximal number of iterations
     * @param nOfThreadsStr
     *            maximal number of threads
     * @param showIterationsStr
     *            if true, then the restored image is displayed after each
     *            iteration
     * @return path to deblurred image or error message
     */
    public static String deconvolveMRNSD(String pathToBlurredImage, String pathToPsf, String pathToDeblurredImage, String preconditionerStr, String preconditionerTolStr, String boundaryStr, String resizingStr, String outputStr, String precisionStr, String stoppingTolStr, String thresholdStr,
            String logConvergenceStr, String maxItersStr, String nOfThreadsStr, String showIterationsStr) {
        double stoppingTol;
        double preconditionerTol;
        boolean showIterations, logConvergence;
        double threshold;
        int maxIters;
        int nOfThreads;
        PreconditionerType preconditioner = null;
        BoundaryType boundary = null;
        ResizingType resizing = null;
        OutputType output = null;
        PrecisionType precision = null;
        ImagePlus imX = null;
        ImagePlus imB = IJ.openImage(pathToBlurredImage);
        if (imB == null) {
            return "Cannot open image " + pathToBlurredImage;
        }
        ImagePlus[][] imPSF = new ImagePlus[1][1];
        imPSF[0][0] = IJ.openImage(pathToPsf);
        if (imPSF[0][0] == null) {
            return "Cannot open image " + pathToPsf;
        }
        ImageProcessor ipB = imB.getProcessor();
        if (ipB instanceof ColorProcessor) {
            return "RGB images are not currently supported.";
        }
        if (imB.getStackSize() > 1) {
            return "For 3D images use Parallel Iterative Deconvolution 3D";
        }
        ImageProcessor ipPSF = imPSF[0][0].getProcessor();
        if (ipPSF instanceof ColorProcessor) {
            return "RGB images are not currently supported.";
        }
        if (imPSF[0][0].getStackSize() > 1) {
            return "For 3D images use Parallel Iterative Deconvolution 3D.";
        }
        try {
            maxIters = Integer.parseInt(maxItersStr);
        } catch (Exception ex) {
            return "nOfIters must be a positive integer";
        }
        if (maxIters < 1) {
            return "nOfIters must be a positive integer";
        }
        for (PreconditionerType elem : PreconditionerType.values()) {
            if (elem.toString().equals(preconditionerStr)) {
                preconditioner = elem;
                break;
            }
        }
        if (preconditioner == null) {
            return "preconditioner must be in " + Arrays.toString(PreconditionerType.values());
        }
        for (BoundaryType elem : BoundaryType.values()) {
            if (elem.toString().equals(boundaryStr)) {
                boundary = elem;
                break;
            }
        }
        if (boundary == null) {
            return "boundary must be in " + Arrays.toString(BoundaryType.values());
        }
        for (ResizingType elem : ResizingType.values()) {
            if (elem.toString().equals(resizingStr)) {
                resizing = elem;
                break;
            }
        }
        if (resizing == null) {
            return "resizing must be in " + Arrays.toString(ResizingType.values());
        }
        for (OutputType elem : OutputType.values()) {
            if (elem.toString().equals(outputStr)) {
                output = elem;
                break;
            }
        }
        if (output == null) {
            return "output must be in " + Arrays.toString(OutputType.values());
        }
        for (PrecisionType elem : PrecisionType.values()) {
            if (elem.toString().equals(precisionStr)) {
                precision = elem;
                break;
            }
        }
        if (precision == null) {
            return "precision must be in " + Arrays.toString(PrecisionType.values());
        }
        try {
            preconditionerTol = Double.parseDouble(preconditionerTolStr);
        } catch (Exception ex) {
            return "preconditionerTol must be a number between 0 and 1 or -1 for auto";
        }
        if ((preconditionerTol != -1) && ((preconditionerTol < 0) || (preconditionerTol > 1))) {
            return "preconditionerTol must be a number between 0 and 1 or -1 for auto";
        }
        try {
            stoppingTol = Double.parseDouble(stoppingTolStr);
        } catch (Exception ex) {
            return "stoppingTol must be a number between 0 and 1 or -1 for auto";
        }
        if ((stoppingTol != -1) && ((stoppingTol < 0) || (stoppingTol > 1))) {
            return "stoppingTol must be a number between 0 and 1 or -1 for auto";
        }
        try {
            threshold = Double.parseDouble(thresholdStr);
        } catch (Exception ex) {
            return "threshold must be a nonnegative number or -1 to disable";
        }
        if ((threshold != -1) && (threshold < 0)) {
            return "threshold must be a nonnegative number or -1 to disable";
        }
        try {
            logConvergence = Boolean.parseBoolean(logConvergenceStr);
        } catch (Exception ex) {
            return "logConvergence must be a boolean value (true or false)";
        }
        try {
            nOfThreads = Integer.parseInt(nOfThreadsStr);
        } catch (Exception ex) {
            return "nOfThreads must be power of 2";
        }
        if (nOfThreads < 1) {
            return "nOfThreads must be power of 2";
        }
        if (!ConcurrencyUtils.isPowerOf2(nOfThreads)) {
            return "nOfThreads must be power of 2";
        }
        try {
            showIterations = Boolean.parseBoolean(showIterationsStr);
        } catch (Exception ex) {
            return "showItrations must be a boolean value (true or false)";
        }
        ConcurrencyUtils.setNumberOfThreads(nOfThreads);
        MRNSDOptions options = new MRNSDOptions((stoppingTol == -1) ? true : false, stoppingTol, (threshold == -1) ? false : true, threshold, logConvergence);
        switch (precision) {
        case DOUBLE:
            MRNSDDoubleIterativeDeconvolver2D dmrnsd = new MRNSDDoubleIterativeDeconvolver2D(imB, imPSF, preconditioner, preconditionerTol, boundary, resizing, output, maxIters, showIterations, options);
            imX = dmrnsd.deconvolve();
            break;
        case SINGLE:
            MRNSDFloatIterativeDeconvolver2D fmrnsd = new MRNSDFloatIterativeDeconvolver2D(imB, imPSF, preconditioner, (float) preconditionerTol, boundary, resizing, output, maxIters, showIterations, options);
            imX = fmrnsd.deconvolve();
            break;
        }
        IJ.save(imX, pathToDeblurredImage);
        return pathToDeblurredImage;
    }

    /**
     * Method used to deconvolve with HyBR from a macro.
     * 
     * @param pathToBlurredImage
     *            path to a blurred image file
     * @param pathToPsf
     *            path to a PSF image file
     * @param pathToDeblurredImage
     *            path to a restored image file
     * @param preconditionerStr
     *            type of preconditioner
     * @param preconditionerTolStr
     *            tolerance for the preconditioner
     * @param boundaryStr
     *            type of boundary conditions
     * @param resizingStr
     *            type of resizing
     * @param outputStr
     *            type of the output image
     * @param precisionStr
     *            type of precision
     * @param thresholdStr
     *            the smallest nonnegative value assigned to the restored image,
     *            if thresholdStr==-1, then the thresholding is disabled
     * @param logConvergenceStr
     *            if true, then the convergence information is displayed after
     *            each iteration
     * @param maxItersStr
     *            maximal number of iterations
     * @param nOfThreadsStr
     *            maximal number of threads
     * @param showIterationsStr
     *            if true, then the restored image is displayed after each
     *            iteration
     * @param innerSolverStr
     *            type of an inner solver
     * @param regMethodStr
     *            type of a regularization
     * @param regParamStr
     *            regularization parameter
     * @param omegaStr
     *            omega parameter for weighted GCV
     * @param reorthStr
     *            if true, then the reorthogonalization of Lanczos subspaces is
     *            performed
     * @param beginRegStr
     *            begin regularization after this iteration
     * @param flatTolStr
     *            tolerance for detecting flatness in the GCV curve as a
     *            stopping criteria
     * @return path to deblurred image or error message
     */
    public static String deconvolveHyBR(String pathToBlurredImage, String pathToPsf, String pathToDeblurredImage, String preconditionerStr, String preconditionerTolStr, String boundaryStr, String resizingStr, String outputStr, String precisionStr, String thresholdStr, String logConvergenceStr,
            String maxItersStr, String nOfThreadsStr, String showIterationsStr, String innerSolverStr, String regMethodStr, String regParamStr, String omegaStr, String reorthStr, String beginRegStr, String flatTolStr) {
        double preconditionerTol, regParam, omega, flatTol;
        boolean showIteration, reorth, logConvergence;
        double threshold;
        int maxIters;
        int nOfThreads;
        int beginReg;
        PreconditionerType preconditioner = null;
        HyBRInnerSolver innerSolver = null;
        HyBRRegularizationMethod regMethod = null;
        BoundaryType boundary = null;
        ResizingType resizing = null;
        OutputType output = null;
        PrecisionType precision = null;
        ImagePlus imX = null;
        ImagePlus imB = IJ.openImage(pathToBlurredImage);
        if (imB == null) {
            return "Cannot open image " + pathToBlurredImage;
        }
        ImagePlus[][] imPSF = new ImagePlus[1][1];
        imPSF[0][0] = IJ.openImage(pathToPsf);
        if (imPSF[0][0] == null) {
            return "Cannot open image " + pathToPsf;
        }
        ImageProcessor ipB = imB.getProcessor();
        if (ipB instanceof ColorProcessor) {
            return "RGB images are not currently supported.";
        }
        if (imB.getStackSize() > 1) {
            return "For 3D images use Parallel Iterative Deconvolution 3D";
        }
        ImageProcessor ipPSF = imPSF[0][0].getProcessor();
        if (ipPSF instanceof ColorProcessor) {
            return "RGB images are not currently supported.";
        }
        if (imPSF[0][0].getStackSize() > 1) {
            return "For 3D images use Parallel Iterative Deconvolution 3D.";
        }
        try {
            maxIters = Integer.parseInt(maxItersStr);
        } catch (Exception ex) {
            return "maxIters must be a positive integer";
        }
        if (maxIters < 1) {
            return "maxIters must be a positive integer";
        }
        for (PreconditionerType elem : PreconditionerType.values()) {
            if (elem.toString().equals(preconditionerStr)) {
                preconditioner = elem;
                break;
            }
        }
        if (preconditioner == null) {
            return "preconditioner must be in " + Arrays.toString(PreconditionerType.values());
        }
        for (BoundaryType elem : BoundaryType.values()) {
            if (elem.toString().equals(boundaryStr)) {
                boundary = elem;
                break;
            }
        }
        if (boundary == null) {
            return "boundary must be in " + Arrays.toString(BoundaryType.values());
        }
        for (ResizingType elem : ResizingType.values()) {
            if (elem.toString().equals(resizingStr)) {
                resizing = elem;
                break;
            }
        }
        if (resizing == null) {
            return "resizing must be in " + Arrays.toString(ResizingType.values());
        }
        for (OutputType elem : OutputType.values()) {
            if (elem.toString().equals(outputStr)) {
                output = elem;
                break;
            }
        }
        if (output == null) {
            return "output must be in " + Arrays.toString(OutputType.values());
        }
        for (PrecisionType elem : PrecisionType.values()) {
            if (elem.toString().equals(precisionStr)) {
                precision = elem;
                break;
            }
        }
        if (precision == null) {
            return "precision must be in " + Arrays.toString(PrecisionType.values());
        }
        try {
            preconditionerTol = Double.parseDouble(preconditionerTolStr);
        } catch (Exception ex) {
            return "preconditionerTol must be a number between 0 and 1 or -1 for auto";
        }
        if ((preconditionerTol != -1) && ((preconditionerTol < 0) || (preconditionerTol > 1))) {
            return "preconditionerTol must be a number between 0 and 1 or -1 for auto";
        }
        try {
            threshold = Double.parseDouble(thresholdStr);
        } catch (Exception ex) {
            return "threshold must be a nonnegative number or -1 to disable";
        }
        if ((threshold != -1) && (threshold < 0)) {
            return "threshold must be a nonnegative number or -1 to disable";
        }
        try {
            logConvergence = Boolean.parseBoolean(logConvergenceStr);
        } catch (Exception ex) {
            return "logConvergence must be a boolean value (true or false)";
        }
        try {
            nOfThreads = Integer.parseInt(nOfThreadsStr);
        } catch (Exception ex) {
            return "nOfThreads must be power of 2";
        }
        if (nOfThreads < 1) {
            return "nOfThreads must be power of 2";
        }
        if (!ConcurrencyUtils.isPowerOf2(nOfThreads)) {
            return "nOfThreads must be power of 2";
        }
        try {
            showIteration = Boolean.parseBoolean(showIterationsStr);
        } catch (Exception ex) {
            return "showItration must be a boolean value (true or false)";
        }
        for (HyBRInnerSolver elem : HyBRInnerSolver.values()) {
            if (elem.toString().equals(innerSolverStr)) {
                innerSolver = elem;
                break;
            }
        }
        if (innerSolver == null) {
            return "innerSolver must be in " + Arrays.toString(HyBRInnerSolver.values());
        }
        for (HyBRRegularizationMethod elem : HyBRRegularizationMethod.values()) {
            if (elem.toString().equals(regMethodStr)) {
                regMethod = elem;
                break;
            }
        }
        if (regMethod == null) {
            return "regMethod method must be in " + Arrays.toString(HyBRRegularizationMethod.values());
        }

        try {
            regParam = Double.parseDouble(regParamStr);
        } catch (Exception ex) {
            return "regParam must be a floating-point number between 0 and 1 or -1 for auto";
        }
        if ((regParam != -1) && ((regParam < 0.0) || (regParam > 1.0))) {
            return "regParam must be a floating-point number between 0 and 1 or -1 for auto";
        }

        try {
            omega = Double.parseDouble(omegaStr);
        } catch (Exception ex) {
            return "omega must be a nonnegative floating-point number";
        }
        if (omega < 0.0) {
            return "omega must be a nonnegative floating-point number";
        }
        try {
            reorth = Boolean.parseBoolean(reorthStr);
        } catch (Exception ex) {
            return "reorth must be a boolean value (true or false)";
        }
        try {
            beginReg = Integer.parseInt(beginRegStr);
        } catch (Exception ex) {
            return "beginReg must be an integer number greater than 1";
        }
        if (beginReg <= 1) {
            return "beginReg must be an integer number greater than 1";
        }

        try {
            flatTol = Double.parseDouble(flatTolStr);
        } catch (Exception ex) {
            return "flatTol must be a nonnegative floating-point number";
        }
        if (flatTol < 0.0) {
            return "flatTol must be a nonnegative floating-point number";
        }

        ConcurrencyUtils.setNumberOfThreads(nOfThreads);
        HyBROptions options = new HyBROptions(innerSolver, regMethod, regParam, omega, reorth, beginReg, flatTol, logConvergence, (threshold == -1) ? false : true, threshold);
        switch (precision) {
        case DOUBLE:
            HyBRDoubleIterativeDeconvolver2D dhybr = new HyBRDoubleIterativeDeconvolver2D(imB, imPSF, preconditioner, preconditionerTol, boundary, resizing, output, maxIters, showIteration, options);
            imX = dhybr.deconvolve();
            break;
        case SINGLE:
            HyBRFloatIterativeDeconvolver2D fhybr = new HyBRFloatIterativeDeconvolver2D(imB, imPSF, preconditioner, (float) preconditionerTol, boundary, resizing, output, maxIters, showIteration, options);
            imX = fhybr.deconvolve();
            break;
        }
        IJ.save(imX, pathToDeblurredImage);
        return pathToDeblurredImage;
    }

    /**
     * Method used to deconvolve with CGLS from a macro.
     * 
     * @param pathToBlurredImage
     *            path to a blurred image file
     * @param pathToPsf
     *            path to a PSF image file
     * @param pathToDeblurredImage
     *            path to a restored image file
     * @param preconditionerStr
     *            type of preconditioner
     * @param preconditionerTolStr
     *            tolerance for the preconditioner
     * @param boundaryStr
     *            type of boundary conditions
     * @param resizingStr
     *            type of resizing
     * @param outputStr
     *            type of the output image
     * @param precisionStr
     *            type of precision
     * @param stoppingTolStr
     *            stopping tolerance, if stoppingTolStr==-1, then the tolerance
     *            is computed automatically
     * @param thresholdStr
     *            the smallest nonnegative value assigned to the restored image,
     *            if thresholdStr==-1, then the thresholding is disabled
     * @param logConvergenceStr
     *            if true, then the convergence information is displayed after
     *            each iteration
     * @param maxItersStr
     *            maximal number of iterations
     * @param nOfThreadsStr
     *            maximal number of threads
     * @param showIterationsStr
     *            if true, then the restored image is displayed after each
     *            iteration
     * @return path to deblurred image or error message
     */
    public static String deconvolveCGLS(String pathToBlurredImage, String pathToPsf, String pathToDeblurredImage, String preconditionerStr, String preconditionerTolStr, String boundaryStr, String resizingStr, String outputStr, String precisionStr, String stoppingTolStr, String thresholdStr,
            String logConvergenceStr, String maxItersStr, String nOfThreadsStr, String showIterationsStr) {
        double stoppingTol;
        double preconditionerTol;
        boolean showIterations, logConvergence;
        double threshold;
        int maxIters;
        int nOfThreads;
        PreconditionerType preconditioner = null;
        BoundaryType boundary = null;
        ResizingType resizing = null;
        OutputType output = null;
        PrecisionType precision = null;
        ImagePlus imX = null;
        ImagePlus imB = IJ.openImage(pathToBlurredImage);
        if (imB == null) {
            return "Cannot open image " + pathToBlurredImage;
        }
        ImagePlus[][] imPSF = new ImagePlus[1][1];
        imPSF[0][0] = IJ.openImage(pathToPsf);
        if (imPSF[0][0] == null) {
            return "Cannot open image " + pathToPsf;
        }
        ImageProcessor ipB = imB.getProcessor();
        if (ipB instanceof ColorProcessor) {
            return "RGB images are not currently supported.";
        }
        if (imB.getStackSize() > 1) {
            return "For 3D images use Parallel Iterative Deconvolution 3D.";
        }
        ImageProcessor ipPSF = imPSF[0][0].getProcessor();
        if (ipPSF instanceof ColorProcessor) {
            return "RGB images are not currently supported.";
        }
        if (imPSF[0][0].getStackSize() > 1) {
            return "For 3D images use Parallel Iterative Deconvolution 3D.";
        }
        try {
            maxIters = Integer.parseInt(maxItersStr);
        } catch (Exception ex) {
            return "nOfIters must be a positive integer";
        }
        if (maxIters < 1) {
            return "nOfIters must be a positive integer";
        }
        for (PreconditionerType elem : PreconditionerType.values()) {
            if (elem.toString().equals(preconditionerStr)) {
                preconditioner = elem;
                break;
            }
        }
        if (preconditioner == null) {
            return "preconditioner must be in " + Arrays.toString(PreconditionerType.values());
        }
        for (BoundaryType elem : BoundaryType.values()) {
            if (elem.toString().equals(boundaryStr)) {
                boundary = elem;
                break;
            }
        }
        if (boundary == null) {
            return "boundary must be in " + Arrays.toString(BoundaryType.values());
        }
        for (ResizingType elem : ResizingType.values()) {
            if (elem.toString().equals(resizingStr)) {
                resizing = elem;
                break;
            }
        }
        if (resizing == null) {
            return "resizing must be in " + Arrays.toString(ResizingType.values());
        }
        for (OutputType elem : OutputType.values()) {
            if (elem.toString().equals(outputStr)) {
                output = elem;
                break;
            }
        }
        if (output == null) {
            return "output must be in " + Arrays.toString(OutputType.values());
        }
        for (PrecisionType elem : PrecisionType.values()) {
            if (elem.toString().equals(precisionStr)) {
                precision = elem;
                break;
            }
        }
        if (precision == null) {
            return "precision must be in " + Arrays.toString(PrecisionType.values());
        }
        try {
            preconditionerTol = Double.parseDouble(preconditionerTolStr);
        } catch (Exception ex) {
            return "preconditionerTol must be a number between 0 and 1 or -1 for auto";
        }
        if ((preconditionerTol != -1) && ((preconditionerTol < 0) || (preconditionerTol > 1))) {
            return "preconditionerTol must be a number between 0 and 1 or -1 for auto";
        }
        try {
            stoppingTol = Double.parseDouble(stoppingTolStr);
        } catch (Exception ex) {
            return "stoppingTol must be a number between 0 and 1 or -1 for auto";
        }
        if ((stoppingTol != -1) && ((stoppingTol < 0) || (stoppingTol > 1))) {
            return "stoppingTol must be a number between 0 and 1 or -1 for auto";
        }
        try {
            threshold = Double.parseDouble(thresholdStr);
        } catch (Exception ex) {
            return "threshold must be a nonnegative number or -1 to disable";
        }
        if ((threshold != -1) && (threshold < 0)) {
            return "threshold must be a nonnegative number or -1 to disable";
        }
        try {
            logConvergence = Boolean.parseBoolean(logConvergenceStr);
        } catch (Exception ex) {
            return "logConvergence must be a boolean value (true or false)";
        }
        try {
            nOfThreads = Integer.parseInt(nOfThreadsStr);
        } catch (Exception ex) {
            return "nOfThreads must be power of 2";
        }
        if (nOfThreads < 1) {
            return "nOfThreads must be power of 2";
        }
        if (!ConcurrencyUtils.isPowerOf2(nOfThreads)) {
            return "nOfThreads must be power of 2";
        }
        try {
            showIterations = Boolean.parseBoolean(showIterationsStr);
        } catch (Exception ex) {
            return "showItration must be a boolean value (true or false)";
        }
        ConcurrencyUtils.setNumberOfThreads(nOfThreads);
        CGLSOptions options = new CGLSOptions((stoppingTol == -1) ? true : false, stoppingTol, (threshold == -1) ? false : true, threshold, logConvergence);
        switch (precision) {
        case DOUBLE:
            CGLSDoubleIterativeDeconvolver2D dcgls = new CGLSDoubleIterativeDeconvolver2D(imB, imPSF, preconditioner, preconditionerTol, boundary, resizing, output, maxIters, showIterations, options);
            imX = dcgls.deconvolve();
            break;
        case SINGLE:
            CGLSFloatIterativeDeconvolver2D fcgls = new CGLSFloatIterativeDeconvolver2D(imB, imPSF, preconditioner, (float) preconditionerTol, boundary, resizing, output, maxIters, showIterations, options);
            imX = fcgls.deconvolve();
            break;
        }
        IJ.save(imX, pathToDeblurredImage);
        return pathToDeblurredImage;
    }

    private final static String version = "1.9";

    private static final String[] methodNames = { "MRNSD", "WPL", "CGLS", "HyBR" };

    private static final String[] shortMethodNames = { "mrnsd", "wpl", "cgls", "hybr" };

    private static final String[] precisionNames = { "Single", "Double" };

    private static final String[] precondNames = { "FFT preconditioner", "None" };

    private static final String[] shortPrecondNames = { "fft", "none" };

    private static final String[] boundaryNames = { "Reflexive", "Periodic", "Zero" };

    private static final String[] resizingNames = { "Auto", "Minimal", "Next power of two" };

    private static final String[] outputNames = { "Same as source", "Byte (8-bit)", "Short (16-bit)", "Float (32-bit)" };

    private static final String[] shortBoundaryNames = { "ref", "per", "zero" };

    private JFrame mainPanel, psfCreatePanel, psfEditPanel, hybrOptionsPanel, mrnsdOptionsPanel, cglsOptionsPanel, wplOptionsFrame;

    private CGLSDoubleIterativeDeconvolver2D dcgls;

    private HyBRDoubleIterativeDeconvolver2D dhybr;

    private MRNSDDoubleIterativeDeconvolver2D dmrnsd;

    private WPLDoubleIterativeDeconvolver2D dwpl;

    private CGLSFloatIterativeDeconvolver2D fcgls;

    private HyBRFloatIterativeDeconvolver2D fhybr;

    private MRNSDFloatIterativeDeconvolver2D fmrnsd;

    private WPLFloatIterativeDeconvolver2D fwpl;

    private ImagePlus imB, imX;

    private ImagePlus[][] imPSF;

    private int[] windowIDs;

    private String[] imageTitles;

    private JComboBox blurChoice, psfChoice, methodChoice, precondChoice, boundaryChoice, resizingChoice, outputChoice, precisionChoice;

    private JTextField itersField, threadsField, precondField;

    private JCheckBox variantPSFCheck, itersCheck, precondCheck;

    private JButton definePSFButton, editPSFButton, optionsButton, deconvolveButton, cancelButton;

    private int psfRows, psfColumns;

    private MRNSDOptions mrnsdOptions = new MRNSDOptions();

    private boolean mrnsdOptionsSet = false;

    private CGLSOptions cglsOptions = new CGLSOptions();

    private boolean cglsOptionsSet = false;

    private HyBROptions hybrOptions = new HyBROptions();

    private boolean hybrOptionsSet = false;

    private WPLOptions wplOptions = new WPLOptions();

    private boolean wplOptionsSet = false;

    private int maxIters, threads;

    private double precTol;

    private ImageListener getImageListener() {
        return this;
    }

    public void run(String arg) {
        if (IJ.versionLessThan("1.35l")) {
            IJ.showMessage("This plugin requires ImageJ 1.35l+");
            return;
        }

        if (!IJ.isJava15()) {
            IJ.showMessage("This plugin requires Sun Java 1.5+");
            return;
        }
        WindowManager.checkForDuplicateName = true;
        ImagePlus.addImageListener(this);
        mainPanel = new MainPanel("Parallel Iterative Deconvolution 2D " + version + " ");
    }

    public void imageClosed(ImagePlus imp) {
        blurChoice.removeItem(imp.getTitle());
        blurChoice.revalidate();
        psfChoice.removeItem(imp.getTitle());
        psfChoice.revalidate();
        if (imX != null) {
            if (imp.getTitle().equals(imX.getTitle())) {
                clean_old_data();
            }
        }
    }

    public void imageOpened(ImagePlus imp) {
        blurChoice.addItem(imp.getTitle());
        blurChoice.revalidate();
        psfChoice.addItem(imp.getTitle());
        psfChoice.revalidate();

    }

    public void imageUpdated(ImagePlus imp) {
    }

    private void clean_old_data() {
        dcgls = null;
        dmrnsd = null;
        dhybr = null;
        dwpl = null;
        fcgls = null;
        fmrnsd = null;
        fhybr = null;
        fwpl = null;
        System.gc();
    }

    private void clean_all() {
        clean_old_data();
        imB = null;
        imPSF = null;
        imX = null;
        windowIDs = null;
        imageTitles = null;
        System.gc();
    }

    private class MRNSDOptionsPanel extends JFrame {

        private static final long serialVersionUID = -8230645552596239730L;

        private JTextField stoppingTolField, thresholdField;

        private JCheckBox stoppingTolCheck, thresholdCheck, logConvergenceCheck;

        public MRNSDOptionsPanel(String name) {
            super(name);
            init();
        }

        private void init() {
            Container pane = getContentPane();
            pane.setLayout(new BoxLayout(pane, BoxLayout.Y_AXIS));

            JPanel stoppingTolPanel = new JPanel();
            stoppingTolPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel stoppingTolLabel = new JLabel("Stopping tolerance:");
            stoppingTolPanel.add(stoppingTolLabel);
            stoppingTolField = new JTextField(new Double(mrnsdOptions.getStoppingTol()).toString(), 5);
            stoppingTolField.addActionListener(new StoppingTolFieldActionListener());
            stoppingTolField.setEnabled(false);
            stoppingTolPanel.add(stoppingTolField);
            stoppingTolCheck = new JCheckBox("Auto");
            stoppingTolCheck.addItemListener(new StoppingTolCheckItemListener());
            stoppingTolCheck.setSelected(mrnsdOptions.getAutoStoppingTol());
            stoppingTolPanel.add(stoppingTolCheck);
            pane.add(stoppingTolPanel);

            JPanel thresholdPanel = new JPanel();
            thresholdPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            thresholdCheck = new JCheckBox("Threshold:  ");
            thresholdCheck.setSelected(mrnsdOptions.getUseThreshold());
            thresholdCheck.addItemListener(new ThresholdCheckItemListener());
            thresholdPanel.add(thresholdCheck);
            thresholdField = new JTextField(new Double(mrnsdOptions.getThreshold()).toString(), 6);
            thresholdField.setEnabled(mrnsdOptions.getUseThreshold());
            thresholdField.addActionListener(new ThresholdFieldActionListener());
            thresholdPanel.add(thresholdField);
            pane.add(thresholdPanel);

            JPanel logConvergencePanel = new JPanel();
            logConvergencePanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            logConvergenceCheck = new JCheckBox("Log convergence");
            logConvergenceCheck.setSelected(mrnsdOptions.getLogConvergence());
            logConvergencePanel.add(logConvergenceCheck);
            pane.add(logConvergencePanel);

            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
            JButton okButton = new JButton("OK");
            okButton.addActionListener(new OkButtonActionListener());
            JButton cancelButton = new JButton("Cancel");
            cancelButton.addActionListener(new CancelButtonActionListener());
            buttonPanel.add(okButton);
            buttonPanel.add(cancelButton);
            pane.add(buttonPanel);
            validate();
            setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
            setResizable(false);
            pack();
            setLocationRelativeTo(null);
            setVisible(true);
        }

        private class ThresholdCheckItemListener implements ItemListener {
            public void itemStateChanged(ItemEvent e) {
                if (thresholdCheck.isSelected()) {
                    thresholdField.setEnabled(true);
                } else {
                    thresholdField.setEnabled(false);
                }
            }
        }

        private class StoppingTolCheckItemListener implements ItemListener {
            public void itemStateChanged(ItemEvent e) {
                if (stoppingTolCheck.isSelected()) {
                    stoppingTolField.setEnabled(false);
                } else {
                    stoppingTolField.setEnabled(true);
                }
            }
        }

        private class StoppingTolFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double stoppingTol = 0;
                try {
                    stoppingTol = Double.parseDouble(stoppingTolField.getText());
                } catch (Exception ex) {
                    IJ.error("Stopping tolerance must be a nonnegative floating-point number.");
                }
                if (stoppingTol < 0.0) {
                    IJ.error("Stopping tolerance must be a nonnegative floating-point number.");
                }
            }
        }

        private class ThresholdFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double threshold = 0;
                try {
                    threshold = Double.parseDouble(thresholdField.getText());
                } catch (Exception ex) {
                    IJ.error("Threshold must be a nonnegative value.");
                }
                if (threshold < 0.0) {
                    IJ.error("Threshold must be a nonnegative value.");
                }
            }
        }

        private boolean assignOptions(MRNSDOptions mrnsdOptions) {
            double stoppingTol = 0;
            if (!stoppingTolCheck.isSelected()) {
                try {
                    stoppingTol = Double.parseDouble(stoppingTolField.getText());
                } catch (Exception ex) {
                    IJ.error("Stopping tolerance must be a nonnegative floating-point number.");
                    return false;
                }
                if (stoppingTol < 0.0) {
                    IJ.error("Stopping tolerance must be a nonnegative floating-point number.");
                    return false;
                }
            }
            double threshold = 0;
            if (!thresholdCheck.isSelected()) {
                try {
                    threshold = Double.parseDouble(thresholdField.getText());
                } catch (Exception ex) {
                    IJ.error("Threshold must be a nonnegative value.");
                    return false;
                }
                if (threshold < 0.0) {
                    IJ.error("Threshold must be a nonnegative value.");
                    return false;
                }
            }
            mrnsdOptions.setAutoStoppingTol(stoppingTolCheck.isSelected());
            mrnsdOptions.setStoppingTol(stoppingTol);
            mrnsdOptions.setUseThreshold(thresholdCheck.isSelected());
            mrnsdOptions.setThreshold(threshold);
            mrnsdOptions.setLogConvergence(logConvergenceCheck.isSelected());
            return true;
        }

        private class OkButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                if (assignOptions(mrnsdOptions) == true) {
                    mrnsdOptionsPanel.setVisible(false);
                    mrnsdOptionsSet = true;
                }
            }
        }

        private class CancelButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                mrnsdOptionsPanel.dispose();
                mrnsdOptionsSet = false;
            }
        }

    }

    private class CGLSOptionsPanel extends JFrame {

        private static final long serialVersionUID = -8230645552596239730L;

        private JTextField stoppingTolField, thresholdField;

        private JCheckBox stoppingTolCheck, thresholdCheck, logConvergenceCheck;

        public CGLSOptionsPanel(String name) {
            super(name);
            init();
        }

        private void init() {
            Container pane = getContentPane();
            pane.setLayout(new BoxLayout(pane, BoxLayout.Y_AXIS));

            JPanel stoppingTolPanel = new JPanel();
            stoppingTolPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel stoppingTolLabel = new JLabel("Stopping tolerance:");
            stoppingTolPanel.add(stoppingTolLabel);
            stoppingTolField = new JTextField(new Double(cglsOptions.getStoppingTol()).toString(), 5);
            stoppingTolField.addActionListener(new StoppingTolFieldActionListener());
            stoppingTolField.setEnabled(false);
            stoppingTolPanel.add(stoppingTolField);
            stoppingTolCheck = new JCheckBox("Auto");
            stoppingTolCheck.addItemListener(new StoppingTolCheckItemListener());
            stoppingTolCheck.setSelected(cglsOptions.getAutoStoppingTol());
            stoppingTolPanel.add(stoppingTolCheck);
            pane.add(stoppingTolPanel);

            JPanel thresholdPanel = new JPanel();
            thresholdPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            thresholdCheck = new JCheckBox("Threshold:  ");
            thresholdCheck.setSelected(cglsOptions.getUseThreshold());
            thresholdCheck.addItemListener(new ThresholdCheckItemListener());
            thresholdPanel.add(thresholdCheck);
            thresholdField = new JTextField(new Double(cglsOptions.getThreshold()).toString(), 6);
            thresholdField.setEnabled(cglsOptions.getUseThreshold());
            thresholdField.addActionListener(new ThresholdFieldActionListener());
            thresholdPanel.add(thresholdField);
            pane.add(thresholdPanel);

            JPanel logConvergencePanel = new JPanel();
            logConvergencePanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            logConvergenceCheck = new JCheckBox("Log convergence");
            logConvergenceCheck.setSelected(cglsOptions.getLogConvergence());
            logConvergencePanel.add(logConvergenceCheck);
            add(logConvergencePanel);

            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
            JButton okButton = new JButton("OK");
            okButton.addActionListener(new OkButtonActionListener());
            JButton cancelButton = new JButton("Cancel");
            cancelButton.addActionListener(new CancelButtonActionListener());
            buttonPanel.add(okButton);
            buttonPanel.add(cancelButton);
            pane.add(buttonPanel);
            validate();
            validate();
            setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
            setResizable(false);
            pack();
            setLocationRelativeTo(null);
            setVisible(true);
        }

        private class ThresholdCheckItemListener implements ItemListener {
            public void itemStateChanged(ItemEvent e) {
                if (thresholdCheck.isSelected()) {
                    thresholdField.setEnabled(true);
                } else {
                    thresholdField.setEnabled(false);
                }
            }
        }

        private class StoppingTolCheckItemListener implements ItemListener {
            public void itemStateChanged(ItemEvent e) {
                if (stoppingTolCheck.isSelected()) {
                    stoppingTolField.setEnabled(false);
                } else {
                    stoppingTolField.setEnabled(true);
                }
            }
        }

        private class StoppingTolFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double stoppingTol = 0;
                try {
                    stoppingTol = Double.parseDouble(stoppingTolField.getText());
                } catch (Exception ex) {
                    IJ.error("Stopping tolerance must be a nonnegative floating-point number.");
                }
                if (stoppingTol < 0.0) {
                    IJ.error("Stopping tolerance must be a nonnegative floating-point number.");
                }
            }
        }

        private class ThresholdFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double threshold = 0;
                try {
                    threshold = Double.parseDouble(thresholdField.getText());
                } catch (Exception ex) {
                    IJ.error("Threshold must be a nonnegative value.");
                }
                if (threshold < 0.0) {
                    IJ.error("Threshold must be a nonnegative value.");
                }
            }
        }

        private boolean assignOptions(CGLSOptions cglsOptions) {
            double stoppingTol = 0;
            if (!stoppingTolCheck.isSelected()) {
                try {
                    stoppingTol = Double.parseDouble(stoppingTolField.getText());
                } catch (Exception ex) {
                    IJ.error("Stopping tolerance must be a nonnegative floating-point number.");
                    return false;
                }
                if (stoppingTol < 0.0) {
                    IJ.error("Stopping tolerance must be a nonnegative floating-point number.");
                    return false;
                }
            }
            double threshold = 0;
            if (!thresholdCheck.isSelected()) {
                try {
                    threshold = Double.parseDouble(thresholdField.getText());
                } catch (Exception ex) {
                    IJ.error("Threshold must be a nonnegative value.");
                    return false;
                }
                if (threshold < 0.0) {
                    IJ.error("Threshold must be a nonnegative value.");
                    return false;
                }
            }
            cglsOptions.setAutoStoppingTol(stoppingTolCheck.isSelected());
            cglsOptions.setStoppingTol(stoppingTol);
            cglsOptions.setUseThreshold(thresholdCheck.isSelected());
            cglsOptions.setThreshold(threshold);
            cglsOptions.setLogConvergence(logConvergenceCheck.isSelected());
            return true;
        }

        private class OkButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                if (assignOptions(cglsOptions) == true) {
                    cglsOptionsPanel.setVisible(false);
                    cglsOptionsSet = true;
                }
            }
        }

        private class CancelButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                cglsOptionsPanel.dispose();
                cglsOptionsSet = false;
            }
        }
    }

    private class HyBROptionsPanel extends JFrame {

        private static final long serialVersionUID = 1L;

        private JComboBox solverChoice, regChoice;

        private JTextField regField, begRegField, omegaField, flatTolField, thresholdField;

        private JCheckBox reorthCheck, thresholdCheck, logConvergenceCheck;

        private final String[] innerSolverNames = { "Tikhonov", "None" };

        private final String[] regMethodNames = { "GCV", "WGCV", "Adaptive WGCV", "None" };

        public HyBROptionsPanel(String name) {
            super(name);
            init();
        }

        private void init() {
            Container pane = getContentPane();
            pane.setLayout(new BoxLayout(pane, BoxLayout.Y_AXIS));

            JPanel regPanel = new JPanel();
            regPanel.setLayout(new BoxLayout(regPanel, BoxLayout.Y_AXIS));
            Border border = new TitledBorder(null, "Regularization options", TitledBorder.LEFT, TitledBorder.TOP);
            regPanel.setBorder(border);

            JPanel regPanelFirstLine = new JPanel();
            regPanelFirstLine.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel regMethodLabel = new JLabel("Method:");
            regPanelFirstLine.add(regMethodLabel);
            regChoice = new JComboBox(regMethodNames);
            regChoice.setSelectedIndex(hybrOptions.getRegMethod().ordinal());
            regChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            regChoice.addActionListener(new RegMethodChoiceActionListener());
            regChoice
                    .setToolTipText("<html>Choose a method for computing a regularization parameter:<br><ul><li>GCV - Generalized Cross-Validation</li><li>WGCV - Weighted Generalized Cross-Validation</li><li>Adaptive WGCV - Adaptive Weighted Generalized Cross-Validation</li><li>None - no method for choosing regularization parameter.<br> You have to enter your own value.</li></ul></html>");
            regPanelFirstLine.add(regChoice);
            JLabel regParamLabel = new JLabel("Parameter:");
            regPanelFirstLine.add(regParamLabel);
            regField = new JTextField(new Double(hybrOptions.getRegParameter()).toString(), 5);
            regField.setToolTipText("<html>Regularization parameter</html>");
            regField.setEnabled(false);
            regField.addActionListener(new RegFieldActionListener());
            regPanelFirstLine.add(regField);
            regPanel.add(regPanelFirstLine);

            JPanel regPanelThirdLine = new JPanel();
            regPanelThirdLine.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel omegaLabel = new JLabel("Omega:");
            regPanelThirdLine.add(omegaLabel);
            omegaField = new JTextField(new Double(hybrOptions.getOmega()).toString(), 5);
            omegaField.setToolTipText("<html>Omega parameter</html>");
            omegaField.setEnabled(false);
            omegaField.addActionListener(new OmegaFieldActionListener());
            regPanelThirdLine.add(omegaField);
            regPanel.add(regPanelThirdLine);

            JPanel regPanelForthLine = new JPanel();
            regPanelForthLine.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel begRegLabel = new JLabel("Begin regularization after this iteration:");
            regPanelForthLine.add(begRegLabel);
            begRegField = new JTextField(new Integer(hybrOptions.getBeginReg()).toString(), 2);
            begRegField.addActionListener(new BegRegFieldActionListener());
            regPanelForthLine.add(begRegField);
            regPanel.add(regPanelForthLine);

            pane.add(regPanel);

            JPanel solverPanel = new JPanel();
            solverPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel solverLabel = new JLabel("Inner solver:");
            solverPanel.add(solverLabel);
            solverChoice = new JComboBox(innerSolverNames);
            solverChoice.setSelectedIndex(hybrOptions.getInnerSolver().ordinal());
            solverChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            solverPanel.add(solverChoice);
            pane.add(solverPanel);

            JPanel flatTolPanel = new JPanel();
            flatTolPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel flatTolLabel = new JLabel("Stopping tolerance:");
            flatTolPanel.add(flatTolLabel);
            flatTolField = new JTextField(new Double(hybrOptions.getFlatTolerance()).toString(), 5);
            flatTolField.setToolTipText("<html>Tolerance for detecting flatness in <br>the GCV curve as a stopping criteria</html>");
            flatTolField.addActionListener(new FlatTolFieldActionListener());
            flatTolPanel.add(flatTolField);
            pane.add(flatTolPanel);

            JPanel reorthPanel = new JPanel();
            reorthPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            reorthCheck = new JCheckBox("Reorthogonalize Lanczos subspaces");
            reorthCheck.setSelected(hybrOptions.getReorthogonalize());
            reorthPanel.add(reorthCheck);
            add(reorthPanel);

            JPanel thresholdPanel = new JPanel();
            thresholdPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            thresholdCheck = new JCheckBox("Threshold:  ");
            thresholdCheck.setSelected(hybrOptions.getUseThreshold());
            thresholdCheck.addItemListener(new ThresholdCheckItemListener());
            thresholdPanel.add(thresholdCheck);
            thresholdField = new JTextField(new Double(hybrOptions.getThreshold()).toString(), 6);
            thresholdField.setEnabled(hybrOptions.getUseThreshold());
            thresholdField.addActionListener(new ThresholdFieldActionListener());
            thresholdPanel.add(thresholdField);
            pane.add(thresholdPanel);

            JPanel logConvergencePanel = new JPanel();
            logConvergencePanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            logConvergenceCheck = new JCheckBox("Log convergence");
            logConvergenceCheck.setSelected(hybrOptions.getLogConvergence());
            logConvergencePanel.add(logConvergenceCheck);
            add(logConvergencePanel);

            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
            JButton okButton = new JButton("OK");
            okButton.addActionListener(new OkButtonActionListener());
            JButton cancelButton = new JButton("Cancel");
            cancelButton.addActionListener(new CancelButtonActionListener());
            buttonPanel.add(okButton);
            buttonPanel.add(cancelButton);
            pane.add(buttonPanel);
            validate();
            setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
            setResizable(false);
            pack();
            setLocationRelativeTo(null);
            setVisible(true);
        }

        private class ThresholdFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double threshold = 0;
                try {
                    threshold = Double.parseDouble(thresholdField.getText());
                } catch (Exception ex) {
                    IJ.error("Threshold must be a nonnegative value.");
                }
                if (threshold < 0.0) {
                    IJ.error("Threshold must be a nonnegative value.");
                }
            }
        }

        private class ThresholdCheckItemListener implements ItemListener {
            public void itemStateChanged(ItemEvent e) {
                if (thresholdCheck.isSelected()) {
                    thresholdField.setEnabled(true);
                } else {
                    thresholdField.setEnabled(false);
                }
            }
        }

        private class OmegaFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double omega = 0;
                try {
                    omega = Double.parseDouble(omegaField.getText());
                } catch (Exception ex) {
                    IJ.error("Omega must be a nonnegative floating-point number.");
                }
                if (omega < 0.0) {
                    IJ.error("Omega must be a nonnegative floating-point number.");
                }
            }
        }

        private class FlatTolFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double flatGCVTol = 0;
                try {
                    flatGCVTol = Double.parseDouble(flatTolField.getText());
                } catch (Exception ex) {
                    IJ.error("GCV tolerance must be a positive floating-point number.");
                }
                if (flatGCVTol <= 0.0) {
                    IJ.error("Omega must be a positive floating-point number.");
                }
            }
        }

        private class RegMethodChoiceActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                if (regChoice.getSelectedIndex() == 1) {// WGCV
                    omegaField.setEnabled(true);
                } else {
                    omegaField.setEnabled(false);
                }
                if (regChoice.getSelectedIndex() == 3) {// None
                    regField.setEnabled(true);
                    flatTolField.setEnabled(false);
                } else {
                    regField.setEnabled(false);
                    flatTolField.setEnabled(true);
                }

            }
        }

        private class RegFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double regParam = 0;
                try {
                    regParam = Double.parseDouble(regField.getText());
                } catch (Exception ex) {
                    IJ.error("Regularization parameter must be a floating-point number between 0 and 1.");
                }
                if ((regParam < 0.0) || (regParam > 1.0)) {
                    IJ.error("Regularization parameter must be a floating-point number between 0 and 1.");
                }
            }
        }

        private class BegRegFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double begReg = 0;
                try {
                    begReg = Integer.parseInt(begRegField.getText());
                } catch (Exception ex) {
                    IJ.error("Iteration number must be an integer number greater than 1.");
                }
                if (begReg <= 1) {
                    IJ.error("Iteration number must be an integer number greater than 1.");
                }
            }
        }

        private class OkButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                if (assignOptions(hybrOptions) == true) {
                    hybrOptionsPanel.setVisible(false);
                    hybrOptionsSet = true;
                }
            }
        }

        private class CancelButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                hybrOptionsPanel.dispose();
                hybrOptionsSet = false;
            }
        }

        private boolean assignOptions(HyBROptions hybrOptions) {
            double regParam;
            try {
                regParam = Double.parseDouble(regField.getText());
            } catch (Exception ex) {
                IJ.error("Regularization parameter must be a floating-point number between 0 and 1.");
                return false;
            }
            if ((regParam < 0.0) || (regParam > 1.0)) {
                IJ.error("Regularization parameter must be a floating-point number between 0 and 1.");
                return false;
            }
            double flatTolerance;
            try {
                flatTolerance = Double.parseDouble(flatTolField.getText());
            } catch (Exception ex) {
                IJ.error("GCV tolerance must be a nonnegative floating-point number.");
                return false;
            }
            if (flatTolerance < 0.0) {
                IJ.error("GCV tolerance must be a nonnegative floating-point number.");
                return false;
            }
            double omega;
            try {
                omega = Double.parseDouble(omegaField.getText());
            } catch (Exception ex) {
                IJ.error("Omega must be a nonnegative floating-point number.");
                return false;
            }
            if (omega < 0.0) {
                IJ.error("Omega must be a nonnegative floating-point number.");
                return false;
            }
            int beginReg;
            try {
                beginReg = Integer.parseInt(begRegField.getText());
            } catch (Exception ex) {
                IJ.error("Iteration number must be an integer number greater than 1.");
                return false;
            }
            if (beginReg <= 1) {
                IJ.error("Iteration number must be an integer number greater than 1.");
                return false;
            }
            double threshold = 0;
            if (thresholdCheck.isSelected()) {
                try {
                    threshold = Double.parseDouble(thresholdField.getText());
                } catch (Exception ex) {
                    IJ.error("Threshold must be a nonnegative value.");
                    return false;
                }
                if (threshold < 0.0) {
                    IJ.error("Threshold must be a nonnegative value.");
                    return true;
                }
            }
            hybrOptions.setInnerSolver(HyBRInnerSolver.values()[solverChoice.getSelectedIndex()]);
            hybrOptions.setRegMethod(HyBRRegularizationMethod.values()[regChoice.getSelectedIndex()]);
            hybrOptions.setRegParameter(regParam);
            hybrOptions.setOmega(omega);
            hybrOptions.setReorthogonalize(reorthCheck.isSelected());
            hybrOptions.setBeginReg(beginReg);
            hybrOptions.setFlatTolerance(flatTolerance);
            hybrOptions.setLogConvergence(logConvergenceCheck.isSelected());
            hybrOptions.setUseThreashold(thresholdCheck.isSelected());
            hybrOptions.setThreashold(threshold);
            return true;
        }
    }

    private class WPLOptionsPanel extends JFrame {

        private static final long serialVersionUID = 1L;

        private JTextField gammaField, filterXYfield, changeThreshPercentField, thresholdField;

        private JCheckBox normalizeCheck, logMeanCheck, antiRingCheck, dbCheck, detectDivergenceCheck, thresholdCheck;

        public WPLOptionsPanel(String name) {
            super(name);
            init();
        }

        private void init() {
            Container pane = getContentPane();
            pane.setLayout(new BoxLayout(pane, BoxLayout.Y_AXIS));

            JPanel normalizePanel = new JPanel();
            normalizePanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            normalizeCheck = new JCheckBox("Narmalize PSF");
            normalizeCheck.setSelected(wplOptions.isNormalize());
            normalizePanel.add(normalizeCheck);
            pane.add(normalizePanel);

            JPanel logMeanPanel = new JPanel();
            logMeanPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            logMeanCheck = new JCheckBox("Log mean pixel value to track convergence");
            logMeanCheck.setSelected(wplOptions.isLogConvergence());
            logMeanPanel.add(logMeanCheck);
            pane.add(logMeanPanel);

            JPanel antiRingPanel = new JPanel();
            antiRingPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            antiRingCheck = new JCheckBox("Perform anti-ringing step");
            antiRingCheck.setSelected(wplOptions.isAntiRing());
            antiRingPanel.add(antiRingCheck);
            pane.add(antiRingPanel);

            JPanel detectDivergencePanel = new JPanel();
            detectDivergencePanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            detectDivergenceCheck = new JCheckBox("Detect divergence");
            detectDivergenceCheck.setSelected(wplOptions.isDetectDivergence());
            detectDivergencePanel.add(detectDivergenceCheck);
            pane.add(detectDivergencePanel);

            JPanel dbPanel = new JPanel();
            dbPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            dbCheck = new JCheckBox("Data (image, psf and result) in dB");
            dbCheck.setSelected(wplOptions.isDB());
            dbPanel.add(dbCheck);
            pane.add(dbPanel);

            JPanel thresholdPanel = new JPanel();
            thresholdPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            thresholdCheck = new JCheckBox("Threshold:  ");
            thresholdCheck.setSelected(wplOptions.isUseThreshold());
            thresholdCheck.addItemListener(new ThresholdCheckItemListener());
            thresholdPanel.add(thresholdCheck);
            thresholdField = new JTextField(Double.toString(wplOptions.getThreshold()), 6);
            thresholdField.setEnabled(wplOptions.isUseThreshold());
            thresholdField.addActionListener(new ThresholdFieldActionListener());
            thresholdPanel.add(thresholdField);
            pane.add(thresholdPanel);

            JPanel gammaPanel = new JPanel();
            gammaPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel gammaLabel = new JLabel("Wiener filter gamma (suggest 0 [<.0001] to turn off, 0.0001-0.1 as tests)");
            gammaPanel.add(gammaLabel);
            gammaField = new JTextField(Double.toString(wplOptions.getGamma()), 6);
            gammaField.addActionListener(new GammaFieldActionListener());
            gammaPanel.add(gammaField);
            pane.add(gammaPanel);

            JPanel filterXYPanel = new JPanel();
            filterXYPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel filterXYlabel = new JLabel("Low pass filter x and y, pixels (suggest 1, 0 to turn off)");
            filterXYlabel.setPreferredSize(gammaLabel.getPreferredSize());
            filterXYPanel.add(filterXYlabel);
            filterXYfield = new JTextField(Double.toString(wplOptions.getFilterXY()), 6);
            filterXYfield.addActionListener(new FilterXYFieldActionListener());
            filterXYPanel.add(filterXYfield);
            pane.add(filterXYPanel);

            JPanel changeThreshPercentPanel = new JPanel();
            changeThreshPercentPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel changeThreshPercentLabel = new JLabel("Terminate iteration if mean delta < x% (suggest 0.01, 0 to turn off)");
            changeThreshPercentLabel.setPreferredSize(gammaLabel.getPreferredSize());
            changeThreshPercentPanel.add(changeThreshPercentLabel);
            changeThreshPercentField = new JTextField(Double.toString(wplOptions.getChangeThreshPercent()), 6);
            changeThreshPercentField.addActionListener(new ChangeThreshPercentFieldActionListener());
            changeThreshPercentPanel.add(changeThreshPercentField);
            pane.add(changeThreshPercentPanel);

            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
            JButton okButton = new JButton("OK");
            okButton.addActionListener(new OkButtonActionListener());
            JButton cancelButton = new JButton("Cancel");
            cancelButton.addActionListener(new CancelButtonActionListener());
            buttonPanel.add(okButton);
            buttonPanel.add(cancelButton);
            pane.add(buttonPanel);
            validate();
            setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
            setResizable(false);
            pack();
            setLocationRelativeTo(null);
            setVisible(true);
        }

        private class ThresholdFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double threshold = 0;
                try {
                    threshold = Double.parseDouble(thresholdField.getText());
                } catch (Exception ex) {
                    IJ.error("Threshold must be a nonnegative value.");
                }
                if (threshold < 0.0) {
                    IJ.error("Threshold must be a nonnegative value.");
                }
            }
        }

        private class GammaFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double gamma = 0;
                try {
                    gamma = Double.parseDouble(gammaField.getText());
                } catch (Exception ex) {
                    IJ.error("Gamma must be a nonnegative value.");
                }
                if (gamma < 0.0) {
                    IJ.error("Gamma must be a nonnegative value.");
                }
            }
        }

        private class FilterXYFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double filterXY = 0;
                try {
                    filterXY = Double.parseDouble(filterXYfield.getText());
                } catch (Exception ex) {
                    IJ.error("Filter x and y must be a nonnegative value.");
                }
                if (filterXY < 0.0) {
                    IJ.error("Filter x and y must be a nonnegative value.");
                }
            }
        }

        private class ChangeThreshPercentFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                double changeThreshPercent = 0;
                try {
                    changeThreshPercent = Double.parseDouble(changeThreshPercentField.getText());
                } catch (Exception ex) {
                    IJ.error("Mean delta must be a nonnegative value.");
                }
                if (changeThreshPercent < 0.0) {
                    IJ.error("Mean delta must be a nonnegative value.");
                }
            }
        }

        private class ThresholdCheckItemListener implements ItemListener {
            public void itemStateChanged(ItemEvent e) {
                if (thresholdCheck.isSelected()) {
                    thresholdField.setEnabled(true);
                } else {
                    thresholdField.setEnabled(false);
                }
            }
        }

        private class OkButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                if (assignOptions(wplOptions) == true) {
                    wplOptionsFrame.setVisible(false);
                    wplOptionsSet = true;
                }
            }
        }

        private class CancelButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                wplOptionsFrame.dispose();
                wplOptionsSet = false;
            }
        }

        private boolean assignOptions(WPLOptions wplOptions) {
            double gamma = 0;
            try {
                gamma = Double.parseDouble(gammaField.getText());
            } catch (Exception ex) {
                IJ.error("Gamma must be a nonnegative value.");
                return false;
            }
            if (gamma < 0.0) {
                IJ.error("Gamma must be a nonnegative value.");
                return false;
            }

            double filterXY = 0;
            try {
                filterXY = Double.parseDouble(filterXYfield.getText());
            } catch (Exception ex) {
                IJ.error("Filter x and y must be a nonnegative value.");
                return false;
            }
            if (filterXY < 0.0) {
                IJ.error("Filter x and y must be a nonnegative value.");
                return false;
            }

            double changeThreshPercent = 0;
            try {
                changeThreshPercent = Double.parseDouble(changeThreshPercentField.getText());
            } catch (Exception ex) {
                IJ.error("Mean delta must be a nonnegative value.");
                return false;
            }
            if (changeThreshPercent < 0.0) {
                IJ.error("Mean delta must be a nonnegative value.");
                return false;
            }

            double threshold = 0;
            if (thresholdCheck.isSelected()) {
                try {
                    threshold = Double.parseDouble(thresholdField.getText());
                } catch (Exception ex) {
                    IJ.error("Threshold must be a nonnegative value.");
                    return false;
                }
                if (threshold < 0.0) {
                    IJ.error("Threshold must be a nonnegative value.");
                    return false;
                }
            }

            wplOptions.setAntiRing(antiRingCheck.isSelected());
            wplOptions.setChangeThreshPercent(changeThreshPercent);
            wplOptions.setDB(dbCheck.isSelected());
            wplOptions.setDetectDivergence(detectDivergenceCheck.isSelected());
            wplOptions.setFilterXY(filterXY);
            wplOptions.setGamma(gamma);
            wplOptions.setLogConvergence(logMeanCheck.isSelected());
            wplOptions.setNormalize(normalizeCheck.isSelected());
            wplOptions.setThreshold(threshold);
            wplOptions.setUseThreshold(thresholdCheck.isSelected());
            return true;
        }
    }

    private class PSFCreatePanel extends JFrame {

        private static final long serialVersionUID = -4018857998511946590L;

        public JTextField rowsField, columnsField;

        public PSFCreatePanel(String name) {
            super(name);
            init();
        }

        private void init() {
            Container pane = getContentPane();
            pane.setLayout(new BoxLayout(pane, BoxLayout.Y_AXIS));
            JPanel labelPanel = new JPanel();
            labelPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel label = new JLabel("Enter the number of PSFs in the form of 2D matrix (rows x columns)");
            labelPanel.add(label);
            pane.add(labelPanel);
            JPanel textPanel = new JPanel();
            textPanel.setLayout(new FlowLayout(FlowLayout.CENTER));
            rowsField = new JTextField("2", 3);
            JLabel times = new JLabel(" x ");
            columnsField = new JTextField("2", 3);
            textPanel.add(rowsField);
            textPanel.add(times);
            textPanel.add(columnsField);
            pane.add(textPanel);
            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
            JButton okButton = new JButton("OK");
            okButton.addActionListener(new OkButtonActionListener());
            JButton cancelButton = new JButton("Cancel");
            cancelButton.addActionListener(new CancelButtonActionListener());
            buttonPanel.add(okButton);
            buttonPanel.add(cancelButton);
            pane.add(buttonPanel);
            validate();
            setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
            setResizable(false);
            pack();
            setLocationRelativeTo(null);
            setVisible(true);
        }

        private class OkButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                if (checkTextFields() == true) {
                    psfCreatePanel.dispose();
                    psfEditPanel = new PSFEditPanel("Edit Spatially Variant PSF");
                }
            }
        }

        private class CancelButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                psfCreatePanel.dispose();
            }
        }

        private boolean checkTextFields() {
            try {
                psfRows = Integer.parseInt(rowsField.getText());
            } catch (Exception ex) {
                IJ.error("Number of rows must be a positive integer number.");
                return false;
            }
            if (psfRows < 1) {
                IJ.error("Number of rows must be a positive integer number.");
                return false;
            }
            try {
                psfColumns = Integer.parseInt(columnsField.getText());
            } catch (Exception ex) {
                IJ.error("Number of columns must be a positive integer number.");
                return false;
            }
            if (psfColumns < 1) {
                IJ.error("Number of columns must be a positive integer number.");
                return false;
            }
            return true;
        }
    }

    private class PSFEditPanel extends JFrame {

        private static final long serialVersionUID = 11865121724333240L;

        private final Dimension totalSize = new Dimension(350, 200);

        JButton[][] buttons;

        JButton okButton, cancelButton;

        int counter;

        public PSFEditPanel(String name) {
            super(name);
            init();
        }

        private void init() {
            Container pane = getContentPane();
            pane.setLayout(new BorderLayout());
            JPanel labelPanel = new JPanel();
            labelPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel label = new JLabel("Define all PSFs");
            labelPanel.add(label);
            pane.add(labelPanel, BorderLayout.PAGE_START);
            JPanel tablePanel = new JPanel();
            counter = psfRows * psfColumns;
            tablePanel.setLayout(new GridLayout(psfRows, psfColumns));
            buttons = new JButton[psfRows][psfColumns];
            for (int r = 0; r < psfRows; r++) {
                for (int c = 0; c < psfColumns; c++) {
                    buttons[r][c] = new JButton("Null");
                    buttons[r][c].setToolTipText("Null");
                    buttons[r][c].addActionListener(new PSFButtonsActionListener());
                    buttons[r][c].setHorizontalAlignment(SwingConstants.LEFT);
                    tablePanel.add(buttons[r][c]);
                }
            }
            pane.add(tablePanel, BorderLayout.CENTER);
            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
            okButton = new JButton("OK");
            okButton.setEnabled(false);
            okButton.addActionListener(new OkButtonActionListener());
            cancelButton = new JButton("Cancel");
            cancelButton.addActionListener(new CancelButtonActionListener());
            buttonPanel.add(okButton);
            buttonPanel.add(cancelButton);
            pane.add(buttonPanel, BorderLayout.PAGE_END);
            setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
            setResizable(true);
            setSize(totalSize);
            setLocationRelativeTo(null);
            setVisible(true);
        }

        private class PSFButtonsActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                OpenDialog od = new OpenDialog("Open file", "");
                String directory = od.getDirectory();
                String name = od.getFileName();
                if (name != null) {
                    String path = directory + name;
                    ((JButton) e.getSource()).setText(path);
                    ((JButton) e.getSource()).setToolTipText(path);
                    ((JButton) e.getSource()).setIcon(new ImageIcon(path));
                    ((JButton) e.getSource()).setAlignmentX(CENTER_ALIGNMENT);
                    ((JButton) e.getSource()).setAlignmentY(CENTER_ALIGNMENT);
                    counter--;
                    if (counter == 0) {
                        okButton.setEnabled(true);
                    }
                }
            }
        }

        private class OkButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                imPSF = new ImagePlus[psfRows][psfColumns];
                Opener o = new Opener();
                for (int r = 0; r < psfRows; r++) {
                    for (int c = 0; c < psfColumns; c++) {
                        imPSF[r][c] = o.openImage(buttons[r][c].getText());
                        ImageProcessor ipPSF = imPSF[r][c].getProcessor();
                        if (ipPSF instanceof ColorProcessor) {
                            IJ.showMessage("RGB images are not currently supported.");
                            return;
                        }
                        if (imPSF[r][c].getStackSize() > 1) {
                            IJ.showMessage("For 3D images use Parallel Iterative Deconvolution 3D");
                            return;
                        }
                    }
                }
                psfEditPanel.setVisible(false);
                editPSFButton.setEnabled(true);
            }
        }

        private class CancelButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                psfEditPanel.setVisible(false);
            }
        }
    }

    private class MainPanel extends JFrame {

        private static final long serialVersionUID = 3975356344081858245L;

        private final Cursor defaultCursor = new Cursor(Cursor.DEFAULT_CURSOR);

        private final Cursor waitCursor = new Cursor(Cursor.WAIT_CURSOR);

        private static final int width = 340;

        public MainPanel(String name) {
            super(name);
            ConcurrencyUtils.setUseJCublas(false);
            windowIDs = WindowManager.getIDList();
            if (windowIDs != null) {
                imageTitles = new String[windowIDs.length];
                for (int i = 0; i < windowIDs.length; i++) {
                    ImagePlus im = WindowManager.getImage(windowIDs[i]);
                    if (im != null)
                        imageTitles[i] = im.getTitle();
                    else
                        imageTitles[i] = "";
                }
            }
            init();
        }

        private void init() {
            ToolTipManager.sharedInstance().setDismissDelay(20000);
            Container pane = getContentPane();
            pane.setLayout(new BoxLayout(pane, BoxLayout.Y_AXIS));
            // --------------------------------------------------------------
            JPanel blurPanel = new JPanel();
            blurPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel blurLabel = new JLabel("Blurred image:");
            blurLabel.setPreferredSize(new Dimension(90, blurLabel.getPreferredSize().height));
            blurPanel.add(blurLabel);
            if (windowIDs != null) {
                blurChoice = new JComboBox(imageTitles);
                blurChoice.setSelectedIndex(0);
            } else {
                blurChoice = new JComboBox();
            }
            blurChoice.setPreferredSize(new Dimension(width, blurChoice.getPreferredSize().height));
            blurChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            blurChoice.addActionListener(new BlurChoiceActionListener());
            blurChoice.setToolTipText("<html>Choose a blurred image.</html>");
            blurPanel.add(blurChoice);
            // --------------------------------------------------------------
            JPanel psfPanel = new JPanel();
            psfPanel.setLayout(new GridLayout(2, 1));
            Border border = new TitledBorder(null, null, TitledBorder.LEFT, TitledBorder.TOP);
            psfPanel.setBorder(border);
            JPanel psfChoicePanel = new JPanel();
            psfChoicePanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel psfLabel = new JLabel("PSF:");
            psfLabel.setPreferredSize(new Dimension(85, psfLabel.getPreferredSize().height));
            psfChoicePanel.add(psfLabel);
            if (windowIDs != null) {
                psfChoice = new JComboBox(imageTitles);
            } else {
                psfChoice = new JComboBox();
            }
            psfChoice.setPreferredSize(new Dimension(width, psfChoice.getPreferredSize().height));
            if (windowIDs != null) {
                if (windowIDs.length > 1) {
                    psfChoice.setSelectedIndex(1);
                } else {
                    psfChoice.setSelectedIndex(0);
                }
            }
            psfChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            psfChoice.addActionListener(new PsfChoiceActionListener());
            psfChoice.setToolTipText("<html>Choose a PSF image.</html>");
            psfChoicePanel.add(psfChoice);
            psfPanel.add(psfChoicePanel);
            JPanel psfVariantPanel = new JPanel();
            psfVariantPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
            variantPSFCheck = new JCheckBox("Spatially variant PSF");
            variantPSFCheck.setSelected(false);
            variantPSFCheck.addItemListener(new PSFCheckItemListener());
            psfVariantPanel.add(variantPSFCheck);
            definePSFButton = new JButton("Define");
            definePSFButton.addActionListener(new DefinePSFButtonActionListener());
            definePSFButton.setEnabled(false);
            psfVariantPanel.add(definePSFButton);
            editPSFButton = new JButton("Edit");
            editPSFButton.setEnabled(false);
            editPSFButton.addActionListener(new EditPSFButtonActionListener());
            psfVariantPanel.add(editPSFButton);
            psfPanel.add(psfVariantPanel);
            // --------------------------------------------------------------
            JPanel methodPanel = new JPanel();
            methodPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel methodLabel = new JLabel("Method:");
            methodLabel.setPreferredSize(new Dimension(90, methodLabel.getPreferredSize().height));
            methodPanel.add(methodLabel);
            methodChoice = new JComboBox(methodNames);
            methodChoice
                    .setToolTipText("<html>Choose a method:<br><ul><li>MRNSD - Modified Residual Norm Steepest Descent.<br>This is a nonnegatively constrained algorithm.</li><li>WPL - Wiener Filter Preconditioned Landweber.<br>This is a nonnegatively constrained algorithm.</li><li>CGLS - Conjugate Gradient for Least Squares.</li><li>HyBR - Hybrid Bidiagonalization Regularization.</li></ul></html>");
            methodChoice.setSelectedIndex(0);
            methodChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            methodChoice.addActionListener(new MethodChoiceActionListener());
            methodPanel.add(methodChoice);
            optionsButton = new JButton("Options");
            optionsButton.addActionListener(new OptionsButtonActionListener());
            methodPanel.add(optionsButton);
            // --------------------------------------------------------------
            JPanel precondPanel = new JPanel();
            precondPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel precondLabel1 = new JLabel("Preconditioner:");
            precondLabel1.setPreferredSize(new Dimension(90, precondLabel1.getPreferredSize().height));
            precondPanel.add(precondLabel1);
            precondChoice = new JComboBox(precondNames);
            precondChoice.setSelectedIndex(0);
            precondChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            precondChoice.addActionListener(new PrecondChoiceActionListener());
            precondChoice.setToolTipText("<html>Choose a preconditioner:<br><ul><li>FFT preconditioner - based on the Fast Fourier Transform.</li><li>None - no preconditioner.</li></ul></html>");
            precondPanel.add(precondChoice);
            JLabel precondLabel2 = new JLabel("Tolerance:");
            precondPanel.add(precondLabel2);
            precondField = new JTextField("0.0", 5);
            precondField.addActionListener(new PrecondFieldActionListener());
            precondField.setToolTipText("<html>A tolerance to \"regularize\" the preconditioner.</html>");
            precondField.setEnabled(false);
            precondPanel.add(precondField);
            precondCheck = new JCheckBox("Auto");
            precondCheck.setToolTipText("<html>Automatic choice of a tolerance for preconditioner<br>(based on the Generalized Cross Validation).</html>");
            precondCheck.setSelected(true);
            precondCheck.addItemListener(new PrecondCheckItemListener());
            precondPanel.add(precondCheck);

            // --------------------------------------------------------------
            JPanel boundaryPanel = new JPanel();
            boundaryPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel boundaryLabel = new JLabel("Boundary:");
            boundaryLabel.setPreferredSize(new Dimension(90, boundaryLabel.getPreferredSize().height));
            boundaryPanel.add(boundaryLabel);
            boundaryChoice = new JComboBox(boundaryNames);
            boundaryChoice.setSelectedIndex(0);
            boundaryChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            boundaryChoice.setToolTipText("<html>Choose boundary conditions.</html>");
            boundaryPanel.add(boundaryChoice);
            // --------------------------------------------------------------
            JPanel resizingPanel = new JPanel();
            resizingPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel resizingLabel = new JLabel("Resizing:");
            resizingLabel.setPreferredSize(new Dimension(90, resizingLabel.getPreferredSize().height));
            resizingPanel.add(resizingLabel);
            resizingChoice = new JComboBox(resizingNames);
            resizingChoice.setSelectedIndex(0);
            resizingChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            resizingChoice.setToolTipText("<html>Choose resizing.</html>");
            resizingPanel.add(resizingChoice);
            // --------------------------------------------------------------
            JPanel outputPanel = new JPanel();
            outputPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel outputLabel = new JLabel("Output:");
            outputLabel.setPreferredSize(new Dimension(90, outputLabel.getPreferredSize().height));
            outputPanel.add(outputLabel);
            outputChoice = new JComboBox(outputNames);
            outputChoice.setSelectedIndex(0);
            outputChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            outputChoice.setToolTipText("<html>Choose a type of deblurred image.</html>");
            outputPanel.add(outputChoice);
            // --------------------------------------------------------------
            JPanel precisionPanel = new JPanel();
            precisionPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel precisionLabel = new JLabel("PrecisionType:");
            precisionLabel.setPreferredSize(new Dimension(90, precisionLabel.getPreferredSize().height));
            precisionPanel.add(precisionLabel);
            precisionChoice = new JComboBox(precisionNames);
            precisionChoice.setSelectedIndex(0);
            precisionChoice.setAlignmentX(Component.LEFT_ALIGNMENT);
            precisionChoice.setToolTipText("<html>Choose precision.</html>");
            precisionPanel.add(precisionChoice);

            // --------------------------------------------------------------
            JPanel itersPanel = new JPanel();
            itersPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel itersLabel = new JLabel("Max number of iterations:");
            itersPanel.add(itersLabel);
            itersField = new JTextField("5", 3);
            itersField.addActionListener(new ItersFieldActionListener());
            itersField.setToolTipText("<html>The maximum number of iterations.</html>");
            itersPanel.add(itersField);
            itersCheck = new JCheckBox("Show iterations");
            itersCheck.setToolTipText("<html>Show restored image after each iteration.</html>");
            itersCheck.setSelected(false);
            itersPanel.add(itersCheck);

            // --------------------------------------------------------------
            JPanel threadsPanel = new JPanel();
            threadsPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
            JLabel threadsLabel = new JLabel("Max number of threads (power of 2):  ");
            threadsPanel.add(threadsLabel);
            ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.getNumberOfThreads());
            threadsField = new JTextField(Integer.toString(ConcurrencyUtils.getNumberOfThreads()), 3);
            threadsField.addActionListener(new ThreadsFieldActionListener());
            threadsPanel.add(threadsField);
            // --------------------------------------------------------------
            JPanel buttonPanel = new JPanel();
            buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
            deconvolveButton = new JButton("Deconvolve");
            deconvolveButton.addActionListener(new DeconvolveButtonActionListener());
            if (windowIDs == null) {
                deconvolveButton.setEnabled(false);
            }
            buttonPanel.add(deconvolveButton);
            cancelButton = new JButton("Cancel");
            cancelButton.addActionListener(new CancelButtonActionListener());
            buttonPanel.add(cancelButton);
            // --------------------------------------------------------------
            pane.add(blurPanel);
            pane.add(psfPanel);
            pane.add(methodPanel);
            pane.add(precondPanel);
            pane.add(boundaryPanel);
            pane.add(resizingPanel);
            pane.add(outputPanel);
            pane.add(precisionPanel);
            pane.add(itersPanel);
            pane.add(threadsPanel);
            pane.add(buttonPanel);
            validate();
            setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
            setResizable(false);
            pack();
            setLocationRelativeTo(null);
            setVisible(true);
        }

        private class BlurChoiceActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                windowIDs = WindowManager.getIDList();
                if (windowIDs != null) {
                    deconvolveButton.setEnabled(true);
                } else {
                    deconvolveButton.setEnabled(false);
                }
            }
        }

        private class PsfChoiceActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                windowIDs = WindowManager.getIDList();
                if (windowIDs != null) {
                    deconvolveButton.setEnabled(true);
                } else {
                    deconvolveButton.setEnabled(false);
                }
            }
        }

        private class MethodChoiceActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                MethodType selMethod = MethodType.values()[methodChoice.getSelectedIndex()];
                if (selMethod == MethodType.WPL) {
                    precondChoice.setEnabled(false);
                    precondCheck.setEnabled(false);
                    precondField.setEnabled(false);
                    variantPSFCheck.setSelected(false);
                    variantPSFCheck.setEnabled(false);
                    definePSFButton.setEnabled(false);
                    editPSFButton.setEnabled(false);
                } else {
                    precondChoice.setEnabled(true);
                    precondCheck.setEnabled(true);
                    if (!precondCheck.isSelected()) {
                        precondField.setEnabled(true);
                    }
                    variantPSFCheck.setEnabled(true);

                }
            }
        }

        private class PrecondChoiceActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                PreconditionerType selPrecond = PreconditionerType.values()[precondChoice.getSelectedIndex()];
                switch (selPrecond) {
                case FFT:
                    precondCheck.setEnabled(true);
                    if (precondCheck.isSelected() == false) {
                        precondField.setEnabled(true);
                    } else {
                        precondField.setEnabled(false);
                    }
                    break;
                case NONE:
                    precondField.setEnabled(false);
                    precondCheck.setEnabled(false);
                    break;
                }
            }
        }

        private class ItersFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                try {
                    maxIters = Integer.parseInt(itersField.getText());
                } catch (Exception ex) {
                    IJ.error("Number of iterations must be a positive integer.");
                }
                if (maxIters < 1) {
                    IJ.error("Number of iterations must be a positive integer.");
                }
            }
        }

        private class ThreadsFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                try {
                    threads = Integer.parseInt(threadsField.getText());
                } catch (Exception ex) {
                    IJ.error("Number of threads must be power of 2.");
                    return;
                }
                if (threads < 1) {
                    IJ.error("Number of threads must be power of 2.");
                    return;
                }
                if (!ConcurrencyUtils.isPowerOf2(threads)) {
                    IJ.error("Number of threads must be power of 2.");
                    return;
                }
                ConcurrencyUtils.setNumberOfThreads(threads);
            }
        }

        private class PSFCheckItemListener implements ItemListener {
            public void itemStateChanged(ItemEvent e) {
                if (variantPSFCheck.isSelected() == true) {
                    definePSFButton.setEnabled(true);
                    editPSFButton.setEnabled(false);
                    psfChoice.setEnabled(false);
                } else {
                    definePSFButton.setEnabled(false);
                    editPSFButton.setEnabled(false);
                    psfChoice.setEnabled(true);
                }
            }
        }

        private class PrecondCheckItemListener implements ItemListener {
            public void itemStateChanged(ItemEvent e) {
                if (precondCheck.isSelected() == true) {
                    precondField.setEnabled(false);
                } else {
                    precondField.setEnabled(true);
                }
            }
        }

        private class OptionsButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                final MethodType selMethod = MethodType.values()[methodChoice.getSelectedIndex()];
                switch (selMethod) {
                case CGLS:
                    if (cglsOptionsSet == true) {
                        cglsOptionsPanel.setVisible(true);
                    } else {
                        cglsOptionsPanel = new CGLSOptionsPanel("CGLS options");
                    }
                    break;
                case MRNSD:
                    if (mrnsdOptionsSet == true) {
                        mrnsdOptionsPanel.setVisible(true);
                    } else {
                        mrnsdOptionsPanel = new MRNSDOptionsPanel("MRNSD options");
                    }
                    break;
                case HyBR:
                    if (hybrOptionsSet == true) {
                        hybrOptionsPanel.setVisible(true);
                    } else {
                        hybrOptionsPanel = new HyBROptionsPanel("HyBR options");
                    }
                    break;
                case WPL:
                    if (wplOptionsSet == true) {
                        wplOptionsFrame.setVisible(true);
                    } else {
                        wplOptionsFrame = new WPLOptionsPanel("WPL options");
                    }
                    break;
                }

            }
        }

        private class PrecondFieldActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                if (precondCheck.isSelected() == true) {
                    precTol = -1;
                } else {
                    try {
                        precTol = Double.parseDouble(precondField.getText());
                    } catch (Exception ex) {
                        IJ.error("Tolerance for preconditioner must be a number between 0 and 1.");
                    }
                    if ((precTol < 0) || (precTol > 1)) {
                        IJ.error("Tolerance for preconditioner must be a number between 0 and 1.");
                    }
                }
            }
        }

        private class DefinePSFButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                psfCreatePanel = new PSFCreatePanel("Create Spatially Variant PSF");
            }
        }

        private class EditPSFButtonActionListener implements ActionListener {

            public void actionPerformed(ActionEvent e) {
                psfEditPanel.setVisible(true);
            }
        }

        private class DeconvolveButtonActionListener implements ActionListener {
            private final Timer timer = new Timer();

            public void actionPerformed(ActionEvent e) {
                Thread deconvolveThread = new Thread(new Runnable() {
                    public void run() {
                        imB = WindowManager.getImage((String) blurChoice.getSelectedItem());
                        if (imB == null) {
                            IJ.error("Image " + (String) blurChoice.getSelectedItem() + " was renamed. Please choose the blurred image again.");
                            windowIDs = WindowManager.getIDList();
                            if (windowIDs != null) {
                                imageTitles = new String[windowIDs.length];
                                for (int i = 0; i < windowIDs.length; i++) {
                                    ImagePlus im = WindowManager.getImage(windowIDs[i]);
                                    if (im != null)
                                        imageTitles[i] = im.getTitle();
                                    else
                                        imageTitles[i] = "";
                                }
                            }
                            blurChoice.removeAllItems();
                            for (int i = 0; i < imageTitles.length; i++) {
                                blurChoice.addItem(imageTitles[i]);
                            }
                            blurChoice.revalidate();
                            return;
                        }
                        ImageProcessor ipB = imB.getProcessor();
                        if (ipB instanceof ColorProcessor) {
                            IJ.showMessage("RGB images are not currently supported.");
                            return;
                        }
                        if (imB.getStackSize() > 1) {
                            IJ.showMessage("For 3D images use Parallel Iterative Deconvolution 3D");
                            return;
                        }
                        if (variantPSFCheck.isSelected() == false) {
                            imPSF = new ImagePlus[1][1];
                            imPSF[0][0] = WindowManager.getImage((String) psfChoice.getSelectedItem());
                            if (imPSF[0][0] == null) {
                                IJ.error("Image " + (String) psfChoice.getSelectedItem() + " was renamed. Please choose the blurred image again.");
                                windowIDs = WindowManager.getIDList();
                                if (windowIDs != null) {
                                    imageTitles = new String[windowIDs.length];
                                    for (int i = 0; i < windowIDs.length; i++) {
                                        ImagePlus im = WindowManager.getImage(windowIDs[i]);
                                        if (im != null)
                                            imageTitles[i] = im.getTitle();
                                        else
                                            imageTitles[i] = "";
                                    }
                                }
                                psfChoice.removeAllItems();
                                for (int i = 0; i < imageTitles.length; i++) {
                                    psfChoice.addItem(imageTitles[i]);
                                }
                                psfChoice.revalidate();
                                return;
                            }
                            ImageProcessor ipPSF = imPSF[0][0].getProcessor();
                            if (ipPSF instanceof ColorProcessor) {
                                IJ.showMessage("RGB images are not currently supported.");
                                return;
                            }
                            if (imPSF[0][0].getStackSize() > 1) {
                                IJ.showMessage("For 3D images use Parallel Iterative Deconvolution 3D.");
                                return;
                            }
                        }
                        if (!checkTextFields())
                            return;
                        setCursor(waitCursor);
                        deconvolveButton.setEnabled(false);
                        cancelButton.setEnabled(false);
                        final MethodType selMethod = MethodType.values()[methodChoice.getSelectedIndex()];
                        final PreconditionerType selPrecond = PreconditionerType.values()[precondChoice.getSelectedIndex()];
                        final BoundaryType selBoundary = BoundaryType.values()[boundaryChoice.getSelectedIndex()];
                        final ResizingType selResizing = ResizingType.values()[resizingChoice.getSelectedIndex()];
                        final OutputType selOutput = OutputType.values()[outputChoice.getSelectedIndex()];
                        final PrecisionType selPrecision = PrecisionType.values()[precisionChoice.getSelectedIndex()];
                        clean_old_data();
                        timer.reset().start();
                        switch (selPrecision) {
                        case DOUBLE:
                            switch (selMethod) {
                            case CGLS:
                                dcgls = new CGLSDoubleIterativeDeconvolver2D(imB, imPSF, selPrecond, precTol, selBoundary, selResizing, selOutput, maxIters, itersCheck.isSelected(), cglsOptions);
                                imX = dcgls.deconvolve();
                                if (selPrecond != PreconditionerType.NONE) {
                                    precondField.setText(String.format("%.4f", dcgls.getPreconditioner().getTolerance()));
                                }
                                dcgls = null;
                                timer.stop();
                                break;
                            case MRNSD:
                                dmrnsd = new MRNSDDoubleIterativeDeconvolver2D(imB, imPSF, selPrecond, precTol, selBoundary, selResizing, selOutput, maxIters, itersCheck.isSelected(), mrnsdOptions);
                                imX = dmrnsd.deconvolve();
                                if (selPrecond != PreconditionerType.NONE) {
                                    precondField.setText(String.format("%.4f", dmrnsd.getPreconditioner().getTolerance()));
                                }
                                dmrnsd = null;
                                timer.stop();
                                break;
                            case HyBR:
                                dhybr = new HyBRDoubleIterativeDeconvolver2D(imB, imPSF, selPrecond, precTol, selBoundary, selResizing, selOutput, maxIters, itersCheck.isSelected(), hybrOptions);
                                imX = dhybr.deconvolve();
                                if (selPrecond != PreconditionerType.NONE) {
                                    precondField.setText(String.format("%.4f", dhybr.getPreconditioner().getTolerance()));
                                }
                                dhybr = null;
                                timer.stop();
                                break;
                            case WPL:
                                dwpl = new WPLDoubleIterativeDeconvolver2D(imB, imPSF[0][0], selBoundary, selResizing, selOutput, maxIters, itersCheck.isSelected(), wplOptions);
                                imX = dwpl.deconvolve();
                                dwpl = null;
                                timer.stop();
                                break;
                            }
                            break;
                        case SINGLE:
                            switch (selMethod) {
                            case CGLS:
                                fcgls = new CGLSFloatIterativeDeconvolver2D(imB, imPSF, selPrecond, (float) precTol, selBoundary, selResizing, selOutput, maxIters, itersCheck.isSelected(), cglsOptions);
                                imX = fcgls.deconvolve();
                                if (selPrecond != PreconditionerType.NONE) {
                                    precondField.setText(String.format("%.4f", fcgls.getPreconditioner().getTolerance()));
                                }
                                fcgls = null;
                                timer.stop();
                                break;
                            case MRNSD:
                                fmrnsd = new MRNSDFloatIterativeDeconvolver2D(imB, imPSF, selPrecond, (float) precTol, selBoundary, selResizing, selOutput, maxIters, itersCheck.isSelected(), mrnsdOptions);
                                imX = fmrnsd.deconvolve();
                                if (selPrecond != PreconditionerType.NONE) {
                                    precondField.setText(String.format("%.4f", fmrnsd.getPreconditioner().getTolerance()));
                                }
                                fmrnsd = null;
                                timer.stop();
                                break;
                            case HyBR:
                                fhybr = new HyBRFloatIterativeDeconvolver2D(imB, imPSF, selPrecond, (float) precTol, selBoundary, selResizing, selOutput, maxIters, itersCheck.isSelected(), hybrOptions);
                                imX = fhybr.deconvolve();
                                if (selPrecond != PreconditionerType.NONE) {
                                    precondField.setText(String.format("%.4f", fhybr.getPreconditioner().getTolerance()));
                                }
                                fhybr = null;
                                timer.stop();
                                break;
                            case WPL:
                                fwpl = new WPLFloatIterativeDeconvolver2D(imB, imPSF[0][0], selBoundary, selResizing, selOutput, maxIters, itersCheck.isSelected(), wplOptions);
                                imX = fwpl.deconvolve();
                                fwpl = null;
                                timer.stop();
                                break;
                            }
                            break;
                        }
                        if (selPrecond != PreconditionerType.NONE) {
                            imX.setTitle(WindowManager.makeUniqueName(imB.getShortTitle() + "_deblurred_" + shortMethodNames[methodChoice.getSelectedIndex()] + "_" + shortPrecondNames[precondChoice.getSelectedIndex()] + precondField.getText() + "_"
                                    + shortBoundaryNames[boundaryChoice.getSelectedIndex()]));
                        } else {
                            imX.setTitle(WindowManager.makeUniqueName(imB.getShortTitle() + "_deblurred_" + shortMethodNames[methodChoice.getSelectedIndex()] + "_none_" + shortBoundaryNames[boundaryChoice.getSelectedIndex()]));
                        }
                        if (selMethod == MethodType.WPL) {
                            imX.setTitle(WindowManager.makeUniqueName(imB.getShortTitle() + "_deblurred_" + shortMethodNames[methodChoice.getSelectedIndex()] + "_" + shortBoundaryNames[boundaryChoice.getSelectedIndex()]));
                        }
                        if (itersCheck.isSelected() == false) {
                            imX.show();
                        } else {
                            blurChoice.removeItem("(deblurred)");
                            psfChoice.removeItem("(deblurred)");
                            blurChoice.addItem(imX.getTitle());
                            psfChoice.addItem(imX.getTitle());
                            blurChoice.revalidate();
                            psfChoice.revalidate();
                        }
                        IJ.showStatus(timer.toString());
                        setCursor(defaultCursor);
                        deconvolveButton.setEnabled(true);
                        cancelButton.setEnabled(true);
                    }

                });
                deconvolveThread.setUncaughtExceptionHandler(new DefaultExceptionHandler());
                deconvolveThread.start();
            }
        }

        private class DefaultExceptionHandler implements Thread.UncaughtExceptionHandler {

            public void uncaughtException(Thread t, Throwable e) {
                StringWriter sw = new StringWriter();
                PrintWriter pw = new PrintWriter(sw, true);
                e.printStackTrace(pw);
                pw.flush();
                sw.flush();
                IJ.log(sw.toString());
                mainPanel.dispose();
                ImagePlus.removeImageListener(getImageListener());
                clean_all();
            }

        }

        private class CancelButtonActionListener implements ActionListener {
            public void actionPerformed(ActionEvent e) {
                mainPanel.dispose();
                ImagePlus.removeImageListener(getImageListener());
                clean_all();
            }
        }

        private boolean checkTextFields() {
            try {
                maxIters = Integer.parseInt(itersField.getText());
            } catch (Exception ex) {
                IJ.error("Number of iterations must be a positive integer.");
                return false;
            }
            if (maxIters < 1) {
                IJ.error("Number of iterations must be a positive integer.");
                return false;
            }
            if (precondCheck.isSelected() == true) {
                precTol = -1;
            } else {
                try {
                    precTol = Double.parseDouble(precondField.getText());
                } catch (Exception ex) {
                    IJ.error("Tolerance for preconditioner must be a number between 0 and 1.");
                    return false;
                }
                if ((precTol < 0) || (precTol > 1)) {
                    IJ.error("Tolerance for preconditioner must be a number between 0 and 1.");
                    return false;
                }
            }
            try {
                threads = Integer.parseInt(threadsField.getText());
            } catch (Exception ex) {
                IJ.error("Number of threads must be power of 2.");
                return false;
            }
            if (threads < 1) {
                IJ.error("Number of threads must be power of 2.");
                return false;
            }
            if (!ConcurrencyUtils.isPowerOf2(threads)) {
                IJ.error("Number of threads must be power of 2.");
                return false;
            }
            ConcurrencyUtils.setNumberOfThreads(threads);
            return true;
        }
    }

    public static void main(String args[]) {

        new ImageJ();
        IJ.open("D:\\Research\\Images\\astronaut-blur.png");
        IJ.open("D:\\Research\\Images\\astronaut-psf.png");
        IJ.runPlugIn("edu.emory.mathcs.restoretools.iterative.ParallelIterativeDeconvolution2D", null);
    }

}
