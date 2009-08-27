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
import ij.io.Opener;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;

import cern.colt.matrix.tdouble.algo.solver.HyBRInnerSolver;
import cern.colt.matrix.tdouble.algo.solver.HyBRRegularizationMethod;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
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
 * Benchmark for Parallel Iterative Deconvolution 2D
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class Benchmark2D {

    private static final String path = "/home/pwendyk/images/iterative/";

    private static final String blur_image = "astronaut4096-blur.png";

    private static final String psf_image = "astronaut4096-psf.png";

    private static final BoundaryType boundary = BoundaryType.REFLEXIVE;

    private static final ResizingType resizing = ResizingType.AUTO;

    private static final OutputType output = OutputType.FLOAT;

    private static final int NITER = 10;

    private static final int MAXITER = 5;

    private static final double PREC_TOL = 0.001; //-1 for auto

    private static final String format = "%.2f";

    public static void benchmarkDoubleCGLS2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking DoubleCGLS2D using " + threads + " threads");
        CGLSOptions options = new CGLSOptions(false, 0, false, 0, false);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            CGLSDoubleIterativeDeconvolver2D cgls = new CGLSDoubleIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.NONE, 0, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = cgls.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            cgls = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time: " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("DoubleCGLS2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER);
    }

    public static void benchmarkDoublePCGLS2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking DoublePCGLS2D using " + threads + " threads");
        CGLSOptions options = new CGLSOptions(false, 0, false, 0, false);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            CGLSDoubleIterativeDeconvolver2D cgls = new CGLSDoubleIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.FFT, PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = cgls.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            cgls = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time (tol =  " + PREC_TOL + "): " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("DoublePCGLS2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER, PREC_TOL);
    }

    public static void benchmarkDoubleMRNSD2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking DoubleMRNSD2D using " + threads + " threads");
        MRNSDOptions options = new MRNSDOptions(false, 0, false, 0, false);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            MRNSDDoubleIterativeDeconvolver2D mrnsd = new MRNSDDoubleIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.NONE, PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = mrnsd.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            mrnsd = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time: " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("DoubleMRNSD2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER);
    }

    public static void benchmarkDoublePMRNSD2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking DoublePMRNSD2D using " + threads + " threads");
        MRNSDOptions options = new MRNSDOptions(false, 0, false, 0, false);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            MRNSDDoubleIterativeDeconvolver2D mrnsd = new MRNSDDoubleIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.FFT, PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = mrnsd.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            mrnsd = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time (tol =  " + PREC_TOL + "): " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("DoublePMRNSD2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER, PREC_TOL);
    }

    public static void benchmarkDoubleHyBR2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking DoubleHyBR2D using " + threads + " threads");
        HyBROptions options = new HyBROptions(HyBRInnerSolver.TIKHONOV, HyBRRegularizationMethod.ADAPTWGCV, 0, 0, false, 2, 0, false, false, 0);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            HyBRDoubleIterativeDeconvolver2D hybr = new HyBRDoubleIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.NONE, PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = hybr.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            hybr = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time: " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("DoubleHyBR2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER);
    }

    public static void benchmarkDoublePHyBR2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking preconditioned DoubleHyBR2D using " + threads + " threads");
        HyBROptions options = new HyBROptions(HyBRInnerSolver.TIKHONOV, HyBRRegularizationMethod.ADAPTWGCV, 0, 0, false, 2, 0, false, false, 0);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            HyBRDoubleIterativeDeconvolver2D hybr = new HyBRDoubleIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.FFT, PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = hybr.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            hybr = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time (tol =  " + PREC_TOL + "): " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("DoublePHyBR2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER, PREC_TOL);
    }

    public static void benchmarkDoubleWPL2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking DoubleWPL2D using " + threads + " threads");
        WPLOptions options = new WPLOptions(0, 1.0, 1.0, true, false, false, 0, false, false, false, 0);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            WPLDoubleIterativeDeconvolver2D wpl = new WPLDoubleIterativeDeconvolver2D(blurImage, psfImage[0][0], boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = wpl.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            wpl = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time: " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("DoubleWPL2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER);
    }

    public static void benchmarkDoublePWPL2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking DoublePWPL2D using " + threads + " threads");
        WPLOptions options = new WPLOptions(0.01, 1.0, 1.0, true, false, false, 0, false, false, false, 0);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            WPLDoubleIterativeDeconvolver2D wpl = new WPLDoubleIterativeDeconvolver2D(blurImage, psfImage[0][0], boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = wpl.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            wpl = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time (tol =  " + PREC_TOL + "): " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("DoublePWPL2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER, PREC_TOL);
    }

    public static void benchmarkFloatCGLS2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking FloatCGLS2D using " + threads + " threads");
        CGLSOptions options = new CGLSOptions(false, 0, false, 0, false);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            CGLSFloatIterativeDeconvolver2D cgls = new CGLSFloatIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.NONE, 0, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = cgls.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            cgls = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time: " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("FloatCGLS2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER);
    }

    public static void benchmarkFloatPCGLS2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking FloatPCGLS2D using " + threads + " threads");
        CGLSOptions options = new CGLSOptions(false, 0, false, 0, false);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            CGLSFloatIterativeDeconvolver2D cgls = new CGLSFloatIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.FFT, (float) PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = cgls.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            cgls = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time (tol =  " + PREC_TOL + "): " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("FloatPCGLS2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER, PREC_TOL);
    }

    public static void benchmarkFloatMRNSD2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking FloatMRNSD2D using " + threads + " threads");
        MRNSDOptions options = new MRNSDOptions(false, 0, false, 0, false);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            MRNSDFloatIterativeDeconvolver2D mrnsd = new MRNSDFloatIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.NONE, (float) PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = mrnsd.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            mrnsd = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time: " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("FloatMRNSD2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER);
    }

    public static void benchmarkFloatPMRNSD2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking FloatPMRNSD2D using " + threads + " threads");
        MRNSDOptions options = new MRNSDOptions(false, 0, false, 0, false);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            MRNSDFloatIterativeDeconvolver2D mrnsd = new MRNSDFloatIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.FFT, (float) PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = mrnsd.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            mrnsd = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time (tol =  " + PREC_TOL + "): " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("FloatPMRNSD2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER, PREC_TOL);
    }

    public static void benchmarkFloatHyBR2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking FloatHyBR2D using " + threads + " threads");
        HyBROptions options = new HyBROptions(HyBRInnerSolver.TIKHONOV, HyBRRegularizationMethod.ADAPTWGCV, 0, 0, false, 2, 0, false, false, 0);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            HyBRFloatIterativeDeconvolver2D hybr = new HyBRFloatIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.NONE, (float) PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = hybr.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            hybr = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time: " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("FloatHyBR2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER);
    }

    public static void benchmarkFloatPHyBR2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking preconditioned FloatHyBR2D using " + threads + " threads");
        HyBROptions options = new HyBROptions(HyBRInnerSolver.TIKHONOV, HyBRRegularizationMethod.ADAPTWGCV, 0, 0, false, 2, 0, false, false, 0);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            HyBRFloatIterativeDeconvolver2D hybr = new HyBRFloatIterativeDeconvolver2D(blurImage, psfImage, PreconditionerType.FFT, (float) PREC_TOL, boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = hybr.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            hybr = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time (tol =  " + PREC_TOL + "): " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("FloatPHyBR2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER, PREC_TOL);
    }

    public static void benchmarkFloatWPL2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking FloatWPL2D using " + threads + " threads");
        WPLOptions options = new WPLOptions(0, 1.0, 1.0, true, false, false, 0, false, false, false, 0);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            WPLFloatIterativeDeconvolver2D wpl = new WPLFloatIterativeDeconvolver2D(blurImage, psfImage[0][0], boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = wpl.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            wpl = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time: " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("FloatWPL2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER);
    }

    public static void benchmarkFloatPWPL2D(int threads) {
        ConcurrencyUtils.setNumberOfThreads(threads);
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + blur_image);
        ImagePlus[][] psfImage = new ImagePlus[1][1];
        psfImage[0][0] = o.openImage(path + psf_image);
        double av_time_deblur = 0;
        long elapsedTime_deblur = 0;
        System.out.println("Benchmarking FloatPWPL2D using " + threads + " threads");
        WPLOptions options = new WPLOptions(0.01, 1.0, 1.0, true, false, false, 0, false, false, false, 0);
        for (int i = 0; i < NITER; i++) {
            elapsedTime_deblur = System.nanoTime();
            WPLFloatIterativeDeconvolver2D wpl = new WPLFloatIterativeDeconvolver2D(blurImage, psfImage[0][0], boundary, resizing, output, MAXITER, false, options);
            ImagePlus imX = wpl.deconvolve();
            elapsedTime_deblur = System.nanoTime() - elapsedTime_deblur;
            av_time_deblur = av_time_deblur + elapsedTime_deblur;
            wpl = null;
            imX = null;
            
        }
        blurImage = null;
        psfImage = null;
        System.out.println("Average execution time (tol =  " + PREC_TOL + "): " + String.format(format, av_time_deblur / 1000000000.0 / (double) NITER) + " sec");
        writeResultsToFile("FloatPWPL2D_" + threads + "_threads.txt", (double) av_time_deblur / 1000000000.0 / (double) NITER, PREC_TOL);
    }

    public static void writeResultsToFile(String filename, double time_deblur) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(filename));
            out.write(new Date().toString());
            out.newLine();
            out.write("Number of processors: " + ConcurrencyUtils.getNumberOfThreads());
            out.newLine();
            out.write("deblur time: ");
            out.write(String.format(format, time_deblur));
            out.write(" seconds");
            out.newLine();
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void writeResultsToFile(String filename, double time_deblur, double prec_tol) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(filename));
            out.write(new Date().toString());
            out.newLine();
            out.write("Number of processors: " + ConcurrencyUtils.getNumberOfThreads());
            out.newLine();
            out.write("deblur time (tol = " + prec_tol + "): ");
            out.write(String.format(format, time_deblur));
            out.write(" seconds");
            out.newLine();
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        // benchmarkDoubleCGLS2D(1);
        // 
        // benchmarkDoubleCGLS2D(2);
        // 
        // benchmarkDoubleCGLS2D(4);
        // 
        // benchmarkDoubleCGLS2D(8);
        // 
        // benchmarkDoublePCGLS2D(1);
        // 
        // benchmarkDoublePCGLS2D(2);
        // 
        // benchmarkDoublePCGLS2D(4);
        // 
        // benchmarkDoublePCGLS2D(8);
        // 
        // benchmarkDoubleMRNSD2D(1);
        // 
        // benchmarkDoubleMRNSD2D(2);
        // 
        // benchmarkDoubleMRNSD2D(4);
        // 
        // benchmarkDoubleMRNSD2D(8);
        // 
        // benchmarkDoublePMRNSD2D(1);
        // 
        // benchmarkDoublePMRNSD2D(2);
        // 
        // benchmarkDoublePMRNSD2D(4);
        // 
        // benchmarkDoublePMRNSD2D(8);
        // 
        // benchmarkDoubleHyBR2D(1);
        // 
        // benchmarkDoubleHyBR2D(2);
        // 
        // benchmarkDoubleHyBR2D(4);
        // 
        // benchmarkDoubleHyBR2D(8);
        // 
        // benchmarkDoublePHyBR2D(1);
        // 
        // benchmarkDoublePHyBR2D(2);
        // 
        // benchmarkDoublePHyBR2D(4);
        // 
        // benchmarkDoublePHyBR2D(8);
        // 
        //
        // benchmarkFloatCGLS2D(1);
        // 
        // benchmarkFloatCGLS2D(2);
        // 
        // benchmarkFloatCGLS2D(4);
        // 
        // benchmarkFloatCGLS2D(8);
        // 
        benchmarkFloatPWPL2D(1);
        
        benchmarkFloatPWPL2D(2);
        
        benchmarkFloatPWPL2D(4);
        
        benchmarkFloatPWPL2D(8);
        
        benchmarkFloatPHyBR2D(1);
        
        benchmarkFloatPHyBR2D(2);
        
        benchmarkFloatPHyBR2D(4);
        
        benchmarkFloatPHyBR2D(8);
        
        benchmarkFloatPCGLS2D(1);
        
        benchmarkFloatPCGLS2D(2);
        
        benchmarkFloatPCGLS2D(4);
        
        benchmarkFloatPCGLS2D(8);
        
        // benchmarkFloatMRNSD2D(1);
        // 
        // benchmarkFloatMRNSD2D(2);
        // 
        // benchmarkFloatMRNSD2D(4);
        // 
        // benchmarkFloatMRNSD2D(8);
        // 
        benchmarkFloatPMRNSD2D(1);
        
        benchmarkFloatPMRNSD2D(2);
        
        benchmarkFloatPMRNSD2D(4);
        
        benchmarkFloatPMRNSD2D(8);
        //				
        //		benchmarkFloatHyBR2D(1);
        //		
        //		benchmarkFloatHyBR2D(2);
        //		
        //		benchmarkFloatHyBR2D(4);
        //		
        //		benchmarkFloatHyBR2D(8);
        System.exit(0);

    }
}
