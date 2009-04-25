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

import ij.ImageJ;
import ij.ImagePlus;
import ij.io.Opener;
import cern.colt.Timer;
import edu.emory.mathcs.restoretools.Enums.OutputType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.BoundaryType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.PreconditionerType;
import edu.emory.mathcs.restoretools.iterative.IterativeEnums.ResizingType;
import edu.emory.mathcs.restoretools.iterative.cgls.CGLSDoubleIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.cgls.CGLSOptions;
import edu.emory.mathcs.restoretools.iterative.hybr.HyBRDoubleIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.hybr.HyBROptions;
import edu.emory.mathcs.restoretools.iterative.mrnsd.MRNSDDoubleIterativeDeconvolver2D;
import edu.emory.mathcs.restoretools.iterative.mrnsd.MRNSDDoubleIterativeDeconvolver3D;
import edu.emory.mathcs.restoretools.iterative.mrnsd.MRNSDOptions;

/**
 * Tests of iterative methods.
 * 
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * 
 */
public class TestsIterative {
    private static final String path = "d:\\Research\\Images\\";

    public static void testPMRNSD2D_variant() {
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + "star_cluster\\" + "star-cluster-blur.fit");
        blurImage.show();
        ImagePlus psfImage00 = o.openImage(path + "star_cluster\\" + "mpsf00.fit");
        ImagePlus psfImage01 = o.openImage(path + "star_cluster\\" + "mpsf01.fit");
        ImagePlus psfImage02 = o.openImage(path + "star_cluster\\" + "mpsf02.fit");
        ImagePlus psfImage03 = o.openImage(path + "star_cluster\\" + "mpsf03.fit");
        ImagePlus psfImage04 = o.openImage(path + "star_cluster\\" + "mpsf04.fit");
        ImagePlus psfImage05 = o.openImage(path + "star_cluster\\" + "mpsf05.fit");
        ImagePlus psfImage06 = o.openImage(path + "star_cluster\\" + "mpsf06.fit");
        ImagePlus psfImage07 = o.openImage(path + "star_cluster\\" + "mpsf07.fit");
        ImagePlus psfImage08 = o.openImage(path + "star_cluster\\" + "mpsf08.fit");
        ImagePlus psfImage09 = o.openImage(path + "star_cluster\\" + "mpsf09.fit");
        ImagePlus psfImage10 = o.openImage(path + "star_cluster\\" + "mpsf10.fit");
        ImagePlus psfImage11 = o.openImage(path + "star_cluster\\" + "mpsf11.fit");
        ImagePlus psfImage12 = o.openImage(path + "star_cluster\\" + "mpsf12.fit");
        ImagePlus psfImage13 = o.openImage(path + "star_cluster\\" + "mpsf13.fit");
        ImagePlus psfImage14 = o.openImage(path + "star_cluster\\" + "mpsf14.fit");
        ImagePlus psfImage15 = o.openImage(path + "star_cluster\\" + "mpsf15.fit");
        ImagePlus psfImage16 = o.openImage(path + "star_cluster\\" + "mpsf16.fit");
        ImagePlus psfImage17 = o.openImage(path + "star_cluster\\" + "mpsf17.fit");
        ImagePlus psfImage18 = o.openImage(path + "star_cluster\\" + "mpsf18.fit");
        ImagePlus psfImage19 = o.openImage(path + "star_cluster\\" + "mpsf19.fit");
        ImagePlus psfImage20 = o.openImage(path + "star_cluster\\" + "mpsf20.fit");
        ImagePlus psfImage21 = o.openImage(path + "star_cluster\\" + "mpsf21.fit");
        ImagePlus psfImage22 = o.openImage(path + "star_cluster\\" + "mpsf22.fit");
        ImagePlus psfImage23 = o.openImage(path + "star_cluster\\" + "mpsf23.fit");
        ImagePlus psfImage24 = o.openImage(path + "star_cluster\\" + "mpsf24.fit");

        ImagePlus[][] PSF = new ImagePlus[5][5];
        PSF[0][0] = psfImage04;
        PSF[0][1] = psfImage09;
        PSF[0][2] = psfImage14;
        PSF[0][3] = psfImage19;
        PSF[0][4] = psfImage24;

        PSF[1][0] = psfImage03;
        PSF[1][1] = psfImage08;
        PSF[1][2] = psfImage13;
        PSF[1][3] = psfImage18;
        PSF[1][4] = psfImage23;

        PSF[2][0] = psfImage02;
        PSF[2][1] = psfImage07;
        PSF[2][2] = psfImage12;
        PSF[2][3] = psfImage17;
        PSF[2][4] = psfImage22;

        PSF[3][0] = psfImage01;
        PSF[3][1] = psfImage06;
        PSF[3][2] = psfImage11;
        PSF[3][3] = psfImage16;
        PSF[3][4] = psfImage21;

        PSF[4][0] = psfImage00;
        PSF[4][1] = psfImage05;
        PSF[4][2] = psfImage10;
        PSF[4][3] = psfImage15;
        PSF[4][4] = psfImage20;

        MRNSDOptions options = new MRNSDOptions();
        Timer t = new Timer().start();
        MRNSDDoubleIterativeDeconvolver2D pmrnsd = new MRNSDDoubleIterativeDeconvolver2D(blurImage, PSF, PreconditionerType.FFT, -1, BoundaryType.REFLEXIVE, ResizingType.AUTO, OutputType.SAME_AS_SOURCE, 39, false, options);
        t.stop();
        System.out.println("Constructor: " + t);
        t.reset().start();
        ImagePlus imX = pmrnsd.deconvolve();
        t.stop();
        System.out.println("Deblur: " + t);
        System.out.println("niter=" + imX.getProperty("niter"));
        imX.show();
    }

    public static void testPCGLS2D_variant() {
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + "star_cluster\\" + "star-cluster-blur.fit");
        blurImage.show();
        ImagePlus psfImage00 = o.openImage(path + "star_cluster\\" + "mpsf00.fit");
        ImagePlus psfImage01 = o.openImage(path + "star_cluster\\" + "mpsf01.fit");
        ImagePlus psfImage02 = o.openImage(path + "star_cluster\\" + "mpsf02.fit");
        ImagePlus psfImage03 = o.openImage(path + "star_cluster\\" + "mpsf03.fit");
        ImagePlus psfImage04 = o.openImage(path + "star_cluster\\" + "mpsf04.fit");
        ImagePlus psfImage05 = o.openImage(path + "star_cluster\\" + "mpsf05.fit");
        ImagePlus psfImage06 = o.openImage(path + "star_cluster\\" + "mpsf06.fit");
        ImagePlus psfImage07 = o.openImage(path + "star_cluster\\" + "mpsf07.fit");
        ImagePlus psfImage08 = o.openImage(path + "star_cluster\\" + "mpsf08.fit");
        ImagePlus psfImage09 = o.openImage(path + "star_cluster\\" + "mpsf09.fit");
        ImagePlus psfImage10 = o.openImage(path + "star_cluster\\" + "mpsf10.fit");
        ImagePlus psfImage11 = o.openImage(path + "star_cluster\\" + "mpsf11.fit");
        ImagePlus psfImage12 = o.openImage(path + "star_cluster\\" + "mpsf12.fit");
        ImagePlus psfImage13 = o.openImage(path + "star_cluster\\" + "mpsf13.fit");
        ImagePlus psfImage14 = o.openImage(path + "star_cluster\\" + "mpsf14.fit");
        ImagePlus psfImage15 = o.openImage(path + "star_cluster\\" + "mpsf15.fit");
        ImagePlus psfImage16 = o.openImage(path + "star_cluster\\" + "mpsf16.fit");
        ImagePlus psfImage17 = o.openImage(path + "star_cluster\\" + "mpsf17.fit");
        ImagePlus psfImage18 = o.openImage(path + "star_cluster\\" + "mpsf18.fit");
        ImagePlus psfImage19 = o.openImage(path + "star_cluster\\" + "mpsf19.fit");
        ImagePlus psfImage20 = o.openImage(path + "star_cluster\\" + "mpsf20.fit");
        ImagePlus psfImage21 = o.openImage(path + "star_cluster\\" + "mpsf21.fit");
        ImagePlus psfImage22 = o.openImage(path + "star_cluster\\" + "mpsf22.fit");
        ImagePlus psfImage23 = o.openImage(path + "star_cluster\\" + "mpsf23.fit");
        ImagePlus psfImage24 = o.openImage(path + "star_cluster\\" + "mpsf24.fit");

        ImagePlus[][] PSF = new ImagePlus[5][5];
        PSF[0][0] = psfImage04;
        PSF[0][1] = psfImage09;
        PSF[0][2] = psfImage14;
        PSF[0][3] = psfImage19;
        PSF[0][4] = psfImage24;

        PSF[1][0] = psfImage03;
        PSF[1][1] = psfImage08;
        PSF[1][2] = psfImage13;
        PSF[1][3] = psfImage18;
        PSF[1][4] = psfImage23;

        PSF[2][0] = psfImage02;
        PSF[2][1] = psfImage07;
        PSF[2][2] = psfImage12;
        PSF[2][3] = psfImage17;
        PSF[2][4] = psfImage22;

        PSF[3][0] = psfImage01;
        PSF[3][1] = psfImage06;
        PSF[3][2] = psfImage11;
        PSF[3][3] = psfImage16;
        PSF[3][4] = psfImage21;

        PSF[4][0] = psfImage00;
        PSF[4][1] = psfImage05;
        PSF[4][2] = psfImage10;
        PSF[4][3] = psfImage15;
        PSF[4][4] = psfImage20;

        CGLSOptions options = new CGLSOptions();
        Timer t = new Timer().start();
        CGLSDoubleIterativeDeconvolver2D cgls = new CGLSDoubleIterativeDeconvolver2D(blurImage, PSF, PreconditionerType.FFT, -1, BoundaryType.REFLEXIVE, ResizingType.AUTO, OutputType.SAME_AS_SOURCE, 39, false, options);
        t.stop();
        System.out.println("Constructor: " + t);
        t.reset().start();
        ImagePlus imX = cgls.deconvolve();
        t.stop();
        System.out.println("Deblur: " + t);
        System.out.println("niter=" + imX.getProperty("niter"));
        imX.show();
    }

    public static void testPHyBR2D_variant() {
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + "star_cluster\\" + "star-cluster-blur.fit");
        blurImage.show();
        ImagePlus psfImage00 = o.openImage(path + "star_cluster\\" + "mpsf00.fit");
        ImagePlus psfImage01 = o.openImage(path + "star_cluster\\" + "mpsf01.fit");
        ImagePlus psfImage02 = o.openImage(path + "star_cluster\\" + "mpsf02.fit");
        ImagePlus psfImage03 = o.openImage(path + "star_cluster\\" + "mpsf03.fit");
        ImagePlus psfImage04 = o.openImage(path + "star_cluster\\" + "mpsf04.fit");
        ImagePlus psfImage05 = o.openImage(path + "star_cluster\\" + "mpsf05.fit");
        ImagePlus psfImage06 = o.openImage(path + "star_cluster\\" + "mpsf06.fit");
        ImagePlus psfImage07 = o.openImage(path + "star_cluster\\" + "mpsf07.fit");
        ImagePlus psfImage08 = o.openImage(path + "star_cluster\\" + "mpsf08.fit");
        ImagePlus psfImage09 = o.openImage(path + "star_cluster\\" + "mpsf09.fit");
        ImagePlus psfImage10 = o.openImage(path + "star_cluster\\" + "mpsf10.fit");
        ImagePlus psfImage11 = o.openImage(path + "star_cluster\\" + "mpsf11.fit");
        ImagePlus psfImage12 = o.openImage(path + "star_cluster\\" + "mpsf12.fit");
        ImagePlus psfImage13 = o.openImage(path + "star_cluster\\" + "mpsf13.fit");
        ImagePlus psfImage14 = o.openImage(path + "star_cluster\\" + "mpsf14.fit");
        ImagePlus psfImage15 = o.openImage(path + "star_cluster\\" + "mpsf15.fit");
        ImagePlus psfImage16 = o.openImage(path + "star_cluster\\" + "mpsf16.fit");
        ImagePlus psfImage17 = o.openImage(path + "star_cluster\\" + "mpsf17.fit");
        ImagePlus psfImage18 = o.openImage(path + "star_cluster\\" + "mpsf18.fit");
        ImagePlus psfImage19 = o.openImage(path + "star_cluster\\" + "mpsf19.fit");
        ImagePlus psfImage20 = o.openImage(path + "star_cluster\\" + "mpsf20.fit");
        ImagePlus psfImage21 = o.openImage(path + "star_cluster\\" + "mpsf21.fit");
        ImagePlus psfImage22 = o.openImage(path + "star_cluster\\" + "mpsf22.fit");
        ImagePlus psfImage23 = o.openImage(path + "star_cluster\\" + "mpsf23.fit");
        ImagePlus psfImage24 = o.openImage(path + "star_cluster\\" + "mpsf24.fit");

        ImagePlus[][] PSF = new ImagePlus[5][5];
        PSF[0][0] = psfImage04;
        PSF[0][1] = psfImage09;
        PSF[0][2] = psfImage14;
        PSF[0][3] = psfImage19;
        PSF[0][4] = psfImage24;

        PSF[1][0] = psfImage03;
        PSF[1][1] = psfImage08;
        PSF[1][2] = psfImage13;
        PSF[1][3] = psfImage18;
        PSF[1][4] = psfImage23;

        PSF[2][0] = psfImage02;
        PSF[2][1] = psfImage07;
        PSF[2][2] = psfImage12;
        PSF[2][3] = psfImage17;
        PSF[2][4] = psfImage22;

        PSF[3][0] = psfImage01;
        PSF[3][1] = psfImage06;
        PSF[3][2] = psfImage11;
        PSF[3][3] = psfImage16;
        PSF[3][4] = psfImage21;

        PSF[4][0] = psfImage00;
        PSF[4][1] = psfImage05;
        PSF[4][2] = psfImage10;
        PSF[4][3] = psfImage15;
        PSF[4][4] = psfImage20;

        HyBROptions options = new HyBROptions();
        options.setLogConvergence(true);
//        options.setFlatTolerance(0);
//        options.setUseThreashold(false);
        Timer t = new Timer().start();
        HyBRDoubleIterativeDeconvolver2D hybr = new HyBRDoubleIterativeDeconvolver2D(blurImage, PSF, PreconditionerType.FFT, -1, BoundaryType.REFLEXIVE, ResizingType.AUTO, OutputType.SAME_AS_SOURCE, 39, false, options);
        t.stop();
        System.out.println("Constructor: " + t);
        t.reset().start();
        ImagePlus imX = hybr.deconvolve();
        t.stop();
        System.out.println("Deblur: " + t);
        System.out.println("niter=" + imX.getProperty("niter"));
        imX.show();
    }

    public static void testPMRNSD3D_variant() {
        Opener o = new Opener();
        ImagePlus blurImage = o.openImage(path + "head\\head_blur.tif");
        blurImage.show();
        ImagePlus psfImage000 = o.openImage(path + "head\\head_psf.tif");
        ImagePlus psfImage001 = o.openImage(path + "head\\head_psf.tif");
        ImagePlus psfImage010 = o.openImage(path + "head\\head_psf.tif");
        ImagePlus psfImage011 = o.openImage(path + "head\\head_psf.tif");
        ImagePlus psfImage100 = o.openImage(path + "head\\head_psf.tif");
        ImagePlus psfImage101 = o.openImage(path + "head\\head_psf.tif");
        ImagePlus psfImage110 = o.openImage(path + "head\\head_psf.tif");
        ImagePlus psfImage111 = o.openImage(path + "head\\head_psf.tif");

        ImagePlus[][][] PSF = new ImagePlus[2][2][2];
        PSF[0][0][0] = psfImage000;
        PSF[0][0][1] = psfImage001;
        PSF[0][1][0] = psfImage010;
        PSF[0][1][1] = psfImage011;
        PSF[1][0][0] = psfImage100;
        PSF[1][0][1] = psfImage101;
        PSF[1][1][0] = psfImage110;
        PSF[1][0][0] = psfImage111;

        MRNSDOptions options = new MRNSDOptions();
        Timer t = new Timer().start();
        MRNSDDoubleIterativeDeconvolver3D pmrnsd = new MRNSDDoubleIterativeDeconvolver3D(blurImage, PSF, PreconditionerType.FFT, -1, BoundaryType.REFLEXIVE, ResizingType.AUTO, OutputType.SAME_AS_SOURCE, 39, false, options);
        t.stop();
        System.out.println("Constructor: " + t);
        t.reset().start();
        ImagePlus imX = pmrnsd.deconvolve();
        t.stop();
        System.out.println("Deblur: " + t);
        System.out.println("niter=" + imX.getProperty("niter"));
        imX.show();
    }

    public static void main(String[] args) {
        new ImageJ();
        testPHyBR2D_variant();
    }

}
