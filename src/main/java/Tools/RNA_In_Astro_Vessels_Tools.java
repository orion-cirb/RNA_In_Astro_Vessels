package Tools;

import Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import Tools.Cellpose.CellposeTaskSettings;
import Tools.StardistOrion.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import static ij.IJ.setMinAndMax;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.*;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Voxel3D;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.Measure2Colocalisation;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.geom2.measurementsPopulation.PairObjects3DInt;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import net.haesleinhuepf.clijx.imagej2.ImageJ2Tubeness;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */



public class RNA_In_Astro_Vessels_Tools {
    
private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
public Calibration cal = new Calibration();
public float pixVol = 0;
public double minNucVol= 50;
public double maxNucVol = 1000; 
public double minCellVol= 100;
public double maxCellVol = Double.MAX_VALUE; 
public double minFociVol = 0.03;
public double maxFociVol = 10;

public double bg = 0;
public double stdBg = 0;

// Stardist
public Object syncObject = new Object();
public File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
protected String stardistFociModel = "fociRNA-1.2.zip";
public String stardistOutput = "Label Image"; 
public final double stardistPercentileBottom = 0.2;
public final double stardistPercentileTop = 99.8;
public final double stardistProbThresh = 0.2;
public final double stardistOverlayThresh = 0.35;



// Cellpose
private final String cellPoseEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\CellPose\\" : "/opt/miniconda3/envs/cellpose/";
private final String cellposeModelPath = IJ.isWindows()? System.getProperty("user.home")+"\\.cellpose\\models\\" : "";
public final String cellPoseAstroModel = cellposeModelPath+"cyto2_Iba1_microglia";
public final int cellPoseAstroDiameter = 120;
public final String cellPoseNucModel = "cyto2";
public final int cellPoseNucDiameter = 80;

      
private final CLIJ2 clij2 = CLIJ2.getInstance();

 
 /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    
    /**
     * Rescale to 16bits
     */
    public void to16bits(ImagePlus img) { 
        setMinAndMax(img,0, 65535);
        IJ.run(img, "16-bit", "");
    }
    
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
   
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
 /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs+1];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        channels[chs] = "None";
        return(channels);     
    }
    
    
    
    /** Vessel detection
     * Apply remove outliers
     * Threshold with Li
     * @param img
     * @param roi
     * @return 
     */
    public Object3DInt vesselsDetection(ImagePlus img, Roi roi) {
        ImagePlus imgOutliers = new Duplicator().run(img);
        IJ.run(imgOutliers, "Remove Outliers", "block_radius_x=10 block_radius_y=10 standard_deviations=1 stack");
        ImagePlus median_filter = median_filter(imgOutliers, 4, 4);
        flush_close(imgOutliers);
        threshold(median_filter, "Li", true);
        if (roi != null)
            clearOutSide(median_filter, roi);
        ImagePlus imgMin = min_filter(median_filter, 1, 1);
        flush_close(median_filter);
        ImagePlus imgMax = max_filter(imgMin, 1, 1);
        flush_close(imgMin);
        imgMax.setCalibration(cal);
        Object3DInt vesselObj = new Object3DInt(ImageHandler.wrap(imgMax), 255);
        flush_close(imgMax);
        return(vesselObj);
    }

    /**
     * Dialog ask for channels order and if needed spatial calibration
     * @param channels
     * @param showCal
     * @return ch
     */
    public String[] dialog(String[] channels) {
        String[] channelsName = {"DAPI", "Vessel", "Astro", "Foci"}; 
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
          
        gd.addMessage("Channels", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        int index = 0;
        for (String chName: channelsName) {
            gd.addChoice(chName + ": ", channels, channels[index]);
            index++;
        }
        gd.addMessage("Foci parameters", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        gd.addNumericField("Min vol (µm3): ", minFociVol, 3);
        gd.addNumericField("Max vol (µm3): ", maxFociVol, 3);
        // Calibration
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY pixel size : ", cal.pixelWidth, 3);
        gd.addNumericField("Z pixel size : ", cal.pixelDepth, 3);
        
        gd.showDialog();
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();
        if(gd.wasCanceled())
            ch = null;
        minFociVol = gd.getNextNumber();
        maxFociVol = gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();        
        pixVol = (float) (cal.pixelWidth*cal.pixelHeight*cal.pixelDepth);
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
    /**
     * Median filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus median_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       ImagePlus imgMed = clij2.pull(imgCLMed);
       clij2.release(imgCLMed);
       return(imgMed);
    } 
    
    /**
     * Max filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus max_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMax = clij2.create(imgCL);
       clij2.maximum3DBox(imgCL, imgCLMax, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       return(clij2.pull(imgCLMax));
    } 
    
     /**
     * Min filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus min_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMin = clij2.create(imgCL);
       clij2.minimum3DBox(imgCL, imgCLMin, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       return(clij2.pull(imgCLMin));
    } 
    
     /**
     * Clij2 Tubeness
     */
    private ImagePlus clij_Tubeness(ImagePlus img, float sigma) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLTube = clij2.create(imgCL);
        ImageJ2Tubeness ij2Tubeness = new ImageJ2Tubeness();
        ij2Tubeness.imageJ2Tubeness(clij2, imgCL, imgCLTube, sigma, 0f, 0f, 0f);
        ImagePlus imgTube = clij2.pull(imgCLTube);
        clij2.release(imgCL);
        clij2.release(imgCLTube);
        return(imgTube);
    }

    
    /**
     * Find Z with max intensity in stack
     * @param img
     * @return z
     */
    
    private int find_max(ImagePlus img) {
        double max = 0;
        int zmax = 0;
        for (int z = 1; z <= img.getNSlices(); z++) {
            ImageProcessor ip = img.getStack().getProcessor(z);
            ImageStatistics statistics = new ImageStatistics().getStatistics(ip, ImageStatistics.MEAN, img.getCalibration());
            double meanInt = statistics.mean;
            if (meanInt > max) {
                max = meanInt;
                zmax = z;
            }
        }
        return(zmax);
    }
    
    
    
    /**
     * Find astrocyte region
     * create 3D object from roi along Z axis
     * 
     */
    public Object3DInt findAstroRegion(ImagePlus img, Roi roi) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        ImagePlus imgCrop = new Duplicator().run(img);
        imgCrop.setRoi(poly);
        clearOutSide(imgCrop, poly);
        for (int z = 1; z <= imgCrop.getNSlices(); z++) {
            ImageProcessor ip = imgCrop.getImageStack().getProcessor(z);
            ip.setRoi(poly);
            ip.setColor(1);
            ip.fill(poly);
        }
        Object3DInt astroRegion = new Object3DInt(ImageHandler.wrap(imgCrop), 1);
        flush_close(imgCrop);
        return(astroRegion);        
    }
    
    
    /**
     * Threshold images and fill holes
     * @param img
     * @param thMed
     * @param fill 
     */
    public void threshold(ImagePlus img, String thMed, boolean fill) {
        //  Threshold and binarize
        img.setZ(find_max(img));
        img.updateAndDraw();
        IJ.setAutoThreshold(img, thMed + " dark");
        Prefs.blackBackground = false;
        IJ.run(img, "Convert to Mask", "method="+thMed+" background=Dark");
        if (fill) {
            IJ.run(img, "Fill Holes", "stack");
        }
    }
    
    
    public Objects3DIntPopulation getPopFromImage(ImagePlus img, Calibration cal) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        return pop;
    }
    
    
     /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(poly);
            ip.setBackgroundValue(0);
            ip.fillOutside(poly);
        }
        img.deleteRoi();
        img.updateAndDraw();
    }
    
    
    /**
     * Find min foci volume
     */
    public double findMinFociVolume(Objects3DIntPopulation fociPop) {
        DescriptiveStatistics dotVolStats = new DescriptiveStatistics();
        for (Object3DInt obj : fociPop.getObjects3DInt()) {
            dotVolStats.addValue(new MeasureVolume(obj).getVolumeUnit());
        }
        double minVol = dotVolStats.getMin();
        return(minVol);
    }
    
    
    /**
     * Measure object volume
     */
    public double measureObjVolume(Object3DInt obj) {
        double volObj = new MeasureVolume(obj).getVolumeUnit();
        return(volObj);
    }
    
     /**
     * Measure sum of pop volume
     */
    public double measurePopVolume(Objects3DIntPopulation pop) {
        double sumVol = 0;
        for(Object3DInt obj: pop.getObjects3DInt())
            sumVol += measureObjVolume(obj);
        return(sumVol);
    }  
    
    
     /**
     * Remove objects with size outside of given range
     */
    public void sizeFilterPop(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
    }
    
    
    /**
     * Apply StarDist 2D slice by slice
     * Label detections in 3D
     */
   public Objects3DIntPopulation stardistDetection(ImagePlus img, Roi roi) throws IOException{
       ImagePlus imgIn = new Duplicator().run(img);
       // StarDist
       File starDistModelFile = new File(modelsPath+File.separator+stardistFociModel);
       StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
       star.loadInput(imgIn);
       star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistProbThresh, stardistOverlayThresh, stardistOutput);
       star.run();
       
       // Label detections in 3D
       ImagePlus imgLabels = star.associateLabels();
       imgLabels.setCalibration(cal); 
       if (roi != null)
           clearOutSide(imgLabels, roi);
       
       // Get objects as a population of objects
       Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgLabels));  
       System.out.println(pop.getNbObjects()+" Stardist detections");
       
       // Filter objects
       sizeFilterPop(pop, minFociVol, maxFociVol);
       System.out.println(pop.getNbObjects()+ " detections remaining after size filtering");
       flush_close(img);
       flush_close(imgLabels);
       return(pop);
    }
    
     /**
     * Look for all 3D cells in a Z-stack: 
     * - apply CellPose in 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch threshold parameters
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus img, Roi roi, String model, int diam, double minVol, double maxVol){
        ImagePlus imgDup = new Duplicator().run(img);
        CellposeTaskSettings settings = new CellposeTaskSettings(model, 1, diam, cellPoseEnvDirPath);
        settings.setStitchThreshold(0.5); 
        settings.useGpu(true);
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgDup);
        ImagePlus imgOut = cellpose.run();
        if (roi != null)
            clearOutSide(imgOut, roi);
        imgOut.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        System.out.println(pop.getNbObjects() + " detections");
        sizeFilterPop(pop, minVol, maxVol);
        System.out.println(pop.getNbObjects() + " detections filtered out by size");
        flush_close(imgDup);
        flush_close(imgOut);
        return(pop);
    }
    
   
    
    /**
     * Find astrocyte nucleus 
     * Take nucleus than have max coloc with astrocyte
     * @param nucPop
     * @param astroObj
     * @return nucAstro
     */
    
    public Object3DInt findAstroNuc( Objects3DIntPopulation nucPop, Objects3DIntPopulation astroPop) {
        if (nucPop.getNbObjects() == 0 && astroPop.getNbObjects() == 0) 
            return(null);
        MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(nucPop, astroPop);
        List<PairObjects3DInt> allPairs = coloc.getAllPairs(true);
        return(allPairs.get(0).getObject3D2());
    }
    
   
    
    /**
     * Find foci population inside astro region
     * @param astroRegion
     * @param fociPop
     * @return fociInAstro
     */
    
    public Objects3DIntPopulation findFociInAstro(Object3DInt astroRegion, Objects3DIntPopulation fociPop) {
        Objects3DIntPopulation fociInAstro = new Objects3DIntPopulation();
        for (Object3DInt fociObj: fociPop.getObjects3DInt()) {
            Voxel3D centroid = new MeasureCentroid(fociObj).getCentroidAsVoxel();
            if (astroRegion.contains(new VoxelInt(centroid.getRoundX(), centroid.getRoundY(), centroid.getRoundZ(), 255)))
                fociInAstro.addObject(fociObj);
        }
        return(fociInAstro);
    }
    
    
    
    /**
     * Label object
     * @param obj
     * @param img 
     * @param fontSize 
     */
    public void labelObject(Object3DInt obj, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        
        BoundingBox bb = obj.getBoundingBox();
        int x = bb.xmin;
        int y = bb.ymin;
        img.setSlice(new MeasureCentroid(obj).getCentroidAsPoint().getRoundZ());
        ImageProcessor ip = img.getProcessor();
        ip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
        ip.setColor(255);
        ip.drawString(Integer.toString((int)obj.getLabel()), x, y);
        img.updateAndDraw();
    }
    
    
   /**
    * draw population
    * nucleus in blue
    * soma in green
    * foci in red
    * vessel in grey
    * @param somaObj
    * @param fociPop
    * @param vesselObj
    * @param img 
    * @param imgDir
    * @param imgName
    * @param roi
    */
    public void drawObjects(Object3DInt somaObj, Object3DInt vesselObj, Objects3DIntPopulation fociPop, ImagePlus img, String imgDir, String imgName, Roi roi) {  
        ImageHandler imgObjSoma = ImageHandler.wrap(img).createSameDimensions();
        somaObj.drawObject(imgObjSoma, 255); 
        ImageHandler imgObjFoci = ImageHandler.wrap(img).createSameDimensions();
        fociPop.drawInImage(imgObjFoci);
        ImageHandler imgObjVessel = ImageHandler.wrap(img).createSameDimensions();
        vesselObj.drawObject(imgObjVessel, 255);
       
        // save image for objects population
        ImagePlus[] imgColors = {imgObjFoci.getImagePlus(), imgObjSoma.getImagePlus(), null, imgObjVessel.getImagePlus()};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(img.getCalibration());
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(imgDir + imgName + "_Astro-"+roi.getName()+"-Objects.tif"); 
        flush_close(imgObjects);
    }
    
    
   
// Flush and close images
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    } 
         
    
    
}
