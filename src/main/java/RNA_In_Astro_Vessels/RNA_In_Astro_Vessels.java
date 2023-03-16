package RNA_In_Astro_Vessels;



import Tools.RNA_In_Astro_Vessels_Tools;
import ij.IJ;
import static ij.IJ.setMinAndMax;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.Region;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */



public class RNA_In_Astro_Vessels implements PlugIn {
    
    private String imageDir;
    private final boolean canceled = false;
    public  String outDirResults = "";
    private BufferedWriter outPutResults;
    
    
   
            
    Tools.RNA_In_Astro_Vessels_Tools tools = new RNA_In_Astro_Vessels_Tools();
            
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        
       
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose Directory Containing image Files...");
            if (imageDir == null) {
                return;
            }
            // Find images with file_ext extension ics, ics2
            String file_ext = "ics2";
            ArrayList<String> imageFiles = tools.findImages(imageDir, file_ext);
            if (imageFiles.size() == 0) {
                IJ.showStatus("Error", "No images found with " + file_ext + " extension");
                file_ext = "ics";
                IJ.showStatus("Trying to find images with " + file_ext + " extension");
                imageFiles = tools.findImages(imageDir, file_ext);
                if (imageFiles.size() == 0) {
                    IJ.showMessage("Error", "No images found with " + file_ext + " extension");
                    return;
                }
            }

            // create output folder
            outDirResults = imageDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write headers for results file
            FileWriter fileResults = new FileWriter(outDirResults + "Astro_results.xls", false);
            outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("ImageName\tRoi Name\tAstrocyte region volume (µm3)\tAstrocyte region foci number\tSoma volume (µm3)\tSoma astrocyte foci number\t"
                    + "Vessel volume (µm3)\tFoci in vessels\n");
            outPutResults.flush();
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.cal = tools.findImageCalib(meta);
            
             // Find channel names
            String[] chsName = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Channels dialog
            String[] channels = tools.dialog(chsName);
            if (channels == null) {
                IJ.showStatus("Plugin cancelled");
                return;
            } 
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                // Find ROI file
                String roi_file = imageDir+rootName+".zip";
                if (!new File(roi_file).exists()) {
                    roi_file = imageDir+rootName+".roi";
                    if (!new File(roi_file).exists()) {
                        IJ.showStatus("No ROI file found !");
                        return;
                    }
                }

                // find rois
                RoiManager rm = new RoiManager(false);
                rm.runCommand("Open", roi_file);
                // for each roi open image and crop
                for (Roi roi : rm.getRoisAsArray()) {
                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setId(f);
                    options.setSplitChannels(true);
                    options.setCrop(true);
                    Region reg = new Region(roi.getBounds().x, roi.getBounds().y, roi.getBounds().width, roi.getBounds().height);
                    options.setCropRegion(0, reg);
                    options.doCrop();
                    
                
                    /**
                    * Open channels
                    */
                    
                    // Nucleus channel
                    int indexCh = ArrayUtils.indexOf(chsName, channels[0]);
                    System.out.println("Opening Nucleus channel");
                    ImagePlus imgNuc = BF.openImagePlus(options)[indexCh];
                    // rescale to 16bits (Huygens images)
                    if (imgNuc.getBitDepth() == 32) {
                        setMinAndMax(imgNuc,0, 65535);
                        IJ.run(imgNuc, "16-bit", "");
                    }
                    
                    
                    // Find nucleus
                    Objects3DIntPopulation nucPop = tools.cellposeDetection(imgNuc, roi, tools.cellPoseNucModel, tools.cellPoseNucDiameter, 
                            tools.minNucVol, tools.maxNucVol);
                    System.out.println("--- Roi = "+roi.getName()+"---");
                    System.out.println("Nucleus number = "+nucPop.getNbObjects());                    
                        
                    // Astrocyte channel
                    // Find soma
                    IJ.showStatus("Opening Astrocyte channel");
                    indexCh = ArrayUtils.indexOf(chsName, channels[2]);
                    ImagePlus imgAstro = BF.openImagePlus(options)[indexCh];
                    if (imgAstro.getBitDepth() == 32)
                        tools.to16bits(imgAstro);
                    Objects3DIntPopulation astroPop = tools.cellposeDetection(imgAstro, roi, tools.cellPoseAstroModel, tools.cellPoseAstroDiameter, 
                        tools.minCellVol, tools.maxCellVol);
                        
                    if (astroPop.getNbObjects() != 0) { 
                        
                        // Find nucleus associated to astro cell
                        Object3DInt astroSoma = tools.findAstroNuc(nucPop, astroPop);
                        System.out.println("Soma astrocyte "+astroSoma.getLabel()+" found associated with nucleus");
                        
                        // Find Astrocyte region
                        Object3DInt astroRegion = tools.findAstroRegion(imgAstro, roi);
                        
                         // Vessels channel
                        IJ.showStatus("Opening vessels channel");
                        indexCh = ArrayUtils.indexOf(chsName, channels[1]);
                        ImagePlus imgVessels = BF.openImagePlus(options)[indexCh];
                        if (imgVessels.getBitDepth() == 32)
                            tools.to16bits(imgVessels);
                        Object3DInt vesselObj = tools.vesselsDetection(imgVessels, roi);
                        
                        // Dots channel
                        IJ.showStatus("Opening Dots channel");
                        indexCh = ArrayUtils.indexOf(chsName, channels[3]);
                        ImagePlus imgDots = BF.openImagePlus(options)[indexCh];
                        imgDots.setTitle(rootName+"_Dots");
                        if (imgDots.getBitDepth() == 32)
                            tools.to16bits(imgDots);
                         Objects3DIntPopulation fociPop = tools.stardistDetection(imgDots, roi);
     
                        
                        // Find foci population inside Astro region
                        Objects3DIntPopulation fociInAstroRegionPop = tools.findFociInAstro(astroRegion, fociPop);
                        System.out.println("Foci found in astrocyte region = "+fociInAstroRegionPop.getNbObjects());
                        // Find min foci volume
                        double minFociVol = tools.findMinFociVolume(fociPop);
                        System.out.println("Min foci volume (µm3) =" + minFociVol);
                                                
                        // Find foci population inside soma astrocyte
                        Objects3DIntPopulation fociInSomaAstroPop = tools.findFociInAstro(astroSoma, fociPop);
                        System.out.println("Foci found in soma astrocyte = "+fociInSomaAstroPop.getNbObjects());
                        
                        // Find foci population inside vessel
                        Objects3DIntPopulation fociInVesselPop = tools.findFociInAstro(vesselObj, fociPop);
                        System.out.println("Foci found in vessel = "+fociInVesselPop.getNbObjects());
                
                        // Write parameters
                        double astroRegionVol = tools.measureObjVolume(astroRegion);
                        int fociRegion = (int)Math.round(tools.measurePopVolume(fociInAstroRegionPop) / minFociVol);
                        double astroSomaVol = tools.measureObjVolume(astroSoma);
                        int fociSoma = (int)Math.round(tools.measurePopVolume(fociInSomaAstroPop) / minFociVol);
                        double vesselVol = tools.measureObjVolume(vesselObj);
                        int fociVessel = (int)Math.round(tools.measurePopVolume(fociInVesselPop) / minFociVol);
                        
                        outPutResults.write(rootName+"\t"+roi.getName()+"\t"+astroRegionVol+"\t"+fociRegion+"\t"+astroSomaVol+"\t"+fociSoma
                            +"\t"+vesselVol+"\t"+fociVessel+"\n");
                        outPutResults.flush();
                        
                        // save images objets
                        tools.drawObjects(astroSoma, vesselObj, fociPop, imgVessels, outDirResults, rootName, roi);
                        tools.flush_close(imgNuc);
                        tools.flush_close(imgAstro);
                        tools.flush_close(imgDots);
                        tools.flush_close(imgVessels);
                    }
                    else 
                       System.out.println("No astrocyte cell found !");
                }
            }
            outPutResults.close();
        }
        catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(RNA_In_Astro_Vessels.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done ...");
    }
}
