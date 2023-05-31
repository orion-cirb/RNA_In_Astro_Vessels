# RNA_In_Astro_Vessels

* **Developed for:** Marc
* **Team:** Cohen-Salmon
* **Date:** March 2023
* **Software:** Fiji


### Images description

3D images taken with a x63 objective, then deconvoluted.

4 channels:
  1. *405:* Nuclei
  2. *488:* Vessels
  3. *561:* Astrocytes
  4. *642:* Foci
  
A *.roi* or *.zip* file containing ROI(s) must be provided with each image.

### Plugin description

In each ROI:
* Detect nuclei with Cellpose
* Detect astrocytes somas with Cellpose
* Only keep astrocytes somas colocalizing with a nucleus
* Detect vessels with Remove outliers (to remove microglia cells processes) + Median filtering + thresholding + closing
* Find foci with Stardist
* Count foci in vessels, astrocytes somas and in the rest of the ROI
* Give vessels ans astrocytes somas volume + foci number in each compartment as foci total volume/minimum foci volume

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Cellpose** conda environment + *cyto* model
* **Stardist** conda environment + *RNA-foci-1.2* model

### Version history

Version 1 released on March 16, 2023.
