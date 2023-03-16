# RNA_In_Astro_Vessels

* **Developed for:** Marc
* **Team:** Cohen-Salmon
* **Date:** Mars 2023
* **Software:** Fiji


### Images description

3D images taken with a x63 objective, then deconvoluted

4 channels:
  1. *405:* Nuclei
  2. *488:* Vessels
  3. *561:* astrocytes
  4. *642:* foci
  
A *.roi* or *.zip* file containing ROI(s) must be provided with each image.

### Plugin description

In each ROI:
* Detect nuclei with Cellpose
* Detect astrocyte soma with Cellpose
* Find astrocyte soma max colocalization with nucleus
* Detect vessels with a remove outliers to remove microgila cells processes, Median filtering + thresholding + close filtering
* Find foci with stardist RNA-foci-1.2 model
* Detect foci in outlined roi, in vessel and in soma astrocyte
* Give objects volume and foci number as fovi vol / min foci vol
* 

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Cellpose** conda environment + *cyto* model
* **Stardist** 

### Version history

Version 1 released on March 16, 2023.
