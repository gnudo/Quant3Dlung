# Quantitative lung data analysis tutorial

## 1. Introduction

The following package contains the quantitative morphological and topological analysis tools for high-resolution 3D lung image data as introduced in Lovric _et al._ [PLoS ONE **12** (9), e0183979 (2017), DOI: [10/cddm](http://doi.org/cddm)]. The present version is still a little "patchy": it allows to reproduce all the results from the mentioned article, but also requires a number of manual steps to be performed for a proper setup. However, improvements are on the way.

This tutorial is divided into three parts, whereas the 1st part (*Segmentation*) produces data that can either be used by the **2nd** (*Thickness map analysis*) or the **3rd** part (*Curvature*).


## 2. Segmentation

Let's assume the volume data is saved within single *tif*-files in a directory called `<data-folder-X>`. To segment the volume data according to the method described in [the paper](http://doi.org/cddm), follow the steps below:
1.   Go to `segmentation` subfolder
2.   Go to `Matlab` subfolder.
3.   Compile the Matlab c-function for your platform by running in Matlab: `mex morpho_operation.c`
   4.   Run `segmenter.m` script in Matlab with: `segmenter('<data-folder-X>/','minpeakheight','minpeakprominence','maxwidth','false')`
   5.   Alternatively, run from the terminal with: `matlab -nodisplay -r "segmenter <data-folder-X>/ 75 10 15 false;quit";`
6.   Adjust the variables *minpeakheight*, *minpeakprominence* and *maxwidth* according to your data and re-run the script. The data will be placed in differently named directories. The variable *minpeakprominence* is usually the most sensitive one.
7.   The last argument (*true* or *false*) defines whether the segmenter script saves all intermediate steps or not. The usage is optional - we standardly use *false*.

After that, the data can be found in `<data-folder-X>_mphXX_mppYY_mwZZ` subfolder. Finally, it is necessary to run a 3D connected component analysis to remove "small pixel noise"-parts. To do this:

1.   Run `connected_component.m` script in Matlab with: `connected_component.m('<data-folder-X>_mphXX_mppYY_mwZZ/')`
2.   Alternatively, run from the terminal with: `matlab -nodisplay -r "connected_component <data-folder-X>_mphXX_mppYY_mwZZ/"`
3.   Finally, there will be a new subfolder found in `<data-folder-X>_mphXX_mppYY_mwZZ` called *labeled*, which in the following we call `<seg-data-folder-X>`.



## 3. Thickness map analysis

We assume an already segmented (binarized) dataset that is saved as a tif-stack consisting of single tif-files in a folder called `<seg-data-folder-X>`. There are a few pre-conditions to be met before running/automatizing the following on a random SLS-machine such as *GWS-3* (not all steps are necessary on a normal PC):

*   Download and use your own *Fiji*-installation. The best is to store it in your *e-account* and **either** reset the `fiji` symbolic link in the system path **or** replace the *fiji*-command in `thickcalculator.py` with `/path-to-your-Fiji-installation/ImageJ-linux64`
*   Make sure the single-instance listener is switched on (in *Fiji*, goto `Edit`->>`Options`->>`Misc`)
*   Use the *epd_free* python version: `module add xbl/epd_free/7.3-2-2013.06`
*   Make sure that the Python module *psutils* is installed. If you don't have *root* access, install it locally with `pip install --user psutils`


After all these steps have been pre-checked, you should be save to run (also multiple times from several *bash*-environments):

     python thickcalculator.py <seg-data-folder-X>

or we do the following steps manually:
1.   Open Fiji
2.   Load images stack from `<seg-data-folder-X>`: `File`->>`Import`->>`Image sequence`
3.   Calculate *Local Thickness*: `Analyze`->>`Local Thickness`->>`Local Thickness (Complete process)`
4.   At the end save the file in parent directory as a 32bit RAW file with the name `<seg-data-folder-X>__LocThk.raw`

### 3.1. Histogram calculation (binning) -- not used

One way is then to calculate *binned* histograms from the raw 3D data for which the following script `thickmap/Matlab/local_thickness_csv_driver.m` can be used by calling:

    matlab -nodisplay -r "local_thickness_csv_driver <seg-data-folder-X>__LocThk.raw X Y Z;quit";

or by calling from within matlab:

    local_thickness_csv_driver('<seg-data-folder-X>__LocThk.raw','X','Y','Z')

where *X*,*Y* and *Z* are the coordinates of the 3D thickness map data.


### 3.2. Density calculation

The disadvantage of histograms is that the binning is rather arbitrary, the results can be noisy and it can be harder to compare different histograms with each other quantitatively. For this reason we use the method of *kernel density* calculation, applied directly on the `RAW`-thicknessmap data. For this, there are two scripts available for *Matlab* and *R*, whereas *R* offers by far better calculation capabilities/packages and faster calculations.

#### 3.2.1. Matlab -- not used

The *Matlab* script `kerneldensity.m` requires a `RAW`-thicknessmap dataset and the volume sizes (*x*,*y*,*z*) to be given. It then runs the Matlab function `ksdensity` with different *bandwidth* settings. The scripts is only a small feasibility demonstration and not used in the final images.


#### 3.2.2. R

The *R*-script `thick_density.R` follows the approach from the original [paper](http://doi.org/cddm). The approach is to investigate the influence of different conditions on the same dataset and then calculate the confidence levels for each dataset separately so that we can see whether the quantitative results are influenced by *finite size*-effects of *segmentation* parameters. Let's assume we have three lung pressures to investigate (10, 20, 30 cmH2O) and thus three original datasets. Now we run the segmentation algorithm with different parameters and apply the Thickness map analysis on each segmented volume. Thus we end up having multiple volumes per pressure. To do the density calculation in *R* do the following:
1.   Create a helper file `data_info_<X>.txt` (see example) and put all absolute paths to the different datasets, obtained from different conditions
2.   In `thick_density.R`, adapt the path to `data_info.txt`
3.   Adapt the volume size *data_xyz*
4.   Adapt the pixel size
5.   Adapt the density bandwith (standardly: 2 or 4)
6.   Adapt the plots
7.   Run from within *R*

Another script (`thick_density_standalone.R`) has been created in case *R* runs out of memory while doing the density calculation. In that case, the density calculation can be run on another machine with more memory producing a temporary binary file, which in return is used by `thick_density.R` just for plotting.

## 4. Curvature analysis


Again, we start by assuming an already segmented (binarized) dataset that is saved as a tif-stack consisting of single tif-files in a folder called `<seg-data-folder-X>`. As before, there are a few pre-conditions to be met before running/automatizing:


*	 Download and compile MEPP with the given patch files (see Section 4.2.1).
*   Download and use the latest version of *VTK*. The best is download the source and compile it from scratch will all python-bindings.
*   Make sure to export the respective python paths, by adapting the correct folder names (i.e. of the VTK installation) and running: `source curvature/vtkpython-gws-sl6.sh`

### 4.1. Air-to-tissue surface extraction

The air-to-tissue surface extraction is basically conducted in two steps by using *Fiji* (from within Python) in combination with *VTK*. Run the following steps:
1.  In the binarized (segmented) `<seg-data-folder-X>`-data, make sure that *lung tissue* is **white** and *air* is **black** in terms of gray values.
2.  Run from a terminal: `python curvature/Python/createRaw.py <seg-data-folder-X>/ [<destination-folder>/]`
3.  Likewise, run: `python curvature/Python/convertRaw2Ply.py <seg-data-folder-NAME>.raw`

While **Step 2** will create a *RAW*-file `<filename>_X_Y_Z.raw` in either the parent directory of the given `<destination-folder>`, **Step 3** will do the actual surface-extraction and produce a *PLY*-file (in the so-called *Stanford Triangle Format*).


### 4.2. Curvature calculation

The actual curvature calculation is conducted in the [MEPP](https://liris.cnrs.fr/mepp/) software. Since MEPP needs to be modified and *re-compiled* in order to give correct results, first all necessary compilation steps for Scientific Linux 6 are explained. Those are valid for other systems as well such as Ubuntu/Fedora/etc. In a next step the actual curvature calculation is explained.


#### 4.2.1. MEPP compilation steps

**(1)** First install all necessary libraries:

    yum install boost-devel.x86_64 xerces-c-devel.x86_64
    yum install libavutil50.x86_64 gmp-devel.x86_64 mpfr.x86_64 mpfr-devel.x86_64
    yum install libavcodec52.x86_64 libavformat52.x86_64 libavutil49.x86_64 libavutil50.x86_64 libpostproc51.x86_64 libswscale0.x86_64

**(2)** Download [CGAL](http://www.cgal.org/) and build with cmake:

    cmake . -DCMAKE_INSTALL_PREFIX=/<Download-location>/CGAL -DWITH_examples=true -DWITH_demos=true

**(3)** Download libqglviewer and (e.g. libQGLViewer-2.4.0.tar.gz) and compile/build (with `qmake` and `make`), and use QT4 with `/usr/bin/qmake-qt4 PREFIX=/<QGL-Destination>/QGLViewer`

**(4)** Install *cmake-gui* through YUM.

**(5)** Download and configure/patch *MEPP*:

    git clone https://github.com/MEPP-team/MEPP.git
    git-apply correct_curvature.patch  ## correct curvature signs
    git-apply mepp-deactivate-viewer   ## deactivate viewer to speed-up
    git-apply mepp-sl_6.patch          ## export calculated curvature values

**(6)** Compile *MEPP*:

    cd MEPP
    mkdir build; cd build/



#### 4.2.2. Curvature calculation

1.   Open *PLY*-file in MEPP (and watch increase in occupied RAM-memory)
2.   Laplacian smoothing: Goto `Tools/Processing`->`Various Processing`->`Laplacian smoothing`
3.   Set *Deformation factor: 0.5* and *Iteration number: 10*. Additionally enable *Preserve boundaries*.
4.   Run Curvature analysis by: Goto `Analysis/Filtering`->`Curvature`->`Calculate (Normal Cycle Algorithm`
5.   Set *Method: Geodesic* and *Curvature radius: 3.5*. A bigger radius may be chosen when aiming at higher precision for very low curvature values, but a lower *Curvature radius* only makes sense if the mesh precisely matches the volume.
6.   Save the file as *OBJ* for latter processing.

The results are then saved in `<Destination>` as set in ... or in the patch-file, respectively as *CSV* file with the following structure

    X, Y, Z, kappa_1, kappa_2

where *(X,Y,Z)* are the position vectors of each vertex and *kappa_1* and *kappa_2* the minimal and maximal curvature respectively.



### 4.3. Distribution analyses

Since curvatures are calculated per vertex and each vertex is surrounded by a particular surface area patch, the areas of the surface patches need to be taken into account as weighting factors when calculating the 2D distribution of surface curvatures. First, however, the surface patch areas have to be extracted and added to the *CSV*-file. Here, a two-step approach is done:

1.   Create a new *CSV*-file by adding the respective surface area patches: `python curvature/Python/mapCellArea2VertexArea.py <FILE.OBJ> <ORIG-CURV.csv>`
2.   As a result, two new files will be created:
    * `<ORIG-CURV_NEW.csv>` (incorporating the surface patch area per vertex in addition to the two curvature values, as stated above)
    * `<ORIG-CURV_CURV.ply>` (a new *PLY*-file with the curvature values for 3D rendering and visualization).

The two different distribution analyses are then further processed in *R*.



#### 4.3.1. Interface shape distribution (ISD)

ISDs represent 2D probability density plots of the minimumum and maximum curvatures of a surface mesh. They are calculated with the script `ISD_final.R` where following steps/adjustements need to be undertaken:

1.   Adjust settings (range, nbin)
2.   Save as PDF (with Cairo) + set size: 8.97 x 6.76
3.   (Note that if R crashes it just needs to be restarted -- origin of bug unknown)
4.   Run in terminal `convert -density 600 -quality 85 -depth 8 <ISD_plot>.pdf <ISD_plot>.png`
5.   Save from within Inkscape or similar as *EPS*

#### 4.3.2. Gauss and Mean curvatures

Simply run `curvature/GaussMeanCurvature.R` and adjust it like above.

## 5. Validation

The validation includes the creation of synthetic data, i.e. spheres of different size and distribution that are embedded into a pre-defined volume. Subsequently, the thickness map and the curvature analysis can be run on this volumetric dataset and the results compared with the ground truth. In the created datasets, balls have gray value 255 and the bulk material has gray value 0.

### 5.1. Data creation

Run the script `createvolume.m` from the `validation` subfolder and adjust the following:

   * `resolution`: how big should the data cube be
   * `R_balls`: (array of) different ball radii in pixels
   * `p_balls`: (array of) percentual amount of ball counts
   
### 5.2. Thickness map analysis

Run the thickness map analysis as described in Sec. 2. Afterwards run the script `validation_thickness.R` to plot the results. Make sure to check the following points, before executing the script:

   * The thickness map result should be located in the folder and named: `fullValidation_LocThk.raw`.
   * Input the correct data size values in *data_xyz*.
   * Change the variable *px_size* accordingly.

By running the script, the

   1. calculated thickness map values will be loaded and the PDE will be calculated, and
   2. ground truth data will be loaded from `quant_data.txt`.

Finally , both results are plotted for comparison.

### 5.3. Curvature analysis

For verifying the curvature results, the best way is to check whether the *measured* curvatures actually correctly represent the sphere radii as well as the percentual distribution (related to the surfaces of the different spheres radii) within the volume.

For doing so, we first conduct the curvature analysis as described in Sec. 3 (for which we use a curvature radius of 7.5 in the normal cycle algorithm). We then run the script `validation_curvature.R` which:

   * extracts the different sphere radii from the obtained curvature values and calculates the PDE
   * loads the ground truth data and convolves it with a Gaussian in order to mimic uncertainties in the curvature calculation which originates from imperfect meshes (obtained from pixel data).
   

