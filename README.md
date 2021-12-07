# StrikeSlipDef
This programme calculates the displacement due to a propagating dislocation in an elastic half space. StrikeSlipDef then creates cells from this finely discritizes fault trace assigning a width to the slipping surface that is user defined as is the depth. The code then calculates the surface displacment at given location using Okada(1992) subroutine for each cell. Note: the contributions from each event are not summed together, this is done later in the post-processing step.  StrikeSlipDef reads in an input file and produces outputs 3 files which are detailed below. StrikeSlipDef can be found and compiled in the src folder using the makefile.

StrikeSlipDef makes the following assumptions:
- Deformation due to dislocation is purely elastic
- Fault dip is vertical and follows surface trace
- Slip of 1m assumed for both shear and tensile dislocation
- Site of interest is near the longitude 27 Degrees (UTM zone 35)

Conventions:
    x : horizontal displacement (positive in east direction)
    y : horizontal displacement (positive in north direction)
    z : vertical displacement  (positive in up direction)


## Input files

Input file for StrikeSlipDef entitled 'input_file.txt' in this the input parameters are:

depth :             Depth to top of slipping cells
                    Unit : meter

rake :              Orientation of shear dislocation
                    Standard convention applied (i.e.right lateral slip = 180d).
                    Unit= Degrees

slipping_width :    width of slipping zone
                    Unit: meters

fault_file :        name of file containing discritized fault trace
                    Case provided is Marmara fault based on Şengör et al. (2014) discretization

tekr_lon, tekr_lat: location where surface displacement is calculated.
                    Unit: Degrees

## Output
Three files are produced by StrikeSlipDef, all filenames contain the depth of the source
in place of the asterisk  :

obs_tekr_\*.dat     : contains the surface displacement at position  pzn_lon,pzn_lat due to 1m shear dislocation at each source.
                      Format is as follows :
                      line 1 : number of sources
                      line 2 : number of observation points
                      line 3 : spacing between observation points
                      line 4 : lon , lat
                      line 5 : displacement in x direction (values along the line represent the displacement due to different sources )
                      line 6 : displacement in y direction (values along the line represent the displacement due to different sources )
                      line 7 : displacement in z direction (values along the line represent the displacement due to different sources )

obs_tekr_ten_\*.dat  :same format as obs_tekr_\*.dat but displacement is due to tensile dislocations

fault_latlon.txt:   file containing geographical coordinates of the corners for each source cell.
                    Each line is a source, with the first position being the top left (west) and proceeding clockwise, the final value is the strike of the cell:
                    
                    long1, lat1, depth1, long2, lat2, depth2, long3, lat3, depth3, long4, lat4, depth4, strike


## Postprocessing:

Matlab script plot_results.m reads output from StrikeSlipDef and produces figures


## References

Şengör, A. M. C. et al. The geometry of the North Anatolian transform fault
in the Sea of Marmara and its temporal evolution: implications for the development
of intracontinental transform faults. Canadian Journal of Earth Sciences 51, 222-242 (2014)

Okada, Y., 1992, Internal deformation due to shear and tensile faults in a half-space,
Bull. Seism. Soc. Am., 82, 1018-1040.

Georef source code:  https://github.com/andherit/georef, doi : 10.5281/zenodo.840875
