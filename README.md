# Quasi-Optics MatLab Library
 MatLab library used for computing and simulating quasi-optics and standard EM systems and problems. The library is developed as part of the assignments' solutions for the EE4725 Quasi-Optical Systems course in TU Delft.

## Functions and Objects
 This section describes the implemented functions and objects in the library.

### Coordinate Systems Transformations
 In EM problems, the orientation of the E and H-fields depend on the observation point. Furthermore, it is often the case that the E and H-fields must be calculated at many observation points. However, MatLab's existing functions (`cart2sphvec` and `sph2cartvec`) compute for a single observation point at a time. Hence, in this library new functions (`cart2sph_vector` and `sph2cart_vector`) are implemented, in order to compute for many observation points. Finally, functions for converting coordinates between systems are implemented with the goal of improving workspace and code readability.

 #### Cartesian-To-Spherical Vector Conversion (`cart2sph_vector`)
 This function converts vector in cartesian coordinate system to a vector in spherical coordinate system.  
 The function takes M-by-N-by-3 matrix for the cartesian vector, and M-by-N matrix for Theta and Phi. The cartesian matrix's third dimension represents the cartesian unit vectors as x, y, and z respectively. The M-by-N matricies of each cartesian unit vector holds its magnitude for the corresponding vector location at the elevation and azimuth. The function returns M-by-N-by-3 matrix for the spherical vector with third dimensions representing the radial, elevation angle, and azimuth angle unit vectors respectively.

 #### Spherical-To-Cartesian Vector Conversion (`sph2cart_vector`)
 This function converts vector in spherical coordinate system to a vector in cartesian coordinate system.  
 The function takes M-by-N-by-3 matrix for the spherical vector, and M-by-N matrix for Theta and Phi. The spherical matrix's third dimension represents the spherical unit vectors as radial, elevation angle, and azimuth angle vectors respectively. The M-by-N matricies of each spherical unit vector holds its magnitude for the corresponding vector location at the elevation and azimuth. The function returns M-by-N-by-3 matrix for the cartesian vector with third dimensions representing the x, y, and z unit vectors respectively.

 #### Cartesian-To-Spherical Coordinate Conversion (`cart2sph_cord`)
 This function is implemented to simplify the workspace and improve the computation/simulation code's readability.

 Wrap function of cart2sph, which returns the spherical coordinates in a 3D matrix.  
 The function takes the coordinate points separately or in M-by-N-by-3 matrix. For a M-by-N-by-3 matrix input, the third dimension represents  x, y, and z coordinates respectively. For other input formats, the  coordinates are either specified or in the order x, y, and z. The  function returns the spherical coordinates in a M-by-N-by-3 matrix with third dimension representing radial distance, elevation, and azimuth angle respectively.

 #### Spherical-To-Cartesian Coordinate Conversion (`sph2cart_cord`)
 This function is implemented to simplify the workspace and improve the computation/simulation code's readability.

 Wrap function of sph2cart, which returns the spherical coordinates in a 3D matrix.  
 This function takes the coordinate points separately or in M-by-N-by-3 matrix. For a M-by-N-by-3 matrix input, the third dimension reperesents the radial distance, elevation angle, and azimuth angle respectively.  For other input formats, the coordinates are either specified or in the order radial distance, elevation angle, and azimuth angle. The function returns the cartesian coordinates in a M-by-N-by-3 matrix with the third dimension representing x, y, and z respectively.

#### Miscellaneous
 This section describes other implemented miscellaneous functions in the library.

 #### Meshgrid Grids in One Matrix (`meshgrid_comb`)
 This function is implemented to simplify the workspace and improve the computation/simulation code's readability. 

 Wrap function of meshgrid, which returns the grid matricies together in a 3D or 4D matrix.  
 The function takes either 2 or 3 vectors. For 2 vectors, it returns 3D matrix with the third dimension holding the vector grids in the order of parsing them to the function. For 3 vectors, it returns 4D matrix with the forth dimension holding the vector grids in the order of parsing them to the function.

 #### UV Representation (`uv_repr`)
 This function calculates UV coordinate representation.  
 This function takes the elevation and azimuth angle and calculates the UV representation of the provided spherical coordinate points. For a 3D matrix input, the third dimension represents the elevation and azimuth angles respectively. For other input formats, the coordinates are either specified or in the order elevation and azimuth angle. The function returns the UV representation in a M-by-N-by-2 matrix with third dimension representing U and V respectively.

 #### Magnitude Normalization (`norm_magnitude`)
 This function normalizes the magnitude to the maximum magnitude.  
 The function takes the magnitude (scalar, vector, or matrix) and returns the normalized to the maximum magnitude. The output can be in dB or in linear scale. If scale is not specified (by 'dB' or 'linear'), the default is 'dB' scale.
