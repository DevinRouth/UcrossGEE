// Mixture Tuned Match Filtering (MTMF) Classification (MTMF_V.1.0)

/*This script is maintained by the UHPSI Lab (http://highplainsstewardship.com/) at Yale University's School of Forestry and Environmental 
Studies. It performs the Mixture Tuned Matched Filtering (MTMF) classification.*/


/*-------------------------------------------------------------------*/
// MTMF Process Overview
/*-------------------------------------------------------------------*/

/* Introduction:
///////pasted from MF explanation on line 587 of MTMF script:/*###### 
The Mixture Tuned Matched Filtering process has two main inputs: an image and an endmember (either field/lab collected or a pure pixel from 
the image). The algorithm returns two values for each pixel of the input image: a Matched Filtering score and a Mixture Tuning score (otherwise 
known as an infeasibility score).  

The Matched Filtering score represents how closely a pixel matches the endmember on a 0 to 1 range where 0 is least like the 
endmember and 1 indicates a strong match. This score also represents relative subpixel abundance of the endmember in each pixel, with the decimal
value between 0 and 1 translating to a subpixel abundance percentage, e.g. an MF score of 0.23 would equate to a predicted subpixel endmember 
abundance of 23%. Note: though it is mathematically possible to compute MF scores that are <0 or >1, this occurrence is rare and usually 
indicates low spectral variance in the data. (Values <0 or >1 do not allow for the normal subpixel abundance interpretation.)

The Infeasibility Score, derived from the mixture tuning process, ranges from 0 to an indefinite maximum value. Mathematically, infeasibility scores
are the geometric distance 
The likelihood of a false positive increases as the infeasibility 
increases.  Therefore, pixels with a high MF score and a low infeasibility score are those most likely to contain the 
endmember.
//////
*/


/*-------------------------------------------------------------------*/
// MTMF Algorithm Table of Contents
/*-------------------------------------------------------------------*/

/*
1. User Defined Options...................................................(lines 56-146)
	1.1. Define the Study Area............................................(lines 63-78)
    1.2. Define the Endmember.............................................(lines 82-107)
	1.3. Choose Input Imagery.............................................(lines 108-126)
	1.4. Format Input Imagery.............................................(lines 108-126)
2. Formatting Reserves (i.e. palettes)....................................(lines 127-138)
3. Image Preparation......................................................(lines 150-209)
4. Homogenous Region Calculation and Selection............................(lines 210-320)
5. Perform the Shift Difference Calculation on the Homogenous Region......(lines 321-352)
6. Finalize the MNF Transformation........................................(lines 353-439)
7. Display the first three MNF bands and chart the eigenvalues............(lines 440-461)
8. Derive the MF Scores...................................................(lines 
9. Derive Infeasibility Scores............................................(lines
Addendum 1: Optional Context..............................................(lines
Addendum 2: Optional Charts...............................................(lines
*/




/*-------------------------------------------------------------------*/
//Section 1: User Defined Options
/*-------------------------------------------------------------------*/
// Choose location to display in playground
Map.setCenter(-106.54251, 44.56938, 12);
Map.setStyle('satellite');


/*--------------------------------------------------------*/
// Section 1.3 Choose Input Data

// ****Choose data from the GEE Data Catalog or load your own data from MapsEngine****

// If using a single image...
//		- Insert ID for a single image of your study area
var singleImage = ee.Image('LC8_L1T/LC80440342013154LGN00');

var Ref_2014153 = ee.Image('GME/images/08039105425737821391-14148998190422555071');
var prepmultiply = Ref_2014153.multiply(10000);
var mainImage153 = prepmultiply.mask(prepmultiply).divide(10000);

// Define the image on which you would like to perform the analysis
var chosenImageprep = mainImage153;

// ****Input the known resolution of you dataset (in meters)***
var resolution = ee.Number(30);



/*--------------------------------------------------------*/
// Section 1.4 Format Input Data

// ****Select bands for a spectral subset
//    Define them, pairwise (inclusive), in the array below.
var brarray = [
	[15, 20]
];
// N.B Ideally, all bands from an image would be used. However, if you're performing the analysis on imagery with many bands, 
// (e.g. hyperspectral data) consider subsetting the data to lessen the computational intensity of the algorithm. Furthermore,
// if particular bands of a dataset are known to be "bad", you can choose not to select them (and thus remove them from the algorithm).

/*
E.g.    [[1,10]] would select bands 1-10;
		
        [[2,5],
        [8,12],
        [16,20]] would select bands 2-5,8-12, and 16-20;
        
        [[1,1],
        [2,2],
        [3,3]] would select bands 1, 2, and 3;
*/


/*--------------------------------------------------------*/
// 1.3. Band Selection for deriving the Shift Difference area
/*
Below is a separate band subset selection array used to subset bands for use in the shift difference area derivation.
Ideally, the shift difference would be performed using all of the selected bands from above, but given the computational intensity
of the calculation within GEE the script times-out if it is run on more than a few bands. The current suggestion is to select 
the true color bands. (Most analysts performing an MNF transformation would select a shift difference area by visually inspecting 
the image with these bands; as such, the true color bands for LandSat 8 appear below as the default values.)

If you would like to make changes, define them in a pairwise fashion (inclusive) in the array below.
*/

var noisebrarray = [
	[15, 15],
	[17, 17],
	[20, 20]
];




/*--------------------------------------------------------*/
// Section 1.1 Define the Study Area

// ****If you would like to import your own polygon from a fusion table, define it here...
var importedStudyArea = ee.FeatureCollection('ft:1GxT2Q22KIGDjLFY1J11IX8lTAQPIoNyvBsG2fX0');

// ****If you would you like to draw your own polygon, enter the coordinates of the vertices below...
var drawnStudyArea = ee.FeatureCollection(ee.Feature(ee.Geometry.Polygon(
	[
		[-89.4795, 29.4654],
		[-89.4164, 29.4552],
		[-89.4823, 29.19],
		[-89.5489, 29.2056]
	])));


// ****Choose study area variable
var studyarea = importedStudyArea; // Choose either importedStudyArea or drawnStudyArea (as above)
// Map.addLayer(studyarea, {opacity:0}, 'Study Area', false);


/*--------------------------------------------------------*/
// Section 1.2 Define Endmember

// Non-Image Based Endmembers
// We provide options for field/lab-collected endmembers in addition to image-based endmembers. Though we do not yet 
// include a pixel purity index to select the purest pixels for image-based endmember selection, we hope to include 
// this in the future.


// ****If using a an image-based endmember (drawn or from a training region)...input 1 below
// ****If using a custom endmember....................................input 2 below

var endmemberchoice = 1;


// ****If you chose 1, proceed to line 101. If you chose 2, proceed to line 115.

// If you would like to import your own polygon from a fusion table, define it here...
var importedtrainingregion = ee.FeatureCollection('ft:1zsmCtvfmAag_0A27MwXkEA6e3474yUm3uEp_oUIp');

// If you would you like to draw your own polygon to generate an endmember, enter the coordinates of the vertices below...
var drawnTrainingRegion = ee.FeatureCollection(ee.Feature(ee.Geometry.Polygon(
	[
		[-89.45849, 29.36348],
		[-89.4554, 29.36116],
		[-89.45832, 29.35885],
		[-89.46158, 29.36027]
	])));

// ****Choose training region variable as indicated above
var trainingregion = importedtrainingregion; // Choose either importedtrainingregion or drawnTrainingRegion (as above)


// ****If you chose 2...

// Input the reflectance values of the endmember in the list below. Be sure to match the number and order 
// of the input data bands from the imagery or the algorithm will not operate correctly.
//-89.478707, 29.328061
var customEndmember = ee.Array([4155, 4272, 4099, 4364, 4427, 4403 // Endmember reflectance values at each band from first image
	// Endmember reflectance values at each band from next image
	// Etc.
]); // Etc.
// Note: the number of rows in this list should correspond to the number of images being used; the number of values in each 
// line should correspond to the number of bands selected for the analysis.




// ~~~

///////////////////
// What's your infeasibility score threshold?
/* This number may vary greatly between different data and classifications.  The default below corresponds with 
the default image and endmember inputs.  Typically, the threshold can range from close to zero to ≈20.  The UHPSI 
lab is in the midst of developing an automatic threshold selection algorithm, which would operate regardless of 
the data inputs and without any addition user knowledge, to help achieve the best results in identifying presence 
and abundance of the endmember.*/
var userinfthreshold = 0.05;





/*-------------------------------------------------------------------*/
// Section 2. Formatting Reserves
/*-------------------------------------------------------------------*/
// Preferred color schemes can be designed here for later use

var paletteRMS = ['3366FF', '000033', 'ffffff', '3366CC', '000000'];
var palette_blue2purp = ['ff00ff', '006eff', '00ffd0', '459619'];
var MFpalette = ['8c2a04', '800000', 'ff0000', 'ffA500', 'ffff00', '808000', '008000', '00ff00'];




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
///////////////////////NO USER INPUT AFTER THIS POINT///////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////




/*-------------------------------------------------------------------*/
// Section 3: Image Preparation
/*-------------------------------------------------------------------*/

// This function performs the band range selection using the arrays defined in the section above.
var bandrangecat = function(inputarray) {

	var openlist = [];

	var numofranges = inputarray.length - 1;

	for (var i = 0; i <= numofranges; i++) {

		var start = inputarray[i][0];
		var end = inputarray[i][1];

		for (var j = start; j <= end; j++) {
			openlist.push(j);
		}
	}
	return openlist
};

// Use the helper function from above to form a concatenated list of image bands.
var catbands = bandrangecat(brarray);

// Format the concatenated list of image bands as an array.
var catbandsarray = ee.Array(catbands).subtract(1).toList();

// Use the helper function from above to form a concatenated list of noise bands.
var catnoisebands = bandrangecat(noisebrarray);

// Format the concatenated list of noise bands as an array.
var catnoisebandsarray = ee.Array(catnoisebands).subtract(1).toList();


/*--------------------------------------------------------*/
// Format input data


// Spectrally subset input data
var originalImage = chosenImageprep.select(catbandsarray).clip(studyarea.geometry());
Map.addLayer(originalImage, null, 'Original Image Clipped');

// Get band names of bands in use
var bands = originalImage.bandNames();

// Compute the total number of image bands and store this number
var numberofbands = bands.length();

// Subset the image for the noise shift difference calculations
var SubsetNoiseImage = chosenImageprep.select(catnoisebandsarray).clip(studyarea.geometry());

// Retrieve a list of noise band names
var noisebands = SubsetNoiseImage.bandNames();

// Compute the total number of noise bands and store this number
var numberofnoisebands = noisebands.length();



/*-------------------------------------------------------------------*/
// Section 4: Homogenous Region Calculation and Selection
/*-------------------------------------------------------------------*/


// Format a kernel to use when finding the homogenous area
var square_kernel = ee.Kernel.square(resolution.multiply(3), "meters");

// Find standard deviation for each pixel and its determined neighborhood
var stdDev = SubsetNoiseImage.reduceNeighborhood(ee.Reducer.stdDev(), square_kernel);
// Map.addLayer(stdDev, {min:0, max:0.1}, 'Std Dev');
// Map.addLayer(stdDev.select([1]), {min:0, max:1000,palette:['ff00ff','006eff', '00ffd0', '459619']}, 'Std Dev');

// Compute the quadratic mean (root mean square / RMS) of the neighborhood standard deviation through all bands for each pixel
// and then sum these values for each pixel
var RMS = stdDev.multiply(stdDev).reduce(ee.Reducer.sum()).divide(SubsetNoiseImage.bandNames().length()).sqrt();

// Find and store the minimum and maximum variance values (and their range) within the area of interest
var RMSDict = RMS
	.reduceRegion({
		reducer: ee.Reducer.minMax(),
		geometry: studyarea.geometry(),
		scale: resolution
	});
var dictrmsmax = ee.Number(RMSDict.get('sum_max'));
var dictrmsmin = ee.Number(RMSDict.get('sum_min'));
var dictrange = dictrmsmax.subtract(dictrmsmin);


/*--------------------------------------------------------*/
// Find the area with the lowest variance

// Display the quadratic mean layer using the computed minimum and maximum
var RMS_vis = RMS.visualize({
	min: dictrmsmin,
	max: dictrmsmax,
	palette: palette_blue2purp
});
// addToMap(RMS_vis,{},'RMS Full Range', false);

// Define the threshold of how much variance is acceptable for consideration for a shift difference area
var desiredpercentage = 0.04;
var percentofrange = dictrange.multiply(desiredpercentage);
var bottompercent = dictrmsmin.add(percentofrange);
var threshold = RMS.select(['sum']).lt(bottompercent);

// Display pixels that have variance below the defined threshold
var lowRMS = threshold.mask(threshold);
// Map.addLayer(lowRMS, {min:0, max:1, palette:palette_blue2purp}, 'Pixels with low RMS', false);


/*--------------------------------------------------------*/
// Make each area an island

// Define the kernel used in the connectedComponents() call
var kernel_fixed = ee.Kernel.fixed(3, 3, [
	[1, 1, 1],
	[1, 0, 1],
	[1, 1, 1]
]);

// Connect all pixels that meet the threshold defined above
var connected = lowRMS.connectedComponents(kernel_fixed, 256);
// Map.addLayer(connected,{}, 'Connected Low RMS Islands');

// Determine the minimum number of pixels in the areas of interest
var pixelmin = ee.Algorithms.If(numberofbands.lte(50), ee.Number(50), numberofbands);

// Compute the number of connected pixels in each island
var count = connected.connectedPixelCount(ee.Array(pixelmin).multiply(4), true);
// Map.addLayer(count, null, 'Connected Count');

// Reproject the layer and disregard all areas below the minimum pixel threshold
var precount = count.reproject(
	originalImage.projection().crs(),
	null,
	resolution).gt(ee.Number(pixelmin));
// Map.addLayer(precount, null, 'Connected Count Reprojected');

// Run a mask to leave only the islands of interest then restore the bands
var countFV = precount.mask(precount);
var RTV = countFV.addBands(countFV);

// Turn each island into a vector feature and filter out islands that do not
// meet the minimum area requirement
var islandsprep = RTV.reduceToVectors({
	reducer: ee.Reducer.first(),
	geometry: studyarea.geometry(),
	geometryType: 'polygon',
	scale: resolution,
	bestEffort: true
});
// Map.addLayer(islandsprep, {color:'e51bc7'}, 'Islands', false);

var islandsWithArea = islandsprep.map(function(f) {
	return f.set('area', f.geometry().area(5))
});
// Map.addLayer(withArea);

var islands = islandsWithArea.filterMetadata('area', 'greater_than', resolution.multiply(resolution).multiply(numberofbands));
// Map.addLayer(islands, {color: '900000'}, 'Big Islands');


/*--------------------------------------------------------*/
// Sum each region's variance then sort the regions

//Find the total variance within each box and choose the box with the least variance
var buffvariance = RMS.reduceRegions(islands, ee.Reducer.sum(), resolution);
// Map.addLayer(buffvariance,{color:'00ff00'},'Buffered Sums of RMS St.d.');

// Select the buffered area with the lowest variance to get shift difference area
var buffvarImage = buffvariance.limit(1, 'sum', true);
// Map.addLayer(buffvarImage, {color:'ff0000'}, 'Final Buffered Area');




/*-------------------------------------------------------------------*/
// Section 5: Perform the Shift Difference Calculation
/*-------------------------------------------------------------------*/

// Clip the image to the chosen shift difference area
var kernelarea = originalImage.clip(buffvarImage.geometry());
// Map.addLayer(kernelarea,{},'Noise Kernel Area',false);


// Define kernels that link a pixel to its neighbors immediately above and to the left of it
var kernel_left = ee.Kernel.fixed(3, 3, [
	[0, 0, 0],
	[1, 0, 0],
	[0, 0, 0]
]);

var kernel_up = ee.Kernel.fixed(3, 3, [
	[0, 1, 0],
	[0, 0, 0],
	[0, 0, 0]
]);

// Create a layer stack of niehgboring pixel values in order to perform math between pixel values
var kernelimage_left = kernelarea.neighborhoodToBands(kernel_left);
var kernelimage_up = kernelarea.neighborhoodToBands(kernel_up);

var diff_left = kernelimage_left.subtract(kernelarea);
var diff_up = kernelimage_up.subtract(kernelarea);

// Find average difference between pixels for the whole shift difference area
var diff = diff_left.add(diff_up).divide(2).clip(kernelarea.geometry());
// Map.addLayer(diff, {min:-100, max:100}, 'diff');



/*-------------------------------------------------------------------*/
// Section 6: Finalize the MNF Transformation
/*-------------------------------------------------------------------*/

// Find the covariance matrix of the finalized shift difference area
var covardict = diff.toArray().reduceRegion(ee.Reducer.covariance(), null, resolution, null, null, false, 800000000);

// Convert the covariance matrix into an array
var noisecovariancematrix = ee.Array(covardict.get('array'));
// print("Noise Coavariance Matrix", noisecovariancematrix);

// Decompose the matrix into eigenvalues and eigenvectors
var eigendecomplist = noisecovariancematrix.eigen();


/*----------------------------------------*/
// MNF Matrix

// Use the results of the decomposition to formulate the required matrices for the subsequent mathematics
var eigenvalues = eigendecomplist.slice(1, 0, 1);
var eigenvectors = eigendecomplist.slice(1, 1);

var matrixr = eigenvalues.sqrt().pow(-1).matrixToDiag();
var matrixcmnf = eigenvalues.pow(-1).matrixToDiag();


/*--------------------------------------------------------*/
// Noise-whiten the dataset

// Convert the image to an array
var arrayimage = originalImage.toArray();
// Map.addLayer(arrayimage);

// Find the mean value of the bands in the whole image
var meanimage = originalImage.reduceRegion('mean', studyarea.geometry(), resolution, null, null, false, 800000000);

// Make an array from the image’s band means for each band
var meanarray = ee.Array(meanimage.values(bands));

// Mean correct the image
var meancenteredimage = arrayimage.subtract(meanarray);
// Map.addLayer(meancenteredimage, {}, "Mean Centered Image");


// Multiply the mean centered image by the noise eigenvectors then scale the data by the noise standard deviation values
var nwarrayimage = meancenteredimage.arrayRepeat(1, 1).arrayTranspose()
	.matrixMultiply(eigenvectors.transpose()) //took out transpose()
	.matrixMultiply(matrixr);
// Map.addLayer(nwarrayimage, {min:0, max:10}, "Noise Whitened Array Image");

var nwcovardict = nwarrayimage.arrayProject([1]).reduceRegion(ee.Reducer.covariance(), studyarea.geometry(), resolution, null, null, false, 800000000);

// Eigendecompose the covariance matrix
var nwcovariancematrix = ee.Array(nwcovardict.get('array'));
var nweigendecomplist = nwcovariancematrix.eigen();

// Retrieve the eigenvalues and eigenvectors for each MNF transformed band
var nweigenvectors = nweigendecomplist.slice(1, 1);
var nweigenvalues = nweigendecomplist.slice(1, 0, 1);

// Finalize the MNF Transformation by multiplying the second eigenvector matrix by the noise-whitened data
var mnfdata = nwarrayimage.matrixMultiply(nweigenvectors.transpose()); //took out .transpose()
// Map.addLayer(mnfdata);

// Use a map function to retrieve the band names for the image
var bl = ee.List.sequence(1, numberofbands, 1);
var fbands = bl.map(function(n) {
	return ee.String('Band ').cat(ee.Number(n).int());
});

// Flatten the array image back into a normal image and add results to map
var mnfimage = mnfdata.arrayFlatten([
	['MNF Transformed'], fbands
], " ");
// Map.addLayer(mnfimage, {min:0, max:10}, "MNF Transformed Data");

// **End of the MNF Transformation**

var mnfcovardict = mnfdata.arrayProject([1]).reduceRegion(ee.Reducer.covariance(), null, resolution, null, null, false, 3000000000);
var mnfcovariance = ee.Array(mnfcovardict.get('array'));
// print(mnfcovariance);
var mnfeigendecomp = mnfcovariance.eigen();
var mnfeigenvalues = mnfeigendecomp.slice(1, 0, 1);



// /*-------------------------------------------------------------------*/
// // Section 7: Display MNF bands and chart the eigenvalues
// /*-------------------------------------------------------------------*/

// var eigenValueArray = ee.Array(mnfeigenvalues).repeat(0,1);

// var charty = Chart.array.values(eigenValueArray, 0, bl).setSeriesNames(["Eigen Values"]);
// charty = charty.setOptions({
//   title: 'Eigenvalues For MNF Bands',
//   hAxis: {
//     title: 'MNF Bands'
//   },
//   vAxis: {
//     title: 'Eigenvalue'
//   },
//   lineWidth: 1,
//   pointSize: 4,
//   series: {
//     0: {color: 'darkgreen'}
//   }
// });
// print(charty);




/*-------------------------------------------------------------------*/
// Section VI: Transform the Endmember into MNF Space
/*-------------------------------------------------------------------*/

// Prep then project the endmember vector into MNF space
/*Note: the endmember array and the mean array must be adjusted so as to have 
  the appropriate dimensions in order for the linear algebra to operate correctly.*/



/*--------------------------------------------------------*/
// Subset and clip the endmember data

// Spectral subset of each endmember image in the input collection
var emImage = chosenImageprep.select(catbandsarray).clip(trainingregion.geometry());
// Map.addLayer(emImage, null, 'Bands Selected for Endmember Image');

// Find the mean reflectance value of each band from the training region
var imageEnd = emImage.reduceRegion(ee.Reducer.mean(), trainingregion.geometry(), resolution);
// print(imageEnd);


/*-------------------------------------------------------------------*/
// Section 4: Prepare endmember for analysis
/*-------------------------------------------------------------------*/

// Call on user’s choice of whether to use field or images-based endmember.
var finalprep = ee.Algorithms.If(endmemberchoice == 1, imageEnd.toArray(), customEndmember);
// print('final prep', finalprep);

// Make the unidimensional endmember array a two dimensional array.
var endmember = ee.Array(finalprep).repeat(1, 1);
// print('Final Endmember', endmember);


//Mean center the endmember
var meancenteredendmember = endmember.subtract(meanarray.repeat(1, 1));

//Noise whiten the endmember
var nwendmember = meancenteredendmember
	.transpose()
	.matrixMultiply(eigenvectors.transpose())
	.matrixMultiply(matrixr);


// print('meancenteredendmember',meancenteredendmember.transpose());
// print('eigenvectors',eigenvectors);  //this is used for both in eric one
// print('eigen transposed', eigenvectors.transpose());
// print('meancenteredimage',meancenteredimage.arrayRepeat(1,1).arrayTranspose());

//Finalize the MNF transformation by multiplying the endmember array with the noise-whitened image eigenvectors
var targetspectra = nwendmember.matrixMultiply(nweigenvectors.transpose()).transpose();
// print("MNF Transformed Endmember (Target Spectra)", targetspectra);



/*-------------------------------------------------------------------*/
// Section VII: Derive MF Score
/*-------------------------------------------------------------------*/

var mnfbands = mnfdata.arrayFlatten([
	['MNF Bands'], bands
]).bandNames();

//Find the covariance matrix of the MNF transformed data
var mnfcovardict = mnfdata.toArray().arrayProject([1])
	.reduceRegion(ee.Reducer.covariance(),
		null, null, null, null, false, 800000000);
//print('MNF Covariance', mnfcovardict);

//Decompose the matrix into eigenvalues and eigenvectors 
var mnfcovariancematrix = ee.Array(mnfcovardict.get('array'));
var mnfeigendecomplist = mnfcovariancematrix.eigen();
var mnfeigenvalues = mnfeigendecomplist.slice(1, 0, 1);
// print('MNF Eigen Values', mnfeigenvalues);

//Derive the diagonoal matrix of MNF eigenvalues for the MF calculation
var matrixcmnf = mnfeigenvalues.pow(-1).matrixToDiag();
// print('MNF Matrix', matrixcmnf);


/*----  Calculate Pixel MF Scores -------*/

// Perform the final calculations to compute MF scores at each pixel
var mfscores = mnfdata.matrixMultiply(matrixcmnf
	.matrixMultiply(targetspectra)
	.divide(targetspectra.transpose()
		.matrixMultiply(matrixcmnf)
		.matrixMultiply(targetspectra)
		.get([0, 0])));

Map.addLayer(mfscores.arrayProject([0]).arrayFlatten([
	["MF Score"]
]), {
	min: -1,
	max: 1,
	palette: MFpalette
}, "MF Scores");




/*--------------------------------------------------------*/
// Section VIII: Derive Infeasibility Score
/*--------------------------------------------------------*/

// Perform the calculations to compute Infeasibility scores at each pixel
var inffinal = mnfdata.subtract(mfscores.matrixMultiply(targetspectra.transpose()))
	.matrixMultiply(mnfdata.subtract(mfscores.matrixMultiply(targetspectra.transpose()))
		.arrayTranspose())
	.sqrt()
	.arrayProject([0])
	.arrayFlatten([
		["Infeasibility"]
	])
	.clip(studyarea.geometry())
	.divide(mfscores.matrixMultiply(mnfeigenvalues.sqrt().subtract(
				ee.Array.identity(numberofbands)
				.matrixDiagonal())
			.transpose())
		.multiply(-1)
		.add(mnfeigenvalues.sqrt().transpose()).pow(2).matrixMultiply(mfscores.matrixMultiply(mnfeigenvalues.sqrt().subtract(
					ee.Array.identity(numberofbands)
					.matrixDiagonal())
				.transpose())
			.multiply(-1)
			.add(mnfeigenvalues.sqrt().transpose()).pow(2).arrayTranspose())
		.sqrt()
		.arrayProject([0])
		.arrayFlatten([
			[" "]
		])
		.clip(studyarea.geometry()));

Map.addLayer(inffinal, {
	min: 0,
	max: 0.2,
	palette: ['b6f430', '0c1744']
}, "Infeasibility Scores", false);

// Threshold out high infeasibility scores, which indicate a potential false positive
var infThreshold = inffinal.gt(userinfthreshold);
// Map.addLayer(infThreshold, null, 'Binary Infeasibility Score');

// Create a mask to cover all of the false positives
var lowInf = infThreshold.mask(infThreshold);
Map.addLayer(lowInf, {
	min: 0,
	max: 1,
	palette: ['b6f430', '491010']
}, 'Pixels with high Infeasiblity');


// ~~~









/*--------------------------------------------------------*/
// Addendum 1: Optional Context
/*--------------------------------------------------------*/

// Insert and display option contextual polygons, such as an outline of the study area, training regions used for
// endmember selection, or ground truth plots.

// Paint an outline of the study area
var paintImage = ee.Image(0).mask(0);
Map.addLayer(paintImage.paint(studyarea, '3300ff', 2), null, 'Outline of Study Area');

// Optional display of the training regions used for image based endmember
Map.addLayer(trainingregion, {
	opacity: 0
}, 'Endmember Training Region(s)', false);




// // /*--------------------------------------------------------*/
// // // Addendum 2: Optional Charts
// // /*--------------------------------------------------------*/

// // /*This section takes sample points from different landcover types around the study site and plots their spectral profiles onto 
// // a chart alongside the profiles of the image- or custom-endmember options.  The charts only show bands used in the 
// // classification after the image has been spectrally subset.  In the future, these spectral profiles could be averages of 
// // several points or the average reflectance values of polygonal training regions either drawn in real time or uploaded from 
// // field-collected shapefiles.*/ 

// // // Create random points across study area (default)...
// // var randomPoints = ee.FeatureCollection.randomPoints(studyarea, 3);
// // //Map.addLayer(randomPoints, null, 'Random Chart Points');

// // // Manually select points to chart....
// // //		N.B. To manually select points from certain locations in the study site, clikc the inspector, click on a point on the map,
// // //		and copy and paste the coordinates from the console into the three lines below.
// // var landcoverOne = ee.Feature(ee.Geometry.Point(-106.5184, 44.6161),{'label': 'One'});
// // var landcoverTwo = ee.Feature(ee.Geometry.Point(-106.5083, 44.5747),{'label': 'Two'});
// // var landcoverThree = ee.Feature(ee.Geometry.Point(-106.5060, 44.5777),{'label': 'Three'});

// // // Cast the sample points into a Feature Collection.
// // var manualPoints = ee.FeatureCollection([landcoverOne, landcoverTwo, landcoverThree]);
// // //Map.addLayer(manualPoints, null, 'Maual Chart Points');

// // //****Choose chart points variable
// // var chartPoints = randomPoints; // Choose either randomPoints or manualPoints

// // var collection = ee.ImageCollection.fromImages(originalImage);
// // // Make a list of reflectance values of the points
// // var info = collection.getRegion(chartPoints, resolution);
// // // print(info);

// // // Format the data intolists in order to chart values pulled from the image with values from
// // // the endmember feature collection.
// // var bandList = ee.List(info.get(0)).slice(4);
// // // print('Band Names', bandList);

// // var wavelengthList=  ee.List.sequence(1, bandList.length());
// // var wavelengthArray = ee.Array(wavelengthList).repeat(1,1);
// // // print('Wavelength List', wavelengthArray);

// // var oneList = ee.List(info.get(5)).slice(4);
// // var oneArray = ee.Array(oneList).repeat(1,1);
// // // print('Landcover one', oneArray);

// // var twoList = ee.List(info.get(1)).slice(4);
// // var twoArray = ee.Array(twoList).repeat(1,1);
// // // print('Landcover two', twoArray);

// // var threeList = ee.List(info.get(3)).slice(4);
// // var threeArray = ee.Array(threeList).repeat(1,1);
// // // print('Landcover three', threeArray);

// // // Concatenate the lists of each cover type and cast them into an array.
// // var arraysconcat = ee.Array.cat([oneArray,twoArray,threeArray,endmember],1);
// // // print("Concatenated Arrays", arraysconcat);

// // // Chart the values.
// // var arrayChart = Chart.array.values(
// //     arraysconcat,0, wavelengthArray).setSeriesNames(["landcoverOne","landcoverTwo","landcoverThree","Endmember"]);
// // arrayChart = arrayChart.setOptions({
// //   title: 'Spectral Profiles at three points in the study area',
// //   hAxis: {
// //     title: 'Stacked Image Bands'
// //   },
// //   vAxis: {
// //     title: 'Reflectance Value'
// //   },
// //   lineWidth: 2,
// //   pointSize: 1,
// //   series: {
// //     0: {color: 'darkgreen'},
// //     1: {color: 'lightgreen'},
// //     2: {color: 'lightblue'},
// //     3: {color: 'red'},
// //   }
// // });
// // print(arrayChart);
// // print('If you would like to manually select the');
// // print('location of the above points, see line 448.');
// // // ~~~
