// Minimum Noise Fraction (MNF) Transformation

/*This script is maintained by the UHPSI Lab (http://highplainsstewardship.com/) at Yale University's School of Forestry and Environmental 
Studies. The script performs the Minimum Noise Fraction (MNF) Transformation. The input is a single image, and the output is an image in 
which the data is transformed and the new bands are reordered by image quality (decreasing).*/


/*-------------------------------------------------------------------*/
// MNF Process Overview
/*-------------------------------------------------------------------*/

/* Introduction:
The Minimum Noise Fraction Transformation (MNF) is used to reduce noise in an image. This transformation must first estimate the 
underlying noise extant in a dataset, which is done by at-sensor measurements or by post-collection approximation methods. This 
script uses a post-collection approximation method referred to as a "shift difference calculcation", which was chosen over other 
options given the ability to perform it automatically on any image and without requiring extra input from the user. The algorithm 
will search for a homogenous area within the scene and will compute the differences in reflectance values between neighboring pixels in 
this area; the basis for this calculation is as follows: the most homogenous area in the scene should exhibit variation between 
neighboring pixels due (mostly) to noise, thus the arithmetic differences between pixel reflectance values can be used as a basis for 
approximating the minimum value of noise in the image. (If a less homogenous area is used, the danger exists of over-correcting for noise 
by performing the same shift difference calculation in an area that exhibits greater spectral variation due to changes in the landscape 
rather than noise.)

After noise estimation via the shift difference calculations, the MNF transformation undergoes two PCA-type rotations to complete the 
minimization. These operations are described line-by-line in the script.

The final product of the script is an image with transformed bands ordered by image quality (decreasing). (The final image contains the same 
number of bands that the input image contained.) 

For background, consider exploring the following paper on the maximum noise fraction transformation (a precursor to the current 
formulation of the MNF transformation):
http://ieeexplore.ieee.org/xpl/articleDetails.jsp?reload=true&tp=&arnumber=3001&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D3001

*/


/*-------------------------------------------------------------------*/
// MNF Script Outline
/*-------------------------------------------------------------------*/
/*
1. User Defined Options.................................................(lines 56-146)
	1.1. Choose Input Data................................................(lines 63-78)
    1.2. Spectral Subset Option.........................................(lines 79-104)
	1.3. Band Selection for deriving Shift Difference area................(lines 105-123)
2. Formatting Reserves (i.e. palettes)..................................(lines 124-135)
3. Image Preparation....................................................(lines 147-209)
4. Homogenous Region Calculation and Selection..........................(lines 210-320)
5. Perform the Shift Difference Calculation on the Homogenous Region....(lines 321-352)
6. Perform the required mathematics to finalize the MNF Transformation..(lines 353-439)
7. Display the first three MNF bands and chart the eigenvalues..........(lines 440-461)
*/




/*-------------------------------------------------------------------*/
//Section 1: User Defined Options
/*-------------------------------------------------------------------*/
// Choose location to display in playground
centerMap(-106.54251, 44.56938, 12);


/*--------------------------------------------------------*/
// Section 1.1 Choose Input Data

// ****Choose data from the GEE Data Catalog or load your own data from MapsEngine****
// N.B. Input your dates such that your image of interest is the first within the image collection.
var collection = ee.ImageCollection('LANDSAT/LC8_L1T_8DAY_TOA').filterDate('2014-07-15', '2014-07-31');
addToMap(collection, {}, "Original Imagery");

// ****Input a polygon feature defining your study area****
var studyarea = ee.FeatureCollection('ft:1GxT2Q22KIGDjLFY1J11IX8lTAQPIoNyvBsG2fX0');
// addToMap(studyarea, {opacity:0}, 'Study Area');

// ****Input the known resolution of you dataset (in meters)***
var resolution = ee.Number(30);


/*--------------------------------------------------------*/
// Section 1.2 Spectral Subset Option
/*
On which band range(s) would you like to perform the transformation?
Define them, pairwise (inclusive), in the array below.

Ideally, all bands from an image would be used. However, if you're performing the analysis on imagery with many bands, 
(e.g. hyperspectral data) consider subsetting the data to lessen the computational intensity of the algorithm.

E.g.    [[1,10]] would select bands 1-10;
		
        [[2,5],
        [8,12],
        [16,20]] would select bands 2-5,8-12, and 16-20;
        
        [[1,1],
        [2,2],
        [3,3]] would select bands 1, 2, and 3;

*/

var brarray = [
	[1, 7]
];




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
	[2, 2],
	[3, 3],
	[4, 4]
];




/*-------------------------------------------------------------------*/
// Section 2. Formatting Reserves
/*-------------------------------------------------------------------*/
// Preferred color schemes can be designed here for later use

var paletteRMS = ['3366FF', '000033', 'ffffff', '3366CC', '000000'];

var palette_blue2purp = ['ff00ff', '006eff', '00ffd0', '459619'];




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

// Select the first image from the collection
var chosenimage = ee.Image(collection.first());

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

// Use the halper function from above to form a concatenated list of image bands.
var catbands = bandrangecat(brarray);

// Format the concatenated list of image bands as an array.
var catbandsarray = ee.Array(catbands).subtract(1).toList();

// Use the halper function from above to form a concatenated list of noise bands.
var catnoisebands = bandrangecat(noisebrarray);

// Format the concatenated list of noise bands as an array.
var catnoisebandsarray = ee.Array(catnoisebands).subtract(1).toList();


/*--------------------------------------------------------*/
// Format input data


// Spectrally subset input data
var originalImage = chosenimage.select(catbandsarray).clip(studyarea.geometry());

// Get band names of bands in use
var bands = originalImage.bandNames();

// Compute the total number of image bands and store this number
var numberofbands = bands.length();

// Subset the image for the noise shift difference calculations
var SubsetNoiseImage = chosenimage.select(catnoisebandsarray).clip(studyarea.geometry());

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
// addToMap(stdDev, {min:0, max:0.1}, 'Std Dev');
// addToMap(stdDev.select([1]), {min:0, max:1000,palette:['ff00ff','006eff', '00ffd0', '459619']}, 'Std Dev');

// Compute the quadratic mean (root mean square / RMS) of the neighborhood standard deviation through all bands for each pixel
// and then sum these values for each pixel
var RMS = stdDev.multiply(stdDev).reduce(ee.Reducer.sum()).divide(SubsetNoiseImage.bandNames().length()).sqrt();
// addToMap(RMS, {min:0, max:1}, 'RMS');

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
var desiredpercentage = 0.05;
var percentofrange = dictrange.multiply(desiredpercentage);
var bottompercent = dictrmsmin.add(percentofrange);
var threshold = RMS.select(['sum']).lt(bottompercent);

// Display pixels that have variance below the defined threshold
var lowRMS = threshold.mask(threshold);
// addToMap(lowRMS, {min:0, max:1, palette:palette_blue2purp}, 'Pixels with low RMS', false);


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
// addToMap(connected,{}, 'Connected Low RMS Islands');

// Determine the minimum number of pixels in the areas of interest
var pixelmin = ee.Algorithms.If(numberofbands.lte(50), ee.Number(50), numberofbands);

// Compute the number of connected pixels in each island
var count = connected.connectedPixelCount(ee.Array(pixelmin).multiply(4), true);
// addToMap(count, null, 'Connected Count');

// Reproject the layer and disregard all areas below the minimum pixel threshold
var precount = count.reproject(
	chosenimage.projection().crs(),
	null,
	resolution).gt(ee.Number(pixelmin));
// addToMap(precount, null, 'Connected Count Reprojected');

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
// addToMap(islandsprep, {color:'e51bc7'}, 'Islands', false);

var islandsWithArea = islandsprep.map(function(f) {
	return f.set('area', f.geometry().area(5))
});
// addToMap(withArea);

var islands = islandsWithArea.filterMetadata('area', 'greater_than', resolution.multiply(resolution).multiply(numberofbands));
// addToMap(islands, {color: '900000'}, 'Big Islands');


/*--------------------------------------------------------*/
// Sum each region's variance then sort the rßegions

//Find the total variance within each box and choose the box with the least variance
var buffvariance = RMS.reduceRegions(islands, ee.Reducer.sum(), resolution);
// addToMap(buffvariance,{color:'00ff00'},'Buffered Sums of RMS St.d.');

// Select the buffered area with the lowest variance to get shift difference area
var buffvarImage = buffvariance.limit(1, 'sum', true);
// addToMap(buffvarImage, {color:'ff0000'}, 'Final Buffered Area');




/*-------------------------------------------------------------------*/
// Section 5: Perform the Shift Difference Calculation
/*-------------------------------------------------------------------*/

// Clip the image to the chosen shift difference area
var kernelarea = originalImage.clip(buffvarImage.geometry());
// addToMap(kernelarea,{},'Noise Kernel Area',false);


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
// addToMap(diff, {min:-100, max:100}, 'diff');




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
// addToMap(arrayimage);

// Find the mean value of the bands in the whole image
var meanimage = originalImage.reduceRegion('mean', studyarea.geometry(), resolution, null, null, false, 800000000);

// Make an array from the image’s band means for each band
var meanarray = ee.Array(meanimage.values(bands));

// Mean correct the image
var meancenteredimage = arrayimage.subtract(meanarray);
// addToMap(meancenteredimage, {}, "Mean Centered Image");

// Multiply the mean centered image by the noise eigenvectors then scale the data by the noise standard deviation values
var nwarrayimage = meancenteredimage.arrayRepeat(1, 1).arrayTranspose()
	.matrixMultiply(eigenvectors.transpose())
	.matrixMultiply(matrixr);
// addToMap(nwarrayimage, {min:0, max:10}, "Noise Whitened Array Image");

// Covariate the noise-whitened dataset
var nwimage = nwarrayimage.arrayFlatten([
	['Noise-Whitened'], bands
]);
// addToMap(nwimage);
var nwcovardict = nwarrayimage.arrayProject([1]).reduceRegion(ee.Reducer.covariance(), studyarea.geometry(), resolution, null, null, false, 800000000);

// Eigendecompose the covariance matrix
var nwcovariancematrix = ee.Array(nwcovardict.get('array'));
var nweigendecomplist = nwcovariancematrix.eigen();

// Retrieve the eigenvalues and eigenvectors for each MNF transformed band
var nweigenvectors = nweigendecomplist.slice(1, 1);
var nweigenvalues = nweigendecomplist.slice(1, 0, 1);

// Finalize the MNF Transformation by multiplying the second eigenvector matrix by the noise-whitened data
var mnfdata = nwarrayimage.matrixMultiply(nweigenvectors.transpose());
// addToMap(mnfdata);

// Use a map function to retrieve the band names for the image
var bl = ee.List.sequence(1, numberofbands, 1);
var fbands = bl.map(function(n) {
	return ee.String('Band ').cat(ee.Number(n).int());
});

// Flatten the array image back into a normal image and add results to map
var mnfimage = mnfdata.arrayFlatten([
	['MNF Transformed'], fbands
], " ");
addToMap(mnfimage, {
	min: 0,
	max: 10
}, "MNF Transformed Data");

// **End of the MNF Transformation**

var mnfcovardict = mnfdata.arrayProject([1]).reduceRegion(ee.Reducer.covariance(), null, resolution, null, null, false, 3000000000);
var mnfcovariance = ee.Array(mnfcovardict.get('array'));
// print(mnfcovariance);
var mnfeigendecomp = mnfcovariance.eigen();
var mnfeigenvalues = mnfeigendecomp.slice(1, 0, 1);




/*-------------------------------------------------------------------*/
// Section 7: Display MNF bands and chart the eigenvalues
/*-------------------------------------------------------------------*/

var eigenValueArray = ee.Array(mnfeigenvalues).repeat(0, 1);

var charty = Chart.array.values(eigenValueArray, 0, bl).setSeriesNames(["Eigen Values"]);
charty = charty.setOptions({
	title: 'Eigenvalues For MNF Bands',
	hAxis: {
		title: 'MNF Bands'
	},
	vAxis: {
		title: 'Eigenvalue'
	},
	lineWidth: 1,
	pointSize: 4,
	series: {
		0: {
			color: 'darkgreen'
		}
	}
});
print(charty);
