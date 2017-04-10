// Spectral Angle Mapper (SAM) Classification (SAM_V.1.0)

/*This script is maintained by the UHPSI Lab (http://highplainsstewardship.com/) at Yale University's School of Forestry and Environmental 
Studies. The script performs the Minimum Noise Fraction (MNF) Transformation. The input is a single image, and the output is an image in 
which the data is transformed and the new bands are reordered by image quality (decreasing).*/


////////////////////// <--- Denotes User Input Options


/*-------------------------------------------------------------------*/
//SAM Process Outlined
/*-------------------------------------------------------------------*/
/*
1. Input Image 
    1.1. Should be in reflectance values (not radiance)
2. Select Bands to use
3. Select Endmember
4. MNF Transformation (optional)
    4.1. Image
    4.2. Endmember
5. Subset of MNF transformed bands (not yet included)
6. Spectral Angle Calculation
7. View Results and Select Threshold
    7.1. Threshold in radians (default results are in degrees)
*/




/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
//                        Classification
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

centerMap(-106.54251, 44.56938, 12);



/*-------------------------------------------------------------------*/
//Section I: User Defined Options
/*-------------------------------------------------------------------*/

/*###### 
This section allows the user to select a set of bands for subsetting both the image and the endmember.  It can be likened to 
the “Bad Bands List” in ENVI, where the user provides the  “good bands”, and the program will automatically use those spectral 
bands in any further processing.  The default bands we have omitted match the “no data” bands identified by USGS/NASA as well 
as bands known to cause issues due to atmospheric absorption of moisture. For the images provided here, we have identified an 
additional set of bands that contain striping, cross track, and Smile effect issues, but we have chosen not to subset those out 
currently.  Because of this, both the Smile effect and striping can be seeing in the final results.
######*/

//////////////////////
/* Image Spectral Subset
    On which band ranges would you like to perform the analysis?
    Define them, pairwise (inclusive), in the array below.*/

var brarray = [
	[11, 35],
	[82, 93],
	[95, 97],
	[102, 114],
	[135, 144],
	[147, 163],
	[189, 189],
	[202, 213]
];





/*###### 
There is a separate band selection section used to subset bands for use in the shift difference area derivation. Ideally, the 
shift difference would be performed using all of the selected bands for the image, but given the nature of that calculation 
within GEE, the script freezes if it is run on more than a few bands.  The current  fix for this issue is to select the true 
color bands—most analysts performing an MNF transformation would select a shift difference area by looking at the image with 
these bands.  As such, we replicate an accepted practice.
######*/

////////////////////////
/* Band Selection for Shift Difference Calculation
    Though this operation can be performed on any number of bands, we recommend using the true color bands, which appear 
    below as the default values. If you would like to make changes, define them in a pairwise fashion (inclusive) in the 
    array below. While the script is still in beta, do not select more than three bands.*/
var noisebrarray = [
	[11, 11],
	[21, 21],
	[29, 29]
];

// ~~~ 




/*###### 
Here we provide the option to choose the data on which to perform the operation.  While this analysis is most suited to 
hyperspectral data, such as that from the Hyperion sensor on board the EO-1 satellite, there are currently only three 
Hyperion images available in GEE, two of which we have uploaded through MapsEngine. However, SAM can be applied to any 
imagery in the GEE catalogue.  In the following example, we mosaicked the two images to cover our study site, which makes 
this portion of the script a little messy..
######*/

// Temporary images from which to select:
var Ref_2014153 = ee.Image('GME/images/08039105425737821391-14148998190422555071');
var Ref_2014161 = ee.Image('GME/images/04725300055333973066-04594044960761723980');
var GEESample = ee.Image('SAMPLES/HYPERION_SAMPLE/EO1H0210402010195110KF');

/*Mask non image pixels by increasing their values (which are in reflectance) via multiplication, masking the imagery, then 
dividing the images by the original multiplicative value. If masked without the multiplication first, a partial mask is 
created, and the image appears partially transparent. */
var prepmultiply = Ref_2014153.multiply(10000);
var mainImage153 = prepmultiply.mask(prepmultiply).divide(10000);
//addToMap(mainImage153.select(['b23']), {min:0, max:1, /*palette:['ff0000','0000ff']*/}, 'East 153');

var prepmultiply161 = Ref_2014161.multiply(10000);
var mainImage161 = prepmultiply161.mask(prepmultiply161).divide(10000);
//addToMap(mainImage161.select(['b23']), {min:0, max:1}, 'West 161');

//Create an image collection from the east and west images
var bothSides = ee.ImageCollection.fromImages([mainImage161, mainImage153]);
//addToMap(bothSides.select(['b23']));

//Mosaic East/West Image Collection
var mosaic = bothSides.mosaic();
addToMap(mosaic.select(['b29', 'b21', 'b11']), {
	min: 0,
	max: 0.6
}, 'Full Ranch Mosaic');

/////////////////////
/* Choose Input Image or Image Collection 
    Which image from above would you like to use?*/
var chosenimage = mosaic;

// ~~~




/*###### 
As noted above,the option is provided to insert a chosen endmember.  We provide options for both field/lab-collected 
endmembers and image-based endmembers.  Though we do not yet include a pixel purity index to select the purest pixels for
image-based endmember selection, we hope to include this soon.  Additionally, we'd like to modify the current selection to 
take the average of several points chosen (or what would be averaging the spectra of a training region.
######*/

/////////////////////
/* Endmember options
    Would you like to use an endmember from the spectral library or a image-based endmember?
    Input 1 for "library" or 2 for "pixel" below */
var endmemberchoice = 2;

// Spectral library options:
var UXfield_spurgeav_EO1 = ee.FeatureCollection('ft:1590F-lS4y4aWpvTkkae_830LuQxpP2ZzsOpmFoo7');
var UXlab_sagebouquet_EO1 = ee.FeatureCollection('ft:1_8_qrU_Rkfo-uM5E-vg9pKeF23dlapx_0OsU6122');
var UXfield_sagedryflwr_EO1 = ee.FeatureCollection('ft:1m_A-zKHjUqN-nRLOtF9zPtzqQIXrmp9feOd9pMEc');
var UXlab_testspurge_EO1 = ee.FeatureCollection('ft:1bIMEQ6ZRM1Ow0oqOEEtPnbYrZnmlRf61PpKZQbHW');
//(More options at https://ee-api.appspot.com/e463088961d68e08db29a3ae83b2cd72)

///////////////////
/* If you chose '1', copy and paste the variable name of your desired endmember from the above spectral library in the 
        variable “selectedlibraryendmember” line below.
   If you chose '2', input the coordinates for your pixel.
    ****NOTE! Leave both of the following variables in place (regardless of the fact that only one of them will be used).
*/
var selectedlibraryendmember = UXlab_testspurge_EO1;
var pointcoordinates = ee.Feature.Point(-106.4918, 44.5964);

// ~~~




/*###### 
The product of the SAM classification is a map where each pixel has a value for the angular distance in N-Dimensional space 
(where “N” equals the number of bands selected for the analysis) from the endmember.  Using training regions, a user should 
choose a threshold angle to determine to determine the highest angle value that a pixel can exhibit for it to be classified 
as the endmember type.  This threshold will impact the map’s accuracy assessment. 
######*/

///////////////////
// What's your threshold angle?
/* This is the maximum angle difference (in degrees) to allow for a pixel to classify it as the input endmember.  The 
default is 0.10 radians (or ~5 degrees), but the user can adjust this value based on knowledge of the endmember/data.*/
var userangthreshold = 6;

// ~~~




// Optional contextual layers to add to the map
var ranchoutline = ee.FeatureCollection('ft:1GxT2Q22KIGDjLFY1J11IX8lTAQPIoNyvBsG2fX0');
// Paint the ranch to an image so it can be transparent
var ranchimage = ee.Image(0).mask(0);

// Import field-collected data for reference
var july_waypnts = ee.FeatureCollection('ft:1etFXuQ73DuWMdqMImrqc1iTS1kOUgHy-sUAr8e_n');
var spurgetraining = ee.FeatureCollection('ft:1zsmCtvfmAag_0A27MwXkEA6e3474yUm3uEp_oUIp');

// ~~~




/*-------------------------------------------------------------------*/
// Section II: Helper functions
/*-------------------------------------------------------------------*/

/*###### 
When called, the helper functions loop through the bands of an input image in order to subset the images and endmembers 
(regardless of the band names).  One potential problem with these functions is that it is impossible to count bands in the 
proper order if the user’s input image has already been spectrally subset. A possible solution would be to use a Feature 
Collection filter and the ‘contains’ operator to filter for band names containing '12', '15', etc...
######*/

var bandrangecat = function(inputarray, inputimage) {

	var openlist = [];

	var numofranges = inputarray.length - 1;

	for (var i = 0; i <= numofranges; i++) {

		var start = inputarray[i][0] - 1;
		var end = inputarray[i][1];

		var slicedimage = inputimage.slice(start, end);

		openlist[i] = slicedimage;

	}
	var fullbandrange = ee.Image.cat(openlist);
	return fullbandrange;
};

var bandrangecat2 = function(inputarray) {

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


var catbands = bandrangecat2(brarray);
// print('Concatenated list of bands',catbands);

// ~~~




/*-------------------------------------------------------------------*/
// Section III: Image preparation
/*-------------------------------------------------------------------*/

/*###### 
Here the image needs to be subset to contain the chosen bands from the first section.
######*/

//Run the band helper function to subset the chosen image for the chosen bands.
var OriginalImage = bandrangecat(brarray, chosenimage);
// print(OriginalImage.getInfo());

// Get the band names of the subset image to perform a band count.
var bands = OriginalImage.bandNames();
//print('Chosen Image Bands', bands);

// Count the bands of the subset image.
var numberofbands = bands.length();

//Get the CRS and CRS transform of the chosen image for later use.
var maincrs = chosenimage.getInfo().bands[0].crs;
// print(maincrs);
var maincrs_transform = chosenimage.getInfo().bands[0].crs_transform;
// print(maincrs_transform);

// ~~~




/*-------------------------------------------------------------------*/
// Section IV: Add formatting reserves
/*-------------------------------------------------------------------*/

//Add formatting reserves.

//Make palettes to define subsequent visualization parameters.
var paletteRMS = ['3366FF', '000033', 'ffffff', '3366CC', '000000'];
var palette_blue2purp = ['ff00ff', '006eff', '00ffd0', '459619'];
var SAMpalette = ['00ff00', '008000', '808000', 'ffff00', 'ffA500', 'ff0000', '800000', '8c2a04'];

// ~~~




/*-------------------------------------------------------------------*/
// Section V: Prepare endmember for analysis
/*-------------------------------------------------------------------*/

/*###### 
The endmember needs to be subset to the same bands as the image in order for the analysis to work.  This section prepares 
both feature collection-based and image-based endmembers then determines which will be used based on the above user input.
######*/

// Prep then project the endmember vector into MNF space.
/*Note: the endmember array and the mean array must be adjusted so as to have 
  the appropriate dimensions in order for the linear algebra to operate correctly.*/

//Spectrally subset the library-based endmember using band ranges from the helper functions.
var libendmemberprep1 = selectedlibraryendmember.filter(ee.Filter.inList("Band", catbands)).sort("Band");
//Format reflectance values of the library-based endmember as an array.
var libendmemberprep2 = libendmemberprep1.aggregate_array('Pixel_Value');
var library = libendmemberprep2;

//Spectrally subset the pixel based endmember using bands from the helper function.
//Create an Image Collection to use in getRegion (which only takes collections).
var collection = ee.ImageCollection.fromImages([OriginalImage]);

//Produce a list of reflectance values in the image under the pixel-based endmember.
var pointendmemberprep1 = collection.getRegion(pointcoordinates, 30);
//print(pointendmemberprep1);

//Create a list of reflectance values from variable above.
var pointendmemberprep2 = ee.List(pointendmemberprep1.get(1)).slice(4);
// print('Pixel Endmember Reflectance Values', pointendmemberprep2);
var point = pointendmemberprep2;

//Call on user’s choice of whether to use field or images-based endmember.
var finalprep = ee.Algorithms.If(endmemberchoice == 1, library, point);

//Make the unidimensional endmember array a two dimensional array.
var endmember = ee.Array(finalprep).repeat(1, 1);
//print('Final Endmember', endmember);

// ~~~




/*-------------------------------------------------------------------*/
// Section VI: Spectral Angle Calculation
/*-------------------------------------------------------------------*/


/*###### 
SAM  project each pixel’s and each endmember’s values as a vector in N-dimensional space, where n equals the number of bands.
The angular distance between pixel and endmember vectors indicates how similar a pixel’s spectra is from the endmember 
spectra.

Because the angle is measured between the two vectors (the direction of  which represents the spectral feature or shape of 
the profile), the brightness of the inputs (reflected in the length of the vectors) will have relatively little impact on 
the match estimate.  Additionally, this means a lab-collected endmember can be used in the same manner as an image-collected 
or field-collected endmember.

The spectral angle can have values between 0 and 180 degrees (Pi radians).

Useful References:
      1.  Kruse et al., 1993  (Remote Sensing of the Environment)
          The Spectral Image Processing System (SIPS) - Interactive Visualization and Analysis of Imaging Spectrometer Data.
                                  
      2.  Chang, 2000 (IEEE)
          An Information-Theoretic Approach to Spectral Variability, Similarity, and Discrimination for Hyperspectral Image 
            Analysis. 
                        
      3.  Keshava, 2004 (IEEE Transactions on Geoscience and Remote Sensing)
          Distance Metrics and Band Selection in Hyperspectral Processing With Applications to Material Identification and 
            Spectral Libraries. 
######*/




//Format Image Variables.
//  t = image (test) values
//  r = endmember (reference) values

// Convert the standard imageto an array image.
var tValues = OriginalImage.toArray().arrayRepeat(1, 1);
//print('bigxtranspose',bigxtranspose.getInfo());

// Transpose the array image values.
var tValuesTranspose = tValues.arrayTranspose();
//print('Final Image Input Array', bigx.getInfo());



/*---- Format the Numerator and Denominator of the Vector Angle Formula --------*/

//Transpose the endmember array.
var rValues = endmember.transpose();
//print('Final Endmember Input Array', mu);

//Calculate the dot product of the array image and the endmember..
var tTimesrValuesmagnitude = tValuesTranspose.matrixMultiply(endmember);
//print('Multiplied Matrix', bigxtimesmutranspose);

//Calculate the magnitude of each vector in the array image..
var tValuesmagnitude = tValuesTranspose.matrixMultiply(tValues).sqrt();
//print(bigxmagnitude);

//Calculate the magnitude of the endmember vector..
var rValuesmagnitude = rValues.matrixMultiply(endmember).sqrt();
//print(mumagnitude);



//Set data for numerator.
var numerator = tTimesrValuesmagnitude;
//print('Numerator', numerator);

//Make denominator by multiplying the magnitudes of the array image vectors by the magnitude of the endmember. 
var denominator = tValuesmagnitude.matrixMultiply(rValuesmagnitude);
//print('Denominator', denominator);

//Divide the numerator by the denominator.
var fraction = numerator.divide(denominator);
//print('Fraction', fraction);

/*----- Take the arc-cosine of the fraction --------*/

//Take inverse cosine of the fraction.
var radangle = fraction.acos();
//print('Angle in radians', radangle);

// Convert radians to degrees.
var degangle = radangle.multiply(180).divide(Math.PI);
//addToMap(degangle);

/*------ Format images to map ----------*/

//Clip the final image to the ranchoutline.
var degreeangle = degangle.arrayGet([0, 0]).clip(ranchoutline);
//addToMap(degreeangle);

//Find the maximum and minimum angles in the image.
var reddict = degreeangle.reduceRegion({
	reducer: ee.Reducer.minMax()
		.combine(ee.Reducer.sampleStdDev(), '', true)
		.combine(ee.Reducer.mean(), '', true),
	geometry: ranchoutline.geometry(),
	scale: 30
});
//Turn the dictionary into an array.
var comarrdict = reddict.toArray();
//Get the outputs of the dictionary as individual values
var dictmean = ee.Number(comarrdict.get([1]));
var dictsd = ee.Number(comarrdict.get([3]));
var dictmax = ee.Number(comarrdict.get([0]));
var dictmin = ee.Number(comarrdict.get([2]));
var dictrange = dictmax.subtract(dictmin);

//Add the SAM Rule Image, where each pixel contains the value of its angular distance from the endmember.
var SAM_vis = degreeangle.visualize({
	min: dictmin,
	max: dictmax,
	palette: SAMpalette
});
addToMap(SAM_vis, null, 'Standard Values - Spurge');

//Add the classified image with the user defined threshold.
addToMap(degreeangle.gt(userangthreshold), {
	min: 0,
	max: 1,
	palette: SAMpalette
}, 'Spurge Present', false);


/* ---- Experiments in Scale ----*/
//This section uses the spread of the data to visualize the results.

//Add the SAM Rule Image, where each pixel contains the value of its angular distance from the endmember.
addToMap(degreeangle, {
	min: 0,
	max: 180,
	palette: SAMpalette
}, 'Full Possible Scale - Spurge', false);

//Display the Rule image, which is normalized to the standard deviation on a scale from -1 to 1.
var normdevImage = degreeangle.subtract(dictmean).divide(dictsd);
addToMap(normdevImage, {
	min: -1,
	max: 1,
	palette: SAMpalette
}, 'SAM Norm Dev', false);

//Display the Rule image, normalized to the range of the image’s angle values on a scale from 0 to 1.
var norm01Image = degreeangle.subtract(dictmin).divide(dictrange);
addToMap(norm01Image, {
	min: 0,
	max: 1,
	palette: SAMpalette
}, 'SAM Norm Range', false);


// ~~~

alert('The image displayed is on a red-orange-yellow-green scale, with red being least like the input' +
	' endmember and green being most like the input endmember.');



/*--------------------------------------------------------*/
// Addendum: Optional Context
/*--------------------------------------------------------*/

//addToMap(july_waypnts, {color:'ffffff'}, 'July Waypoints');
//addToMap(spurgetraining, {color:'ffffff'}, 'Spurge Training Regions');
addToMap(ranchimage.paint(ranchoutline, 0 /* color */ , 3 /* width */ ), {}, 'Ranch Outline');



/*--------------------------------------------------------*/
// Addendum: Optional Charts
/*--------------------------------------------------------*/

/*###### 
This section takes sample points from different landcover types around the study site and plots their spectral profiles onto 
a chart alongside the profiles of the pixel-and field-endmember options.  The charts only show bands used in the 
classification after the image has been spectrally subset.  In the future, these spectral profiles could be averages of 
several points or the average reflectance values of polygonal training regions either drawn in real time or uploaded from 
field-collected shapefiles.
######*/

// Sample Points Creation
var Invasive = ee.Feature(
	ee.Geometry.Point(-106.5184, 44.6161), {
		'label': 'Invasive'
	});
var Crop = ee.Feature(
	ee.Geometry.Point(-106.5083, 44.5747), {
		'label': 'Crop'
	});
var Road = ee.Feature(
	ee.Geometry.Point(-106.5060, 44.5777), {
		'label': 'Road'
	});
var Endmember = ee.Feature(pointcoordinates, {
	'label': 'Endmember'
});

//Cast the  Sample Points into a Feature Collection.
var ucrossPoints = ee.FeatureCollection([Invasive, Crop, Road, Endmember]);

var collectionList = [OriginalImage, OriginalImage];
var imageCollection = ee.ImageCollection.fromImages(collectionList);
var info = imageCollection.getRegion(ucrossPoints, 30);

//Format the data intolists in order to chart values pulled from the image with values from
//the endmember feature collection.
var bandList = ee.List(info.get(0)).slice(4);
//print('Band Names', bandList);

var wavelengthList = libendmemberprep1.aggregate_array('Wavelength');
var wavelengthArray = ee.Array(wavelengthList).repeat(1, 1);
//print('Wavelength List', wavelengthArray);

var invasiveList = ee.List(info.get(7)).slice(4);
var invasiveArray = ee.Array(invasiveList).repeat(1, 1);
//print('Invasive Reflectance', invasiveArray);

var cropList = ee.List(info.get(1)).slice(4);
var cropArray = ee.Array(cropList).repeat(1, 1);
//print('Crop Reflectance', cropArray);

var roadList = ee.List(info.get(3)).slice(4);
var roadArray = ee.Array(roadList).repeat(1, 1);
//print('Road Reflectance', roadArray);

var pixelEndList = ee.List(info.get(5)).slice(4);
var pixelEndArray = ee.Array(pixelEndList).repeat(1, 1);
//print('Road Reflectance', roadArray);

var endmemberList = libendmemberprep1.aggregate_array('Pixel_Value');
var endmemberArray = ee.Array(endmemberList).repeat(1, 1);
//print('Endmember Reflectance', endmemberArray);

//Concatenatethe lists of each cover type and cast them into an array.
var arraysconcat = ee.Array.cat([invasiveArray, cropArray, roadArray, pixelEndArray, endmemberArray], 1);
//print("Concatenated Arrays", arraysconcat);

//Chart the values.
var arrayChart = Chart.array.values(
	arraysconcat, 0, wavelengthArray).setSeriesNames(["Invasive", "Crop", "Road", "Image Endmember", "Endmember"]);
arrayChart = arrayChart.setOptions({
	title: 'EO-1 spectra at three points on Ucross Ranch',
	hAxis: {
		title: 'Wavelength'
	},
	vAxis: {
		title: 'Reflectance Value'
	},
	lineWidth: 1,
	pointSize: 4,
	series: {
		0: {
			color: 'darkgreen'
		},
		1: {
			color: 'lightgreen'
		},
		2: {
			color: 'lightblue'
		},
		3: {
			color: 'red'
		},
		4: {
			color: 'darkred'
		}
	}
});
print(arrayChart);

// ~~~
