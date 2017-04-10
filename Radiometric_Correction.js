// Radiometric Correction

// The radiometric correction function below corrects the inputted "yearximage"
// (i.e., any Landsat 5 / 7 / 8 image expecting to be bound by the "studyarea"
// geometry requested by the function) via an adapted "No Change Area Regression"
// protocol established by Elvidge and Yuan (1995). The "yearyimage" (identical
// in image type and preparation to the "yearximage" but from another date) is the 
// base image to which the "yearximage" is normalized / corrected.

// N.B. "ncpixelinclusion" indicates the number of pixels one would like to (approximately)
// include in the NC regression radiometric correction procedure. As well, images inputted
// to this function are assumed to have undergone the Tasseled Cap Transformation (with band
// names being "Brightness", "Greenness", "Wetness", "Fourth", "Fifth", "Sixth")

// For reference, see:
//http://www.researchgate.net/publication/279894684_Relative_radiometric_normalization_of_Landsat_multispectral_scanner_%28MSS%29_data_using_an_automatic_scattergram-controlled_regression}


var radiometriccorrection = function(yearximage, yearyimage, studyarea, ncpixelinclusion) {

	// Create a list variable for the Tasseled Cap Bands (selected from an image that's
	// already undergone the TC transformation)
	var tcbandlist = ["Brightness", "Greenness", "Wetness", "Fourth", "Fifth", "Sixth"];

	// Create a new bandlist for the radiometrically altered images
	var rbandlist = ["RB1", "RB2", "RB3", "RB4", "RB5", "RB6"];

	// Options to print and ensure 
	// print("Check 1",bandselectionlistx);
	// print("Check 2",bandselectionlisty);

	// Select the TC band list from the year y image and uniformly rename them
	var yearybandselected = ee.Image(yearyimage)
		.select(tcbandlist, rbandlist);
	// print(yearybandselected);

	// Select the TC band list from the year x image and uniformly rename them
	var yearxbandselected = ee.Image(yearximage)
		.select(tcbandlist, rbandlist);
	// print(yearxbandselected);




	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	// This subsection prepares the image variables for later mathematical use when
	// finalizing the radiometric correction

	// Compute the total change across all bands at each pixel; in particular, calculate
	// the root mean square of the change at each pixel in order to highlight pixels
	// with significant change across one or two bands
	var totalchangeperpixel = yearybandselected.subtract(yearxbandselected)
		.multiply(yearybandselected.subtract(yearxbandselected))
		.reduce('sum')
		.divide(ee.Image(6))
		.sqrt();
	// Map.addLayer(totalchangeperpixel,{},"Total Change");

	// Compute the minimum and the maximum change, using the variable above, across the
	// the study area
	var minmaxtotalchange = ee.Image(totalchangeperpixel).reduceRegion(ee.Reducer.minMax()
		.unweighted()
		.combine(ee.Reducer.count()
			.unweighted(), '', true), studyarea, 30);
	// print("MinMax", minmaxtotalchange);

	// Normalize the total change values across all pixels
	var normalizedchange = totalchangeperpixel.subtract(ee.Image(ee.Number(minmaxtotalchange.get(
			"sum_min"))))
		.divide(ee.Image(ee.Number(minmaxtotalchange.get("sum_max")))
			.subtract(ee.Image(ee.Number(minmaxtotalchange.get("sum_min")))));
	// Map.addLayer(normalizedchange,{},"Normalized Change Image");
	// print(normalizedchange.reduceRegion(ee.Reducer.minMax().unweighted()
	//                                 .combine(ee.Reducer.count().unweighted(), '', true)
	//                                 ,studyarea,30));

	// Run a reducer to compute the total number of pixels in the image
	var totalpixelcount = minmaxtotalchange.get("sum_count");
	// print(totalpixelcount);

	// Assign the number of no change pixels from the user input section into a variable
	// for use within the function
	var pixelinclusionlimit = ncpixelinclusion;

	// Convert the raw number of pixels defined as "no change" into a percentage of the 
	// total number of image pixels
	var percentagetoseek = ee.Number(pixelinclusionlimit)
		.divide(totalpixelcount)
		.multiply(100);
	// print(percentagetoseek);

	// Acquire the percentile equivalent of the normalized pixel value (with the percentile
	// sought determined the by user input)
	var normchangepercent = normalizedchange.reduceRegion(ee.Reducer.percentile([percentagetoseek]),
			studyarea, 30)
		.get("sum");
	// print(normchangepercent);

	// Create a mask of all pixels except those of interest to the radiometric correction (i.e., the
	// "no change area" pixels)
	var pixelsofinterestmask = ee.Image(normalizedchange.lte(ee.Image(ee.Number(normchangepercent))));
	// Map.addLayer(pixelsofinterestmask);
	// print(pixelsofinterestmask);

	// Apply the mask to the year y image
	var yearymasked = yearybandselected.mask(pixelsofinterestmask);
	// print(yearymasked.reduceRegion(ee.Reducer.minMax().unweighted()
	//                                 .combine(ee.Reducer.count().unweighted(), '', true)
	//                                 ,studyarea,30));

	// Apply the mask to the year x image
	var yearxmasked = yearxbandselected.mask(pixelsofinterestmask);
	// print(yearxmasked.reduceRegion(ee.Reducer.minMax().unweighted()
	//                                 .combine(ee.Reducer.count().unweighted(), '', true)
	//                                 ,studyarea,30));

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	// This subsection prepares the variables for the mathematical computation of the
	// "no change" regression process as described by Elvidge and Yuan (1996).

	// Compute the mean of the year y nochange pixels
	var yearymean = yearymasked.reduceRegion({
		reducer: ee.Reducer.mean(),
		geometry: studyarea,
		scale: 30
	});
	// print(yearymean);

	// Compute the mean, sample variance, and count of the year x no change pixels (i.e.,
	// the image that we're correcting)
	var collectedreduce = yearxmasked.reduceRegion({
		reducer: ee.Reducer.sampleVariance()
			.unweighted()
			.combine(ee.Reducer.mean()
				.unweighted(), '', true)
			.combine(ee.Reducer.count()
				.unweighted(), '', true),
		geometry: studyarea,
		scale: 30
	});
	// print("Collected Reducer Outputs",collectedreduce);

	// Create an array image of the year x mean values for all bands (i.e., a 1x6 array)
	var yearxmeanarray = collectedreduce.toArray(["RB1_mean", "RB2_mean", "RB3_mean", "RB4_mean",
		"RB5_mean", "RB6_mean"
	]);
	// print(yearxmeanarray);


	// Create an array image of the year x variance values for all bands (i.e., a 1x6 array)
	var yearxvariance = collectedreduce.toArray(["RB1_variance", "RB2_variance", "RB3_variance",
		"RB4_variance", "RB5_variance", "RB6_variance"
	]);
	// print(yearxvariance);


	// Acquire the number of pixels across the year x "no change" region
	var yearxpixelcount = collectedreduce.get("RB1_count");
	// print(yearxpixelcount);

	// Create an array image of the year y mean values
	var yearymeanarray = yearymean.toArray();
	// print(yearymeanarray);




	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	// This subsection includes two calls to format the year x and y "no change"
	// images into array images

	var yearxclippedarray = yearxmasked.toArray();
	// Map.addLayer(yearxclippedarray,{},"Year 2 Clipped Array Image");

	var yearyclippedarray = yearymasked.toArray();
	// Map.addLayer(yearyclippedarray,{},"Year 2 Clipped Array Image");




	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	// This section includes calls to compute the deviation values of both x and y
	// "no change" regions, then uses those values to compute the sample covariance
	// across the two regions/images.

	var yearxdeviation = yearxclippedarray.subtract(yearxmeanarray);
	// Map.addLayer(yearxdeviation,{},"Year 1 Deviation");

	var yearydeviation = yearyclippedarray.subtract(yearymeanarray);
	// Map.addLayer(yearydeviation,{},"Year 2 Deviation");

	var samplecovariance = ee.Array(yearxdeviation.multiply(yearydeviation)
		.arrayFlatten([
			['RB1', 'RB2', 'RB3', 'RB4', 'RB5', 'RB6']
		])
		.reduceRegion(ee.Reducer.sum(), studyarea, 30)
		.toArray()
		.divide(yearxpixelcount));
	// print("Sample Covariance", samplecovariance);




	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	// This subsection includes two calls that format the alpha and beta variables
	// according to mathematics presented Elvidge and Yuan (1996)

	var alpha = samplecovariance.divide(yearxvariance);
	// print("Alpha",alpha);

	var beta = yearymeanarray.subtract(alpha.multiply(yearxmeanarray));
	// print("Beta",beta);




	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	// This subsection finalizes the radiometric normalization using the alpha and beta
	// variables in the computation described by Elvidge and Yuan (1996). It includes the
	// final return call for the radiometric normalization function

	var normalizedcloudfreeimage = ee.Image(yearxbandselected.toArray()
		.multiply(alpha)
		.add(beta)
		.arrayFlatten([tcbandlist])
		.copyProperties(ee.Image(yearximage)));
	// print("Normalized Image",normalizedcloudfreeimage);

	return normalizedcloudfreeimage;
};
