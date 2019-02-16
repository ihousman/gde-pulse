// function run(){

//   function addTimeBand(img){
//   var d = ee.Date(img.get('system:time_start'));
//   var years = d.difference(ee.Date('0000-01-01'), 'year');
//   // var y = d.get('year');
//   // var p = d.getFraction('year');
//   y = ee.Image(years).updateMask(img.mask()).rename('year').float().copyProperties(img,['system:time_start'])

//   return img.addBands(y);
// }

//   // Map.addREST(standardTileURLFunction('https://landfire.cr.usgs.gov/arcgis/services/Landfire/US_140/MapServer/WCSServer/'),'Hillshade',false,15, 'Hillshade for draping layers over')
// var mtbs = ee.ImageCollection('users/ianhousman/MTBS/Collection');
// mtbs = mtbs.map(function(img){return img.select([0],['burnSeverity']).byte().updateMask(img.neq(0).and(img.neq(6)))});
// var mtbsForViz = mtbs.map(addTimeBand);

// mtbsForViz = mtbsForViz.qualityMosaic('year');
// Map.addLayer(mtbsForViz.select([0]),{'min':1,'max':6,'palette':'006400,7fffd4,ffff00,ff0000,7fff00,ffffff','opacity':1,'addToLegend':false},'MTBS 1984-2015 Severity Composite',false,null,null,'Burn severity  from Monitoring Trends in Burn Severity project (mtbs.gov). Most recent year burned is displayed. Use right-click charting to chart all years (1984-2015)')
// Map.addLayer(mtbsForViz.select([1]),{'min':1984,'max':2015,'palette':'9400D3,4B0082,00F,0F0,FF0,FF7F00,F00'},'MTBS 1984-2015 Year of Most Recent Fire',false,'','FFF','Year of most recent fire from Monitoring Trends in Burn Severity project (mtbs.gov). Most recent year burned is displayed. Use right-click charting to chart all years (1984-2015)')


  
// }
function run(){

function simpleAddIndices(in_image){
    in_image = in_image.addBands(in_image.normalizedDifference(['nir', 'red']).select([0],['NDVI']));
    in_image = in_image.addBands(in_image.normalizedDifference(['nir', 'swir2']).select([0],['NBR']));
    in_image = in_image.addBands(in_image.normalizedDifference(['nir', 'swir1']).select([0],['NDMI']));
    in_image = in_image.addBands(in_image.normalizedDifference(['green', 'swir1']).select([0],['NDSI']));
  
    return in_image;
}
///////////////////////////////////////////////////////////////////////////////
// Function to compute the Tasseled Cap transformation and return an image
// with the following bands added: ['brightness', 'greenness', 'wetness', 
// 'fourth', 'fifth', 'sixth']
function getTasseledCap(image) {
 
  var bands = ee.List(['blue','green','red','nir','swir1','swir2']);
  // // Kauth-Thomas coefficients for Thematic Mapper data
  // var coefficients = ee.Array([
  //   [0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
  //   [-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
  //   [0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
  //   [-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
  //   [-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
  //   [0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]
  // ]);
  
  //Crist 1985 coeffs - TOA refl (http://www.gis.usu.edu/~doug/RS5750/assign/OLD/RSE(17)-301.pdf)
  var coefficients = ee.Array([[0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303],
                    [-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446],
                    [0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109]
                   //  [-0.2117, -0.0284, 0.1302, -0.1007, 0.6529, -0.7078],
                   //  [-0.8669, -0.1835, 0.3856, 0.0408, -0.1132, 0.2272],
                   // [0.3677, -0.8200, 0.4354, 0.0518, -0.0066, -0.0104]
                   ]);
  // Make an Array Image, with a 1-D Array per pixel.
  var arrayImage1D = image.select(bands).toArray();
  
  // Make an Array Image with a 2-D Array per pixel, 6x1.
  var arrayImage2D = arrayImage1D.toArray(1);
  
  var componentsImage = ee.Image(coefficients)
    .matrixMultiply(arrayImage2D)
    // Get rid of the extra dimensions.
    .arrayProject([0])
    // Get a multi-band image with TC-named bands.
    .arrayFlatten(
      [['brightness', 'greenness', 'wetness']])
    .float();
  
  return image.addBands(componentsImage);
}

///////////////////////////////////////////////////////////////////////////////
// Function to add Tasseled Cap angles and distances to an image.
// Assumes image has bands: 'brightness', 'greenness', and 'wetness'.

function simpleAddTCAngles(image){
  // Select brightness, greenness, and wetness bands
  var brightness = image.select(['brightness']);
  var greenness = image.select(['greenness']);
  var wetness = image.select(['wetness']);
  
  // Calculate Tasseled Cap angles and distances
  var tcAngleBG = brightness.atan2(greenness).divide(Math.PI).rename('tcAngleBG');
  // var tcAngleGW = greenness.atan2(wetness).divide(Math.PI).rename('tcAngleGW');
  // var tcAngleBW = brightness.atan2(wetness).divide(Math.PI).rename('tcAngleBW');
  // var tcDistBG = brightness.hypot(greenness).rename('tcDistBG');
  // var tcDistGW = greenness.hypot(wetness).rename('tcDistGW');
  // var tcDistBW = brightness.hypot(wetness).rename('tcDistBW');
  image = image.addBands(tcAngleBG)
  // .addBands(tcAngleGW)
  //   .addBands(tcAngleBW).addBands(tcDistBG).addBands(tcDistGW)
  //   .addBands(tcDistBW);
  return image;
}
var addYear = function(image) {
  var t = image.get('system:time_start');
  var y = ee.Date(t).get('year');
  var yimg = ee.Image(y).short().rename('year');
  var addyimg = image.addBands(yimg);//yimg.addBands(image)
  return addyimg.float() ;
};
/////////////////////////////////////////////////
//Helper function to join two collections- Source: code.earthengine.google.com
    function joinCollections(c1,c2){
      var MergeBands = function(element) {
        // A function to merge the bands together.
        // After a join, results are in 'primary' and 'secondary' properties.
        return ee.Image.cat(element.get('primary'), element.get('secondary'));
      };

      var join = ee.Join.inner();
      var filter = ee.Filter.equals('system:time_start', null, 'system:time_start');
      var joined = ee.ImageCollection(join.apply(c1, c2, filter));
     
      joined = ee.ImageCollection(joined.map(MergeBands));
      // joined = joined.map(function(img){return img.mask(img.mask().and(img.reduce(ee.Reducer.min()).neq(0)))});
      return joined;
    }
///////////////////////////////////////////
//Function to  add SAVI and EVI
function addSAVIandEVI(img){
  // Add Enhanced Vegetation Index (EVI)
  var evi = img.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': img.select('nir'),
      'RED': img.select('red'),
      'BLUE': img.select('blue')
  }).float();
  img = img.addBands(evi.rename('EVI'));
  
  // Add Soil Adjust Vegetation Index (SAVI)
  // using L = 0.5;
  var savi = img.expression(
    '(NIR - RED) * (1 + 0.5)/(NIR + RED + 0.5)', {
      'NIR': img.select('nir'),
      'RED': img.select('red')
  }).float();
  
  
  ////////////////////////////////////////////////////////////////////////////////
  //NIRv: Badgley, G., Field, C. B., & Berry, J. A. (2017). Canopy near-infrared reflectance and terrestrial photosynthesis. Science Advances, 3, e1602244.
  //https://www.researchgate.net/publication/315534107_Canopy_near-infrared_reflectance_and_terrestrial_photosynthesis
  // NIRv function: ‘image’ is a 2 band stack of NDVI and NIR
  //////////////////////////////////////////////////////////////////////////////////////////
  var NIRv =  img.select(['NDVI']).subtract(0.08)
              .multiply(img.select(['nir']));//.multiply(0.0001))

  img = img.addBands(savi.rename('SAVI')).addBands(NIRv.rename('NIRv'));
  return img;
}
//////////////////////////////////////////////////////
//////////
//Helper to multiply image
function multBands(img,distDir,by){
    var out = img.multiply(ee.Image(distDir).multiply(by));
    out  = out.copyProperties(img,['system:time_start']);
    return out;
  }
function addToImage(img,howMuch){
  img = ee.Image(img);
    var out = img.add(ee.Image(howMuch));
    out  = out.copyProperties(img,['system:time_start'])
              .copyProperties(img);
    return out;
  }
//Function for zero padding
//Taken from: https://stackoverflow.com/questions/10073699/pad-a-number-with-leading-zeros-in-javascript
function pad(n, width, z) {
  z = z || '0';
  n = n + '';
  return n.length >= width ? n : new Array(width - n.length + 1).join(z) + n;
}

//Bring in igdes
// var f = ee.FeatureCollection('users/Shree1175/iGDE_5_2018_V1_joined_ndvi_annDepth')
var f = ee.FeatureCollection('projects/igde-work/igde-data/GDEpulse2018_iGDE_V1_20180802_joined_annual_depth_macro_veg');

var fRast = ee.Image().paint(f,5,1);
var fRandom = f.reduceToImage(['ORIG_FID'],ee.Reducer.first());
fRandom = fRandom.randomVisualizer();

// fRast = fRast.where(fRast.mask(),1)

//Set up the years to filter on- this is hard-coded since its set up oddly
var years = ee.List.sequence(1985,2018);


//Reformat the igdes to have a unique feature per year
var igdeyr = years.getInfo().map(function(yz){
  var fieldName ='Depth'+ yz.toString();
  // var t = f.select([fieldName], ['AvgAnnD'])
  //         .map(function(ft){return ft.set('year',yz)});
  var t = f.select([fieldName], ['AvgAnnD']);
  var depth = t.reduceToImage(['AvgAnnD'], ee.Reducer.first());
  var tID = f.select(['ORIG_FID']).reduceToImage(['ORIG_FID'], ee.Reducer.first());
  t = depth.multiply(-1);
  t = t.updateMask(t.select([0]).gt(-1000))
      // .divide(100)
      .addBands(tID.int64())
      .rename(['Depth-To-Groundwater','ORIG_FID'])
      .set('system:time_start',ee.Date.fromYMD(yz,6,1).millis())
  return t;
});
igdeyr = ee.ImageCollection(igdeyr);
// //Various preset min/max stretches for data visualization
// var vizParamsCO1 = {'min': 0.05,'max': [0.3,0.6,0.35],   'bands':'swir1,nir,red'};
var vizParamsCO2 = {'min': 0.15,'max': [0.35,0.8,0.4],   'bands':'swir1,nir,red', 'gamma': 1.6};
// var vizParamsCO3 = {'min': 0.05,'max': [0.3,0.4,0.4],   'bands':'swir1,nir,red', 'gamma':1.6};
// var vizParamsFalse = {'min': 0.1,'max': [0.3,0.3,0.3],   'bands':'nir,swir1,red'};
// var vizParamsViz = {'min': 0.05, 'max': 0.3,'bands': 'red,green,blue', 'gamma': 1.6};
var vizParams = vizParamsCO2;
var startYear = 1985;
var endYear = 2018;
var startJulian = 190;
var endJulian = 250;
var timeBuffer = 2;
var weights = [1,1,5,1,1];
var compositingMethod = 'medoid'
var ts = ee.ImageCollection('projects/igde-work/raster-data/composite-collection');
var ts = ts
        .map(function(img){return multBands(img,1,0.0001)})
        .map(simpleAddIndices)
        // .map(getTasseledCap)
        // .map(simpleAddTCAngles)
        // .map(addSAVIandEVI)
        // .map(function(img){return multBands(img,1,100)})
//         .map(addYear)
var lt = ee.ImageCollection('projects/igde-work/raster-data/LANDTRENDR-collection')
        .select(['.*_NDMI'])//,'.*_SAVI','.*_EVI'])
        .map(function(img){
          img = ee.Image(img);
          var out =  multBands(img,1,0.0001);
          out = addToImage(out,-1);
          out = ee.Image(out);
          out =  multBands(out,1,20)
          return out;
        });

var first = ee.Image(ts.first())//.multiply(0.01);
var last = ee.Image(ts.sort('system:time_start',false).first())//.multiply(0.01);
Map2.addLayer(first,vizParams,'Landsat Medoid Composite-'+startYear.toString(),false,null,null,'Medoid composite of Landsat data from 1985. Fmask was used for cloud and cloud shadow masking');
Map2.addLayer(last,vizParams,'Landsat Medoid Composite-'+endYear.toString(),false,null,null,'Medoid composite of Landsat data from 2018. Fmask was used for cloud and cloud shadow masking');
// var joined = igdeyr;
// var joined = joinCollections(igdeyr,ts.select(['NDVI','NBR','SAVI','EVI'])) ;
var tsForJoin = ts.select(['NDMI']).map(function(img){
  var img = ee.Image(img);
  var out = addToImage(img,-1);
  out = ee.Image(out);
  out = multBands(out,1,20);
  // out = addToImage(out,-10);
  return out
});

var daymet = ee.ImageCollection('projects/igde-work/raster-data/DAYMET-Collection').select([0],['Annual Daily Average Precip mm'])
            //   .map(function(img){
            //   var img = ee.Image(img);
            //   var out = addToImage(img,-10);
            //   out = ee.Image(out);
            //   // out = multBands(out,1,20);
            //   // out = addToImage(out,-10);
            //   return out
            // })
var joined = joinCollections(igdeyr,tsForJoin);
joined = joinCollections(joined,lt);
joined = joinCollections(joined,daymet);
joined = joined.map(function(img){

  var out = img.reduceConnectedComponents(ee.Reducer.mean(), 'ORIG_FID', 256);
  img = img.select([0,2,3,4])
  // out = out.addBands(img.select([0,1,2]))
  out = ee.Image(out.copyProperties(img,['system:time_start'])).updateMask(img.mask());
  
  return out;
  
})
var depthCount = igdeyr.select([0]).count();
Map2.addLayer(depthCount,{'min':1,'max':20,'palette':'AAA,00F'},'Depth to groundwater years with observation count ('+startYear.toString() + '-'+ endYear.toString()+')',false,'observations',null,'Number of years from 1985-2018 with groundwater data');

Map2.addLayer(igdeyr.select([0]).max(),{'min':-50,'max':0,'palette':'800,080'},'Max Depth to groundwater ('+startYear.toString() + '-'+ endYear.toString()+')',false,'m',null,'Max depth to ground water from 1985-2018 for closest well data to given igde');
Map2.addLayer(igdeyr.select([0]).min(),{'min':-50,'max':0,'palette':'800,080'},'Min Depth to groundwater ('+startYear.toString() + '-'+ endYear.toString()+')',false,'m',null,'Min depth to ground water from 1985-2018 for closest well data to given igde');
Map2.addLayer(igdeyr.select([0]).median(),{'min':-50,'max':0,'palette':'800,080'},'Median Depth to groundwater ('+startYear.toString() + '-'+ endYear.toString()+')',false,'m',null,'Median depth to ground water from 1985-2018 for closest well data to given igde');
Map2.addLayer(fRandom,{opacity:0.8,addToLegend:false},'iGDE Random Color by ID',true,null,null,'Random color assigned by unique ID to visualize what area is being considered an individual iGDE');
Map2.addLayer(fRast,{'min':1,'max':1,'palette':'05F','opacity':0.8,addToLegend:false},'iGDEs',true);
chartCollection = joined;


// ////////////////////////////////////////////////////////////////
// //Function to create list of numbers
// function range(start, stop, step){
//   // start = parseInt(start);
//   // stop = parseInt(stop);
//     if (typeof stop=='undefined'){
//         // one param defined
//         stop = start;
//         start = 0;
//     };
//     if (typeof step=='undefined'){
//         step = 1;
//     };
//     if ((step>0 && start>=stop) || (step<0 && start<=stop)){
//         return [];
//     };
//     var result = [];
//     for (var i=start; step>0 ? i<stop : i>stop; i+=step){
//         result.push(i);
//     };
//     return result;
// };
// ////////////////////////////////////////////////////////////////////////


// var nlcd92 = ee.Image('USGS/NLCD/NLCD1992');
// var nlcd01 = ee.Image('USGS/NLCD/NLCD2001');
// var nlcd06 = ee.Image('USGS/NLCD/NLCD2006');
// var nlcd11 = ee.Image('USGS/NLCD/NLCD2011');

// var cdlList = [
// '1997',
// '1998',
// '1999',
// '2000',
// '2001',
// '2002',
// '2003',
// '2004',
// '2005',
// '2006',
// '2007',
// '2008',
// '2009',
// '2010',
// '2011',
// '2012',
// '2013',
// '2014',
// '2015',
// '2016']

// // cdlList.slice(-1).map(function(c){

// //   if(c === '2005' || c === '2007'){
// //        var cdlA = ee.Image('USDA/NASS/CDL/'+c + 'a').select([0]);
  
// //       var cdlB = ee.Image('USDA/NASS/CDL/'+c + 'b').select([0]);
     
// //       var cdl = cdlA.mask(cdlA.mask().or(cdlB.mask()));//.add(cdl).subtract(cdl)
// //       cdl = cdl.where(cdlB.mask().and(cdlA.mask().not()),cdlB)
      
// //       Map.addLayer(cdl,{'addToLegend' : false},'CDL ' +c,false)
// //   }
// //   else{Map.addLayer(ee.Image('USDA/NASS/CDL/'+c).select([0]),{'addToLegend' : false},'CDL ' +c,false,null,null,'Latest CDL Layer for data visualization.  Right-click charting will plot all available CDL layers.')}
  
// // })

// var cdl = ee.ImageCollection('USDA/NASS/CDL').select([0],['cdl']);



// // var years = range(1984,2015);//[1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1997,1998];

// // years.slice(0,10).map(function(yr){
// //   var mtbs = ee.Image('users/ianhousman/MTBS/mtbs_CONUS_dt_'+yr.toString()+'_20160401')
// //   var mtbsEdge = mtbs.gt(0).focal_max(3).subtract(mtbs.gt(0).focal_min(3));
// // // Map.addLayer(mtbsEdge.mask(mtbs),{'min':0,'max':2,'palette':'500,F00','opacity':0.5,'addToLegend':false},'MTBS '+yr.toString(),false,null,null,'Burn areas > 1000 acres west of the Mississippi and > 500 acres for the year '+yr.toString()+' east of the Mississippi (mtbs.gov)') 
// // Map.addLayer(mtbs.mask(mtbs),{'min':1,'max':6,'palette':'#006400,#7fffd4,#ffff00,#ff0000,#7fff00,#ffffff','opacity':1,'addToLegend':false},'MTBS '+yr.toString(),false,null,null,'Burn areas > 1000 acres west of the Mississippi and > 500 acres for the year '+yr.toString()+' east of the Mississippi (mtbs.gov)') 

// // })



// function addTimeBand(img){
//   var d = ee.Date(img.get('system:time_start'));
//   var years = d.difference(ee.Date('0000-01-01'), 'year');
//   // var y = d.get('year');
//   // var p = d.getFraction('year');
//   y = ee.Image(years).updateMask(img.mask()).rename('year').float().copyProperties(img,['system:time_start'])

//   return img.addBands(y);
// }

// // Map.addREST(standardTileURLFunction('http://server.arcgisonline.com/arcgis/rest/services/Elevation/World_Hillshade/MapServer/tile/'),'Hillshade',false,15, 'Hillshade for draping layers over')

// Map.addREST(function(coord, zoom) {
//                     var tilesPerGlobe = 1 << zoom;
//                     var x = coord.x % tilesPerGlobe;
//                     if (x < 0) {x = tilesPerGlobe+x;}
//                     return "http://api.tiles.mapbox.com/v4/digitalglobe.nako1fhg/" + zoom + "/" + x + "/" + coord.y + ".png?access_token=pk.eyJ1IjoiZGlnaXRhbGdsb2JlIiwiYSI6ImNpczBqYml6bzAzb2Yyb3A2YWI1aWVvc2sifQ.ZWTVvUKerHzA6AX41LS9LQ";
//                 },'Digital Globe Terrain',false,20,'Digital Globe terrain data with vector data')
// // Map.addREST(function(coord, zoom) {
// //                     var tilesPerGlobe = 1 << zoom;
// //                     var x = coord.x % tilesPerGlobe;
// //                     if (x < 0) {x = tilesPerGlobe+x;}
// //                     return "http://b-globalheat.strava.com/tiles/running/color7/" + zoom + "/" + x + "/" + coord.y + ".png?v=6";
// //                 },'Strava Running',false,20,'Digital Globe terrain data with vector data')
// // Map.addREST(function(coord, zoom) {
// //                     var tilesPerGlobe = 1 << zoom;
// //                     var x = coord.x % tilesPerGlobe;
// //                     if (x < 0) {x = tilesPerGlobe+x;}
// //                     return "https://tiles.planet.com/data/v1/PSScene3Band/20161221_024131_0e19/" + zoom + "/" + x + "/" + coord.y + ".png?api_key={0ba0756f4e694818aa8eb32e5f933f21}";
// //                 },'Planet Labs',false,20,'Digital Globe terrain data with vector data')

// // Map.addREST(function(coord, zoom) {
// //                     var tilesPerGlobe = 1 << zoom;
// //                     var x = coord.x % tilesPerGlobe;
// //                     if (x < 0) {x = tilesPerGlobe+x;}
// //                     return "http://b-globalheat.strava.com/tiles/cycling/color7/" + zoom + "/" + x + "/" + coord.y + ".png?v=6";
// //                 },'Strava Biking',false,20,'Digital Globe terrain data with vector data')

// // http://b-globalheat.strava.com/tiles/running/color7/9/265/180.png?v=6
// // Map.addREST(standardTileURLFunction('http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/'),'ESRI World Imagery',false,19, 'ESRI worldwide high resolution imagery')

// // Map.addREST(function(coord, zoom) {
// //                     var tilesPerGlobe = 1 << zoom;
// //                     var x = coord.x % tilesPerGlobe;
// //                     if (x < 0) {x = tilesPerGlobe+x;}
// //                     return "http://d.tiles.mapbox.com/v4/mapquest.satellite/" + zoom + "/" + x + "/" + coord.y + ".png?access_token=pk.eyJ1IjoibWFwcXVlc3QiLCJhIjoiY2Q2N2RlMmNhY2NiZTRkMzlmZjJmZDk0NWU0ZGJlNTMifQ.mPRiEubbajc6a5y9ISgydg";
// //                 },'MapBox Hybrid',false,19,'Mapbox high resolution satellite data with vector data')


// // Map.addREST(function(coord, zoom) {
// //                     var tilesPerGlobe = 1 << zoom;
// //                     var x = coord.x % tilesPerGlobe;
// //                     if (x < 0) {x = tilesPerGlobe+x;}
// //                     return "http://api.tiles.mapbox.com/v4/digitalglobe.nal0mpda/" + zoom + "/" + x + "/" + coord.y + ".png?access_token=pk.eyJ1IjoiZGlnaXRhbGdsb2JlIiwiYSI6ImNpczBqYml6bzAzb2Yyb3A2YWI1aWVvc2sifQ.ZWTVvUKerHzA6AX41LS9LQ";
// //                 },'Digital Globe Recent Imagery',false,19,'Digital Globe recent high resolution satellite data with vector data')

// // Map.addREST(function(coord, zoom) {
// //                     var tilesPerGlobe = 1 << zoom;
// //                     var x = coord.x % tilesPerGlobe;
// //                     if (x < 0) {x = tilesPerGlobe+x;}
// //                     return "https://api.mapbox.com/styles/v1/digitalglobe/civvbr6nf00092jruyo1dkw0o/tiles/256/" + zoom + "/" + x + "/" + coord.y + "?access_token=pk.eyJ1IjoiZGlnaXRhbGdsb2JlIiwiYSI6IjVlZjc2MjUxMDA1NzE1Y2M1MjM5NzU1MWUxYWRmNGU5In0.vHPkCmZRekCR4LQqhvN6oA";
// //                 },'Digital Globe Recent Imagery Year',false,14,'Digital Globe recent high resolution satellite data year. Relevant for images at zoom levels >= 13.  May need to zoom out to see this layer')

// // Map.addREST(function(coord, zoom) {
// //                     var tilesPerGlobe = 1 << zoom;
// //                     var x = coord.x % tilesPerGlobe;
// //                     if (x < 0) {x = tilesPerGlobe+x;}
// //                     return "http://thermal.kk7.ch/php/tile.php?typ=skyways&t=all&z=" + zoom + "&x=" + x+100 + "&y=" + coord.y +100+ ".png";
// //                 },'thermal',false,19,'Digital Globe recent high resolution satellite data with vector data')




// // Map.addREST(function(coord, zoom) {
// //                     var tilesPerGlobe = 1 << zoom;
// //                     var x = coord.x % tilesPerGlobe;
// //                     if (x < 0) {x = tilesPerGlobe+x;}
// //                     return "https://api.tiles.mapbox.com/v4/digitalglobe.p9l3p6p2/" + zoom + "/" + x + "/" + coord.y + ".png?access_token=pk.eyJ1IjoiZGlnaXRhbGdsb2JlIiwiYSI6ImNpczBqYml6bzAzb2Yyb3A2YWI1aWVvc2sifQ.ZWTVvUKerHzA6AX41LS9LQ";
// //                 },'Digital Globe Composite Satellite',false,19,'Digital Globe high resolution satellite data ')

// Map.addREST(function(coord, zoom) {
//                         var quadCode = tileXYZToQuadKey(coord.x, coord.y, zoom)
//                         var center = map.getCenter()
//                         var url = 'http://ecn.t1.tiles.virtualearth.net/tiles/h'+quadCode +'.jpeg?g=5256&mkt=en-US&token=AuuPPqhBI42-Eo2fJ3uRFFrvhzM93YmEIZFbOtqjFPNW2V641ZvU4g53ivnlVT8b'
//                         var metadata = 'http://dev.virtualearth.net/REST/v1/Imagery/BasicMetadata/AerialWithLabels/'+center.lat()+','+center.lng()+'?orientation=0&zoomLevel='+zoom+'&include=ImageryProviders&key=AuuPPqhBI42-Eo2fJ3uRFFrvhzM93YmEIZFbOtqjFPNW2V641ZvU4g53ivnlVT8b'
//                         console.log(coord);
//                         console.log(metadata)
//                         $.getJSON(metadata, function(data) {
//                             var dates = data.resourceSets[0].resources[0];
//                             var startDate = dates.vintageStart;
//                             var endDate = dates.vintageEnd;
//                             console.log(dates);
//                             console.log([startDate,endDate])
//                         });
                        
//                     return url
//                 },'Bing Hybrid',false,19,'Bing high resolution satellite data with vector data')

// Map.addREST(function(coord, zoom) {
//                         var quadCode = tileXYZToQuadKey(coord.x, coord.y, zoom)
//                         // var url = 'http://t0.tiles.virtualearth.net/tiles/a'+quadCode +'.jpeg?g=854&mkt=en-US&token=AuuPPqhBI42-Eo2fJ3uRFFrvhzM93YmEIZFbOtqjFPNW2V641ZvU4g53ivnlVT8b'
//                         var url = 'http://mvexel.dev.openstreetmap.org/bingimageanalyzer/tile.php?t='+quadCode+'&nodepth='
                    
//                     return url
//                 },'Bing Hybrid Image Date',false,19,'Bing image acquisition date for high resolution satellite data with vector data')

// var nlcd = ee.ImageCollection([     nlcd92.addBands(ee.Image(1).mask(ee.Image(0))).set('system:time_start',ee.Date.fromYMD(1992,6,1).millis()),
//                                     nlcd01.set('system:time_start',ee.Date.fromYMD(2001,6,1).millis()),
//                                     nlcd06.set('system:time_start',ee.Date.fromYMD(2006,6,1).millis()),
//                                     nlcd11.set('system:time_start',ee.Date.fromYMD(2011,6,1).millis())])
//                                   .select([0,1],['landcover','impervious']);
// Map.addREST(standardTileURLFunction('http://server.arcgisonline.com/arcgis/rest/services/Specialty/Soil_Survey_Map/MapServer/tile/'),'SSURGO Soils',false,16, 'SSURGO soils data from USDA NRCS (http://server.arcgisonline.com/arcgis/rest/services/Specialty/Soil_Survey_Map/MapServer)')

// var globCover = ee.Image('ESA/GLOBCOVER_L4_200901_200912_V2_3').select([0]);
// Map.addLayer(globCover,{'addToLegend':false},'ESA GlobCover',false,null,null,'ESA global landcover based on ENVISAT at 300m res')

// globCover = ee.ImageCollection([globCover.rename(['GlobCover']).set('system:time_start',ee.Date.fromYMD(2009,6,1))]);
// ///////////////////////////////////////////////////
// //Add NAIP to viewer
// var naipYears = ee.List.sequence(2007,2017).getInfo();
// var naip = ee.ImageCollection("USDA/NAIP/DOQQ").select([0,1,2],['R','G','B']);
// naipYears.map(function(yr){
//   var naipT = naip.filter(ee.Filter.calendarRange(yr,yr,'year'));
//   Map.addLayer(naipT.mosaic(),{'min':25,'max':225,'addToLegend':false},'NAIP ' + yr.toString(),false,'','FFF','The National Agriculture Imagery Program (NAIP) acquired aerial imagery during the '+yr.toString()+' agricultural growing season in the continental U.S.');
// });

// ///////////////////////////////////////////
// var usForestType = ee.Image('users/ianhousman/conus_foresttype1').rename('FIA Forest Type');
// Map.addLayer(usForestType.updateMask(usForestType),{'min':0,'max':995,
//   'palette':'000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,000000,020202,020202,020202,020202,020202,020202,020202,050505,050505,050505,050505,050505,050505,050505,050505,050505,050505,050505,050505,070707,070707,070707,070707,070707,070707,070707,070707,0a0a0a,0a0a0a,0a0a0a,0a0a0a,0a0a0a,0a0a0a,0a0a0a,0a0a0a,0a0a0a,0a0a0a,0a0a0a,0a0a0a,0c0c0c,0c0c0c,0c0c0c,0c0c0c,0c0c0c,0c0c0c,0c0c0c,0c0c0c,0c0c0c,0c0c0c,0c0c0c,0f0f0f,0f0f0f,0f0f0f,0f0f0f,0f0f0f,0f0f0f,0f0f0f,0f0f0f,111111,111111,111111,111111,111111,111111,111111,111111,111111,111111,111111,ccaf89,ccaa87,cca884,cca582,cca082,141414,141414,141414,141414,161616,161616,161616,161616,161616,161616,161616,161616,161616,161616,161616,cc9b7f,cc987c,cc937c,cc8e7a,cc8977,cc8275,cc7c75,191919,1c1c1c,1c1c1c,1c1c1c,1c1c1c,1c1c1c,1c1c1c,1c1c1c,1c1c1c,1c1c1c,1c1c1c,1c1c1c,1c1c1c,cc7772,cc7070,1e1e1e,1e1e1e,1e1e1e,1e1e1e,1e1e1e,212121,212121,212121,212121,212121,212121,212121,212121,212121,212121,212121,212121,212121,cc7075,cc6d77,cc6b7a,cc687c,cc6882,cc6684,cc6689,cc608c,232323,232323,232323,262626,262626,262626,262626,262626,262626,262626,262626,282828,cc6091,cc5e93,cc5b98,cc599b,cc59a0,282828,282828,282828,282828,282828,282828,282828,282828,282828,2b2b2b,2b2b2b,2b2b2b,2b2b2b,2b2b2b,2b2b2b,cc56a3,cc54a8,2d2d2d,2d2d2d,2d2d2d,2d2d2d,2d2d2d,2d2d2d,2d2d2d,2d2d2d,303030,303030,303030,303030,303030,303030,303030,303030,333333,333333,cc51aa,333333,cc51af,cc4fb2,333333,333333,333333,333333,333333,333333,333333,333333,333333,353535,353535,353535,353535,353535,353535,353535,cc4cb7,383838,383838,383838,383838,383838,383838,383838,383838,3a3a3a,3a3a3a,3a3a3a,3a3a3a,3a3a3a,3a3a3a,3a3a3a,3a3a3a,3d3d3d,3d3d3d,3d3d3d,cc49ba,cc49bc,cc44c1,cc44c1,cc42c6,cc3fc9,cc3dc9,cc3dcc,cc3acc,cc38cc,cc35cc,3d3d3d,3f3f3f,3f3f3f,3f3f3f,3f3f3f,3f3f3f,3f3f3f,3f3f3f,3f3f3f,cc33cc,3f3f3f,3f3f3f,3f3f3f,424242,424242,424242,424242,444444,444444,444444,444444,444444,444444,444444,444444,444444,444444,444444,444444,cc33cc,444444,444444,cc30cc,c92dcc,474747,474747,494949,494949,494949,494949,494949,494949,494949,494949,494949,494949,494949,494949,4c4c4c,c92bce,4c4c4c,4c4c4c,4c4c4c,4c4c4c,4c4c4c,4c4c4c,4c4c4c,4c4c4c,4c4c4c,4c4c4c,4f4f4f,4f4f4f,4f4f4f,4f4f4f,4f4f4f,4f4f4f,4f4f4f,515151,515151,c628ce,c128ce,515151,515151,515151,515151,515151,515151,515151,515151,545454,545454,545454,545454,545454,545454,545454,545454,545454,545454,bf28cc,545454,565656,565656,ba26cc,b726cc,b223cc,ad23c9,595959,595959,a823c6,595959,595959,595959,595959,595959,595959,595959,595959,595959,a321c6,5b5b5b,9b21c4,9321c4,5b5b5b,5b5b5b,5b5b5b,5b5b5b,5b5b5b,5b5b5b,5b5b5b,5b5b5b,5b5b5b,5e5e5e,5e5e5e,5e5e5e,5e5e5e,606060,606060,606060,8e21c1,8721c1,7f1ec1,771ec1,701cbf,661cbc,5b1cbc,606060,511cba,606060,606060,606060,636363,636363,636363,636363,666666,666666,666666,666666,666666,666666,666666,666666,666666,666666,666666,666666,686868,686868,686868,686868,686868,686868,686868,686868,686868,686868,686868,686868,6b6b6b,6b6b6b,6b6b6b,6b6b6b,6b6b6b,6b6b6b,6b6b6b,6d6d6d,6d6d6d,6d6d6d,6d6d6d,6d6d6d,6d6d6d,6d6d6d,6d6d6d,6d6d6d,6d6d6d,6d6d6d,6d6d6d,707070,707070,707070,707070,707070,707070,707070,707070,707070,707070,707070,707070,727272,727272,727272,727272,727272,727272,727272,757575,757575,757575,757575,757575,757575,757575,757575,757575,757575,757575,757575,777777,777777,777777,777777,777777,777777,777777,777777,777777,777777,4719ba,3d19ba,3019ba,2316b5,1616b5,1621b5,162db2,143ab2,1444af,144caf,1156ad,1160aa,1168aa,1170aa,1175aa,7c7c7c,7c7c7c,7c7c7c,0f7ca8,0f82a5,7c7c7c,7f7f7f,7f7f7f,7f7f7f,7f7f7f,828282,828282,828282,828282,828282,828282,828282,828282,828282,828282,828282,828282,828282,828282,828282,828282,848484,848484,848484,848484,848484,848484,848484,848484,878787,878787,878787,878787,878787,878787,878787,898989,898989,898989,898989,898989,898989,898989,898989,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8c8c8c,8e8e8e,8e8e8e,8e8e8e,8e8e8e,8e8e8e,8e8e8e,8e8e8e,8e8e8e,919191,919191,919191,919191,919191,919191,919191,939393,939393,939393,939393,939393,0f87a5,0c8ca3,939393,939393,0c8ca3,0c91a3,0c93a0,0c939e,939393,939393,939393,969696,969696,969696,969696,989898,989898,989898,989898,989898,989898,989898,989898,989898,989898,989898,9b9b9b,9b9b9b,9b9b9b,9b9b9b,9b9b9b,9b9b9b,9b9b9b,9b9b9b,9b9b9b,9b9b9b,9b9b9b,9b9b9b,9e9e9e,9e9e9e,9e9e9e,9e9e9e,9e9e9e,9e9e9e,9e9e9e,9e9e9e,a0a0a0,a0a0a0,a0a0a0,a0a0a0,a0a0a0,a0a0a0,a0a0a0,a0a0a0,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a3a3a3,a5a5a5,a5a5a5,a5a5a5,a5a5a5,a5a5a5,a5a5a5,a5a5a5,a5a5a5,a8a8a8,a8a8a8,a8a8a8,a8a8a8,a8a8a8,a8a8a8,a8a8a8,a8a8a8,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,aaaaaa,0c939e,0a969b,0a969b,0a9698,0a9698,0a9396,0a9393,079393,079193,afafaf,afafaf,afafaf,afafaf,afafaf,afafaf,afafaf,b2b2b2,b2b2b2,b2b2b2,b2b2b2,b2b2b2,079191,b2b2b2,b2b2b2,b2b2b2,b2b2b2,b2b2b2,b2b2b2,b5b5b5,b5b5b5,b5b5b5,b5b5b5,b5b5b5,b5b5b5,b5b5b5,b5b5b5,b5b5b5,b5b5b5,b5b5b5,b7b7b7,b7b7b7,b7b7b7,b7b7b7,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bababa,bcbcbc,bcbcbc,bcbcbc,bcbcbc,bcbcbc,bcbcbc,bcbcbc,bcbcbc,bfbfbf,bfbfbf,bfbfbf,bfbfbf,bfbfbf,bfbfbf,bfbfbf,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c1c1c1,c4c4c4,c4c4c4,c4c4c4,c4c4c4,c4c4c4,c4c4c4,c4c4c4,c4c4c4,c6c6c6,c6c6c6,078e8e,058c8c,058c8c,c6c6c6,058c8c,c6c6c6,058989,c9c9c9,058787,c9c9c9,c9c9c9,c9c9c9,c9c9c9,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cccccc,cecece,cecece,cecece,cecece,d1d1d1,d1d1d1,d1d1d1,d1d1d1,d1d1d1,d1d1d1,d1d1d1,d1d1d1,d1d1d1,d1d1d1,d1d1d1,d1d1d1,d3d3d3,d3d3d3,d3d3d3,d3d3d3,d3d3d3,d3d3d3,d3d3d3,d3d3d3,d3d3d3,d3d3d3,d3d3d3,d6d6d6,d6d6d6,d6d6d6,d6d6d6,d6d6d6,d6d6d6,d6d6d6,d6d6d6,d8d8d8,d8d8d8,d8d8d8,d8d8d8,d8d8d8,d8d8d8,d8d8d8,d8d8d8,d8d8d8,d8d8d8,d8d8d8,d8d8d8,dbdbdb,dbdbdb,dbdbdb,dbdbdb,dbdbdb,dbdbdb,dbdbdb,dbdbdb,dbdbdb,dbdbdb,dbdbdb,dddddd,dddddd,dddddd,dddddd,dddddd,dddddd,dddddd,dddddd,e0e0e0,e0e0e0,e0e0e0,e0e0e0,e0e0e0,058784,058484,e0e0e0,028282,e2e2e2,e2e2e2,e2e2e2,e2e2e2,e2e2e2,e2e2e2,028282,02827c,e2e2e2,e2e2e2,e2e2e2,e2e2e2,e2e2e2,e2e2e2,e2e2e2,e5e5e5,027c7c,027c77,007c77,007a75,007770,00756d,e8e8e8,e8e8e8,e8e8e8,e8e8e8,00756b,007266,e8e8e8,e8e8e8,eaeaea,eaeaea,eaeaea,eaeaea,eaeaea,eaeaea,007260,00705b,006d59,eaeaea,eaeaea,eaeaea,ededed,ededed,ededed,ededed,006d51,006b4c,006b44,00683d,006633,f4f202,f4f207,f4f20f,f4f214,f4f21c,f4f421,f4f428,f4f42d,f4f433,f4f43a,f7f43f,f7f447,f7f44c,f7f454,f7f459,f7f460,f7f468,f7f46d,f7f475,f7f47a,f9f482,f9f487,f9f78c,f9f793,f9f79b,006626,00631c,f9f7ad,f9f9b2,f9f9ba,f9f9c1,f9f9c6,f9f9cc,f9f9d3,f9f9db,f9f9e0,f9f9e5,00600c,f9f9f2,006000',
// 'addToLegend':false},'FIA Forest Type',false,null,null,'250m spatial resolution forest type map derived from MODIS data using FIA data')

// usForestType = ee.ImageCollection([usForestType.set('system:time_start',ee.Date.fromYMD(2004,6,1))]);

// var mtbs = ee.ImageCollection('projects/USFS/LCMS-NFS/CONUS-Ancillary-Data/MTBS');
// mtbs = mtbs.map(function(img){return img.select([0],['burnSeverity']).byte().updateMask(img.neq(0).and(img.neq(6)))});
// var mtbsForViz = mtbs.map(addTimeBand);

// mtbsForViz = mtbsForViz.qualityMosaic('year');

// var cdlFirst = ee.Image('USDA/NASS/CDL/2016').select([0]);
// Map.addLayer(cdlFirst,{'addToLegend' : false},'Cropland Data Layer' ,false,null,null,'Latest CDL Layer for data visualization.  Right-click charting will plot all available CDL layers.')





// Map.addLayer(nlcd01.select([1]),{'min':0,'max':90,'palette':'000,555,FF0,F30,F00'},'NLCD 01 Impervious',false,'% impervious','FFF','NLCD 2001 percent impervious- depicts percent of a pixel that is impervious');

// Map.addLayer(nlcd06.select([1]),{'min':0,'max':90,'palette':'000,555,FF0,F30,F00'},'NLCD 06 Impervious',false,'% impervious','FFF','NLCD 2006 percent impervious- depicts percent of a pixel that is impervious');
// Map.addLayer(nlcd11.select([1]),{'min':0,'max':90,'palette':'000,555,FF0,F30,F00'},'NLCD 11 Impervious',false,'% impervious','FFF','NLCD 2011 percent impervious- depicts percent of a pixel that is impervious');

// Map.addLayer(nlcd11.select([2]),{'min':0,'max':90,'palette':'000,0F0'},'NLCD 11 Tree Canopy Cover',false,'% tree canopy cover','FFF','NLCD 2011 percent tree canopy cover- depicts percent of a pixel that has tree canopy cover');


// Map.addLayer(nlcd92.select([0]),{'addToLegend':false},'NLCD 92 Landcover',false);
// Map.addLayer(nlcd01.select([0]),{'addToLegend':false},'NLCD 01 Landcover',false);
// Map.addLayer(nlcd06.select([0]),{'addToLegend':false},'NLCD 06 Landcover',false);
// Map.addLayer(nlcd11.select([0]),{'addToLegend':false},'NLCD 11 Landcover',false);


// /////////////////////////////////////////////////////////////////

// //Bring in ids defoliation data
// var ids = ee.FeatureCollection('users/ianhousman/IDS/IDS_Defol');

// //Specify ids defoliation severities interested in
// var severities = [-1,0,1,2,3,4,11];

// //Denote dca_codes
// //From: https://www.fs.fed.us/foresthealth/technology/pdfs/Appendix_E_DCA_20141030.pdf
// var dca_codes = {
//   10000:'General Insects',
//   11000:'Bark Beetles',
//   12000:'Defoliators',
//   13000: 'Chewing Insects',
//   14000:'Sap Feeding Insects',
//   15000: 'Boring Insects',
//   16000: 'Seed/Cone/Flower/Fruit Insects',
//   17000: 'Gallmaker Insects',
//   18000: 'Insect Predators',
//   19000: 'General Diseases',
//   20000: 'Biotic Damage',
//   21000: 'Root/Butt Diseases',
//   22000: 'Stem Decays/Cankers',
//   23000: 'Parsitic/Epiphytic Plants',
//   24000: 'Decline Complexes/Dieback/Wilts',
//   25000: 'Foliage Diseases',
//   26000: 'Stem Rusts',
//   27000: 'Broom Rusts',
//   28000: 'Terminal, Shoot, and Twig Insects',
//   29000: 'Root Insects',
//   30000: 'Fire',
//   40000: 'Wild Animals',
//   50000: 'Abiotic Damage',
//   60000: 'Competition',
//   70000: 'Human Activities',
//   80000: 'Multi-Damage (Insect/Disease)',
//   85000: 'Plants',
//   90000: 'Other Damages and Symptoms',
//   99999: 'No Damage'
// };
// var damage_codes = {1:'Not Specified',
//   2:    'Mortality',
// 3   :'Crown Discoloration',
// 4   :'Crown Dieback',
// 5   :'Topkill',
// 6   :'Branch Breakage',
// 7   :'Main stem Broken or Uprooted',
// 8   :'Branch flagging',
// 9   :'No damage',
// 11: 'Mortality - Previously Undocumented',
// 12: 'Defoliation < 50% of leaves defoliated',
// 13: 'Defoliation 50-75% of leaves defoliated',
// 14: 'Defoliation > 75% of leaves defoliated',
// 18: 'Other Damage (known)',
// 19: 'Unknown Damage'
// };
// //Function for adding IDS data for a given year
// //Adds year, severity and dca code
// function addIDS(ids,year){
//   //Rasterize severity code
//   // var sevCol = severities.map(function(sev){
//   //   var idsSev = ids.filter(ee.Filter.equals('SEVERITY',sev));
//   //   var sevFill = ee.Image().paint(idsSev,1);
//   //   sevFill = sevFill.where(sevFill.mask(),sev).int16();
//   //   return sevFill.rename(['severity']);
//   // });
//   // sevCol = ee.ImageCollection.fromImages(sevCol).max();
  
  
//   var damCol = Object.keys(damage_codes).map(function(dam){
//     dam = parseInt(dam);
//     var idsDam = ids.filter(ee.Filter.equals('DAMAGE_TYP',dam));
//     var damFill = ee.Image().paint(idsDam,1);
//     damFill = damFill.where(damFill.mask(),dam).int16();
//     return damFill.rename(['damage_type']);
//   });
//   damCol = ee.ImageCollection.fromImages(damCol).max();
  
//   //Rasterize dca code
//   var dcaCol = Object.keys(dca_codes).map(function(dca_code){
//     dca_code = parseInt(dca_code);
//     var idsDCA = ids.filter(ee.Filter.and(ee.Filter.gte('DCA_CODE',dca_code),ee.Filter.lt('DCA_CODE',dca_code + 1000)));
    
//     var dcaFill = ee.Image().paint(idsDCA,1);
//     dcaFill = dcaFill.where(dcaFill.mask(),parseInt(dca_code/1000)).byte();
//     return dcaFill.rename(['dca']);
//   });
//   dcaCol = ee.ImageCollection.fromImages(dcaCol).mode();
  
//   //Rasterize year of ids data
//   var idsImg = ee.Image().paint(ids,1);
 
//   idsImg = idsImg.where(idsImg.mask(),year).rename(['year']);
 
//   // Map.addLayer(idsImg,{'min':startYear,'max':endYear,'palette':'FF0,F00','opacity':0.5},year.toString() + ' IDS',false);
//   return idsImg.addBands(damCol).addBands(dcaCol).int16().set('system:time_start',ee.Date.fromYMD(year,6,1).millis());

// }


// var idsStartYear = 1999;
// var idsEndYear = 2017;
// var idsYears = ee.List.sequence(idsStartYear,idsEndYear).getInfo();

// var idsCollection = idsYears.map(function(yr){
//   var idsYr = ids.filter(ee.Filter.equals('SURVEY_YEA',yr));
//   return addIDS(idsYr,yr);
// });
// idsCollection = ee.ImageCollection(idsCollection);
// var idsCollectionUnmasked = idsCollection;



// var mortStartYear = 1997;
// var mortEndYear = 2015;
// var mortYears = ee.List.sequence(mortStartYear,mortEndYear).getInfo();

// var mortCollection = mortYears.map(function(mortYear){
//   var mort = ee.FeatureCollection('users/ianhousman/IDS/IDS_Mort_' + mortYear.toString());
//   var img = addIDS(mort,mortYear);
//   return img
// })
// mortCollection = ee.ImageCollection(mortCollection);



// Map.addLayer(idsCollection.select([0]).max(),{'min':idsStartYear,'max':idsEndYear,'palette':'FF0,F00','opacity':0.8},'IDS Defoliation Survey Year-Most Recent',false,'','FFF','IDS Defoliation data- year of most recent survey.  Double-click chart to see exactly which years were surveyed');
// Map.addLayer(idsCollection.select([0]).count(),{'min':1,'max':Math.floor((idsEndYear-idsStartYear)/4),'palette':'FF0,F00','opacity':0.8},'IDS Defoliation Survey Count',false,'','FFF','IDS Defoliation data- number of years surveyed.  Double-click chart to see exactly which years were surveyed');
// // Map.addLayer(idsCollection.select([1]).max(),{'min':-1,'max':2,'palette':'F00,FF0','opacity':0.8},'IDS Defoliation Max Severity',false,'','FFF','IDS Defoliation data- severity of defoliation.  Double-click chart to see exactly which years were what severity');

// Map.addLayer(mortCollection.select([0]).max(),{'min':idsStartYear,'max':idsEndYear,'palette':'FF0,F00','opacity':0.8},'IDS Mortality Survey Year-Most Recent',false,'','FFF','IDS Mortality data- year of most recent survey.  Double-click chart to see exactly which years were surveyed');
// Map.addLayer(mortCollection.select([0]).count(),{'min':1,'max':Math.floor((idsEndYear-idsStartYear)/4),'palette':'FF0,F00','opacity':0.8},'IDS Mortality Survey Count',false,'','FFF','IDS Mortality data- number of years surveyed.  Double-click chart to see exactly which years were surveyed');
// // Map.addLayer(mortCollection.select([1]).max(),{'min':-1,'max':2,'palette':'F00,FF0','opacity':0.8},'IDS Mortality Max Severity',false,'','FFF','IDS Mortality data- severity of defoliation.  Double-click chart to see exactly which years were what severity');

// // Map.addLayer(idsCollection.select([2]).mode().randomVisualizer(),{},'IDS DCA Mode Collection',false);
// // Map.addLayer(idsCollectionUnmasked,{},'TS',false);

// /////////////////////////////////////////////////////////////////////////////
// // var startYear = 1985;
// // var endYear = 2016;
// // var changeThresh = 50;

// // var conusChange = ee.ImageCollection('projects/glri-phase3/science-team-outputs/conus-lcms-2018')
// //   .filter(ee.Filter.calendarRange(startYear,endYear,'year'));
// // conusChange = conusChange.map(function(img){
// //   var yr = ee.Date(img.get('system:time_start')).get('year')
// //   var change = img.gt(changeThresh);
// //   var conusChangeYr = ee.Image(yr).updateMask(change).rename(['change']).int16()
// //   return img.mask(ee.Image(1)).addBands(conusChangeYr);
// // })  
 
// // Map.addLayer(conusChange.select(['change']).min(),{'min':startYear,'max':endYear,'palette':'FF0,F00'},'CONUS LCMS',false);



// Map.addLayer(mtbsForViz.select([0]),{'min':1,'max':6,'palette':'006400,7fffd4,ffff00,ff0000,7fff00,ffffff','opacity':1,'addToLegend':false},'MTBS 1984-2015 Severity Composite',false,null,null,'Burn severity  from Monitoring Trends in Burn Severity project (mtbs.gov). Most recent year burned is displayed. Use right-click charting to chart all years (1984-2015)')
// Map.addLayer(mtbsForViz.select([1]),{'min':1984,'max':2015,'palette':'9400D3,4B0082,00F,0F0,FF0,FF7F00,F00'},'MTBS 1984-2015 Year of Most Recent Fire',false,'','FFF','Year of most recent fire from Monitoring Trends in Burn Severity project (mtbs.gov). Most recent year burned is displayed. Use right-click charting to chart all years (1984-2015)')


// var nlcdImperv = nlcd.select([1]);
// var nlcdTCC = ee.ImageCollection([nlcd11.select([2],['tcc']).set('system:time_start',ee.Date.fromYMD(2011,6,1))]);
// nlcd = nlcd.select([0]);

// // var plotsGLRIWGS84 = ee.FeatureCollection('ft:1r6KdUOrAnMHm-4bYIpWEDdpC9XUOFejcMaXxCAmw','geometry');
// // var plotsGLRINAD83 = ee.FeatureCollection('ft:1e6Cx5KWexvglX-B_IGWW8eKDPtAwdhfBe607kliS');
// // var plotsUSUPractice = ee.FeatureCollection('ft:1omDILKXadZhkvbaAidBITEcx_aFGiMMgk0lrWecc');
// // var plotsPilot = ee.FeatureCollection('ft:1TTpFIk3FraazYOw6ksQq__XrJunzE4EhoVskepDo','geometry');
// // var plotsPeugot = ee.FeatureCollection('ft:1aJPJNannk9KwZ1VmESqvbBMrDPbMyA_bV7d0-ckb','geometry');
// // var plotsGLRI_300_600 = ee.FeatureCollection('ft:1gT_9j13htipyVWRhNh3-DMf2cqQyYPa1mWLH2Gzo');
// // var plotsGLRI_600_927 = ee.FeatureCollection('ft:1EfFqz2tnib6l8NIXwC3NQp_5NIevZ2LNFr8PtgzK');
// // var otherPlots = ee.FeatureCollection('ft:1X-pNo3AXoZwg9pRsKRIEJ4RvJdyE8IUmHOq68-LQ');
// // var USUPlots2 = ee.FeatureCollection('ft:1X8gbR-JJY5BQ8M5NgYPiUOjSVdSa4an53d3w60SQ');
// // var rioPlots = ee.FeatureCollection('ft:1ul8Hzhan_NA8kUO7FeqCyhnfctkjJTIcZ-I6HBdQ');
// // var gaPrctPlots = ee.FeatureCollection('ft:1OPAbJIHYV6DdN1NwIwyasmFRP8Ag8i_Lrr9K1c57');
// // var plots = ee.FeatureCollection(gaPrctPlots.merge(otherPlots).merge(rioPlots).merge(USUPlots2).merge(plotsUSUPractice).merge(plotsGLRI_600_927));//.merge(otherPlots).merge(plotsGLRI_300_600).merge(plotsGLRIWGS84).merge(plotsPilot).merge(plotsPeugot));
// // plots = ee.Image();//.paint(plots,2,2)


// // var hardrock = ee.FeatureCollection('ft:1ikNCOJuM-z5zOtKED1eQeDtGlL8n_YzYGgPl7WOQ');
// // Map.addLayer(hardrock,{'palette':'FF0'},'Hardrock',false);


addDynamicToMap('https://fwsprimary.wim.usgs.gov/server/rest/services/Wetlands_Raster/ImageServer/exportImage?f=image&bbox=',
                'https://fwsprimary.wim.usgs.gov/server/rest/services/Wetlands/MapServer/export?dpi=96&transparent=true&format=png8&bbox=',
                8,11,
                'NWI',false,'National Wetlands Inventory data as viewed in https://www.fws.gov/wetlands/Data/Mapper.html from zoom levels >= 8')

// // Map.addLayer(plots,{'min':1,'max':1,'palette':'F00','addToLegend':false},'Plots')

// // var ft = ee.FeatureCollection('ft:1xqJ22rcmOqtNi0LfDus21n-7HFFEUYcCwqUHtVQm');
// // Map.addLayer(ft,{'palette':'F00'},'ft')


// /////////////////////////////////////////////////////////////////////////////////
// //Function to handle empty collections that will cause subsequent processes to fail
// //If the collection is empty, will fill it with an empty image
// function fillEmptyCollections(inCollection,dummyImage){                       
//   var dummyCollection = ee.ImageCollection([dummyImage.mask(ee.Image(0))]);
//   var imageCount = inCollection.toList(1).length();
//   return ee.ImageCollection(ee.Algorithms.If(imageCount.gt(0),inCollection,dummyCollection));

// }
// //////////////////////////////////////////////////////////////////////////


// var years = ee.List.sequence(1984,2017);

// var mtbsPalette = '006400,7fffd4,ffff00,ff0000,7fff00,ffffff';
// // var mtbs = ee.ImageCollection('users/ianhousman/MTBS/Collection');

// mtbs = mtbs.map(function(img){return img.unmask()})
// var dummyCDL = ee.Image(cdl.first()).int16();
// var dummyNLCD = ee.Image(nlcd.first()).int16();
// var dummyMTBS = ee.Image(mtbs.first()).int16();
// var dummyNLCDImperv = ee.Image(nlcdImperv.first()).int16();
// var dummyGlobCover = ee.Image(globCover.first()).int16();
// var dummyTCC = ee.Image(nlcdTCC.first()).int16();
// var dummyUSFT = ee.Image(usForestType.first()).int16();
// var dummyMTBS = ee.Image(mtbs.first()).int16();
// var dummyIDS = ee.Image(idsCollectionUnmasked.first()).int16();
// var dummyMort = ee.Image(mortCollection.first()).int16();
// // var dummyLCMS = ee.Image(conusChange.first()).int16();
// var c =years.map(function(yr){
//   var startDate = ee.Date.fromYMD(yr,1,1);
//   var endDate = ee.Date.fromYMD(yr,12,31);
  
//   var cdlT = cdl.filterDate(startDate,endDate);
//   cdlT = ee.Image(fillEmptyCollections(cdlT,dummyCDL).first()).int16();
  
//   var nlcdImpervT = nlcdImperv.filterDate(startDate,endDate);
//    nlcdImpervT = ee.Image(fillEmptyCollections(nlcdImpervT,dummyNLCDImperv).first()).int16();

//   var nlcdT = nlcd.filterDate(startDate,endDate);
//   nlcdT = ee.Image(fillEmptyCollections(nlcdT,dummyNLCD).first()).int16();
  
//    var mtbsT = mtbs.filterDate(startDate,endDate)
//     mtbsT = ee.Image(fillEmptyCollections(mtbsT,dummyMTBS).first()).int16();

//    var globCoverT = globCover.filterDate(startDate,endDate);
//    globCoverT = ee.Image(fillEmptyCollections(globCoverT,dummyGlobCover).first()).int16();

//    var nlcdTCCT = nlcdTCC .filterDate(startDate,endDate);
//    nlcdTCCT = ee.Image(fillEmptyCollections(nlcdTCCT,dummyTCC).first()).int16();
//   var usFT = usForestType.filterDate(startDate,endDate);
//   usFT = ee.Image(fillEmptyCollections(usFT,dummyUSFT).first()).int16();

//   var idsT = idsCollectionUnmasked.filterDate(startDate,endDate);
//   idsT = ee.Image(fillEmptyCollections(idsT,dummyIDS).first()).int16().select([1,2],['defol_damage','defol_dca']);
 
//   var mortT = mortCollection.filterDate(startDate,endDate);
//   mortT = ee.Image(fillEmptyCollections(mortT,dummyMort).first()).int16().select([1,2],['mort_damage','mort_dca']);
 
//   // var lcmsT = conusChange.filterDate(startDate,endDate);
//   // lcmsT = ee.Image(fillEmptyCollections(lcmsT,dummyLCMS).first()).int16().select([0],['LCMS Change Prob']);
 
//   var out = ee.Image(1).addBands(mtbsT).addBands(cdlT).addBands(nlcdT).addBands(nlcdImpervT).addBands(nlcdTCCT).addBands(globCoverT).addBands(usFT).addBands(idsT).addBands(mortT).int16();//.rename(['MTBS','CDL','NLCD LC','NLCD Imperv','NLCD TCC','GlobCover','FIA Forest Types'])
//   out = out.select([1,2,3,4,5,6,7,8,9,10,11])
//   return out.set('system:time_start',ee.Date.fromYMD(yr,6,1).millis())
  
// })
// c = ee.ImageCollection.fromImages(c);

// chartCollection =c;
///////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////
// function despikeCollection(c){
  
//   c = c.toList(10000,0);
//   var left = c.slice(0,-2);
//   var center = c.slice(1,-1);
//   var right = c.slice(2,null);
  
//   var ids = ee.List.sequence(0,left.length().subtract(1));
//   var despiked =ids.map(function(id){
//     var l =ee.Image(left.get(id));
//     var c = ee.Image(center.get(id));
//     var r = ee.Image(right.get(id));
//     var mean = (l.add(r)).divide(2);

//     // var isFillable = (l.mask().eq(1).or(r.mask().eq(1))).and(c.mask().eq(0));
//     // c = c.unmask(isFillable)
//     // c = c.where(isFillable,mean)
   
//     var d1 = c.subtract(l);
//     var d2 = r.subtract(c);
    
//     var negativeSpike = (d1.lt(-spikeThresh).and(d2.gt(spikeThresh/2)));
//     var positiveSpike = (d1.gt(spikeThresh).and(d2.lt(-spikeThresh)))
//     // .or((d1.gt(spikeThresh).and(d2.lt(-spikeThresh))));
//     // .or((d1.lt(-spikeThresh).and(d2.gt(spikeThresh/2))));
    
    
//     var out = c.where(negativeSpike.or(positiveSpike),mean)
//     .copyProperties(c,['system:time_start']);
//     // center = center.insert(id,out)
//     return out
//   })
//   despiked = ee.List([c.get(0),despiked,c.get(-1)]).flatten()
  
//   despiked = ee.ImageCollection.fromImages(despiked)
  
  
//   return despiked
// }

// var spikeThresh = 0.1;//Threshold for how much time series can change in successive years (i.e. 0.1 means timeseries cannot rise and fall or fall and rise more than that)

// var targetDay = 150;
// var sa = eeBoundsPoly;
// //////////////////////////////////////////////////////
// var sensorBandDictLandsatTOA =ee.Dictionary({'L8' : ee.List([1,2,3,4,5,9,6,'cloud']),
//                         'L7' : ee.List([0,1,2,3,4,5,7,'cloud']),
//                         'L5' : ee.List([0,1,2,3,4,5,6,'cloud']),
//                         'L4' : ee.List([0,1,2,3,4,5,6,'cloud'])
//   });
// var bandNamesLandsatTOA = ee.List(['blue','green','red','nir','swir1','temp','swir2','cloud']);
// //////////////////////////////////////////////////////////////////////////////////
// //Various preset min/max stretches for data visualization
// var vizParamsCO1 = {'min': 0.05,'max': [0.3,0.6,0.35],   'bands':'swir1,nir,red'};
// var vizParamsCO2 = {'min': 0.15,'max': [0.35,0.8,0.4],   'bands':'swir1,nir,red', 'gamma': 1.6};
// var vizParamsCO3 = {'min': 0.05,'max': [0.3,0.4,0.4],   'bands':'swir1,nir,red', 'gamma':1.6};
// var vizParamsFalse = {'min': 0.1,'max': [0.3,0.3,0.3],   'bands':'nir,swir1,red'};
// var vizParamsViz = {'min': 0.05, 'max': 0.3,'bands': 'red,green,blue', 'gamma': 1.6};
// // var vizParams = vizParamsCO2;
// var vizParams = {'min': 0.15,'max': [0.35,0.8,0.4],   'bands':'swir1,nir,red', 'gamma': 1.6};//eval("("+compositeViz+")");//{'min':0.1,'max':0.4,'bands':'swir1,nir,red'}
// var trendVizParams = {'min':-0.03,'max':0.01,'palette':'F00,888,0F0'};//eval('('+trendViz+')');
// var diffVizParams = {'min':-0.3,'max':0.2,'palette':'F00,888,0F0'};//eval('('+diffViz+')');
// //////////////////////////////////////////////////////


// function bustClouds(img){
//   var cs = ee.Algorithms.Landsat.simpleCloudScore(img).select(['cloud'])
//   var m = cs.lt(20);
//   return img.updateMask(m)
  
// }
// function getCloudMask(img){
//   var cs = ee.Algorithms.Landsat.simpleCloudScore(img).select(['cloud'])
//   var m = cs.lt(20);
//   return img.addBands(m.rename(['cloud']))
// }

// function addIndices(img){
//   var ndvi = img.normalizedDifference(['nir','red']).rename('ndvi');
//   var nbr = img.normalizedDifference(['nir','swir2']).rename('nbr');
//   return img.addBands(ndvi).addBands(nbr);

// }
// /////////////////////////////////////////////////////////////////
// //Adds the float year with julian proportion to image
// function addDateBand(img){

//   var d = ee.Date(img.get('system:time_start'))
//   var y = d.get('year')
//   var d = y.add(d.getFraction('year'))
//   // var pt = ee.PixelType('float',1950,2100)
//   var db = ee.Image(d).rename('year').float()
//   return img.addBands(db)
//   // .copyProperties(img)
// }
// ///////////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////
// //Function for finding dark outliers in time series
// //Masks pixels that are dark, and dark outliers
// function simpleTDOM2(c){
//   var shadowSumBands = ['nir','swir1'];
//   var irSumThresh = 0.4;
//   var zShadowThresh = -1.2
//   //Get some pixel-wise stats for the time series
//   var irStdDev = c.select(shadowSumBands).reduce(ee.Reducer.stdDev());
//   var irMean = c.select(shadowSumBands).mean();
//   var bandNames = ee.Image(c.first()).bandNames();
//   print('bandNames',bandNames);
//   //Mask out dark dark outliers
//   c = c.map(function(img){
//     var z = img.select(shadowSumBands).subtract(irMean).divide(irStdDev);
//     var irSum = img.select(shadowSumBands).reduce(ee.Reducer.sum())
//     var m = z.lt(zShadowThresh).reduce(ee.Reducer.sum()).eq(2).and(irSum.lt(irSumThresh)).not();
    
//     return img.updateMask(img.mask().and(m));
//   })
  
//   return c.select(bandNames)
// }
// /////////////////////////////////////////////////////////////////

// function getLandsat(startDate,endDate,startJulian,endJulian){
//   var l5s = ee.ImageCollection('LT5_L1T_TOA')
//               .filterBounds(sa)
//               .filterDate(startDate,endDate)
//               .filter(ee.Filter.calendarRange(startJulian,endJulian))
//               .map(getCloudMask)
//               .select(sensorBandDictLandsatTOA.get('L5'),bandNamesLandsatTOA);
//   var l7s = ee.ImageCollection('LE7_L1T_TOA')
//               .filterBounds(sa)
//               .filterDate(startDate,endDate)
//               .filter(ee.Filter.calendarRange(startJulian,endJulian))
//               .map(getCloudMask)
//               .select(sensorBandDictLandsatTOA.get('L7'),bandNamesLandsatTOA);
//   var l8s = ee.ImageCollection('LC8_L1T_TOA')
//               .filterBounds(sa)
//               .filterDate(startDate,endDate)
//               .filter(ee.Filter.calendarRange(startJulian,endJulian))
//               .map(getCloudMask)
//               .select(sensorBandDictLandsatTOA.get('L8'),bandNamesLandsatTOA);
//   return ee.ImageCollection(l5s.merge(l7s).merge(l8s));
// }

// function getLandsat2(startDate,endDate,startJulian,endJulian){
//   var l5s = ee.ImageCollection('LT5_L1T_TOA')
//               .filterBounds(sa)
//               .filterDate(startDate,endDate)
//               .filter(ee.Filter.calendarRange(startJulian,endJulian))
//               .map(bustClouds)
//               .select(sensorBandDictLandsatTOA.get('L5'),bandNamesLandsatTOA);
//   var l7s = ee.ImageCollection('LE7_L1T_TOA')
//               .filterBounds(sa)
//               .filterDate(startDate,endDate)
//               .filter(ee.Filter.calendarRange(startJulian,endJulian))
//               .map(bustClouds)
//               .select(sensorBandDictLandsatTOA.get('L7'),bandNamesLandsatTOA);
//   var l8s = ee.ImageCollection('LC8_L1T_TOA')
//               .filterBounds(sa)
//               .filterDate(startDate,endDate)
//               .filter(ee.Filter.calendarRange(startJulian,endJulian))
//               .map(bustClouds)
//               .select(sensorBandDictLandsatTOA.get('L8'),bandNamesLandsatTOA);
//   return ee.ImageCollection(l5s.merge(l7s).merge(l8s));
// }



// // Function to apply and correct for radiometry, selecting greenest pixel, closest to target day
// var addqa = function(image) {
//           image = image.rename(['blue','green','B4', 'B5', 'B6', 'B10', 'B7','cloud','ndvi','nbr'])
//           var temp = image.select('B10').focal_min();
//           var cweight = image.metadata('CLOUD_COVER').subtract(100).multiply(-1);
//           var targday = ee.Number(targetDay);
          
          
//           var timestamp = ee.Number(image.get('system:time_start'));
//           var day = ee.Number(ee.Date(timestamp).getRelative('day', 'year'));
//           var time = ee.Number(ee.Date(timestamp).get('minute'));
          
//           var daydif = targday.subtract(day).abs();
//           var daydif2 = ee.Number(365).subtract(daydif);
          
//           var pi = ee.Number(3.14159265359);
//           var pipi = pi.multiply(2);
//           var deg2rad = pi.divide(180);
          
//           var part2 = day
//             .add(284)
//             .divide(36.25);
          
//           var part22 = pipi
//             .multiply(part2);
//           var part222 = part22.sin();
          
//           var dec_angle = deg2rad
//             .multiply(23.45)
//             .multiply(part222);
            
//           //Hour angle
//           //var local_time = 
//           var hour_angle = deg2rad
//             .multiply(-22.5);
          
//           //Per-pixel latitude and longitude
//           //var latlon = ee.Image.pixelLonLat();
//           var lat = ee.Image.pixelLonLat().select('latitude')
//             .multiply(deg2rad);
          
//           //Solar elevation angle
//           var cosh = hour_angle.cos();
//           var cosd = dec_angle.cos();
//           var coslat = lat.cos();
//           var sind = dec_angle.sin();
//           var sinlat = lat.sin();
          
//           var solar_elev1 = coslat
//             .multiply(cosd)
//             .multiply(cosh);
          
//           var solar_elev2 = sinlat
//             .multiply(sind);
            
//           var solar_elev = solar_elev1
//             .add(solar_elev2);
            
//           var toa_cor2 = solar_elev.sin();
          
//           var adj4 = image.select('B4').float().divide(toa_cor2);//.float().add(corb3);
//           var adj5 = image.select('B5').float().divide(toa_cor2);//.add(corb4);
//           var adj6 = image.select('B6').float().divide(toa_cor2);//.float().add(corb5);
//           var adj7 = image.select('B7').float().divide(toa_cor2);//.float().add(corb7);
          
//           var adj4_1 = adj4.multiply(10000);
//           var adj5_1 = adj5.multiply(10000);
//           var adj6_1 = adj6.multiply(10000);
//           var adj7_1 = adj7.multiply(10000);
          
//           //var ndvi = (adj5.subtract(adj4)).divide(adj5.add(adj4));
//           var ndvi = (image.select('B5').subtract(image.select('B4'))).divide(image.select('B5').add(image.select('B4')));
//           var weight = ndvi.multiply(temp);
//           //var weight = ndvi.multiply(1);
//           var cweight2 = weight.multiply(daydif2);//.multiply(cweight);
          
//           var output = image.addBands(adj4_1).addBands(adj5_1).addBands(adj6_1).addBands(temp).addBands(adj7_1).addBands(weight).addBands(cweight2);
          
//           return output.rename(ee.List(['blue','green','red','nir','swir1','temp','swir2','cloud','ndvi','nbr','B4_1','B5_1','B6_1','B10_1','B7_1','B5_2','B5_3']));
//           // end of the function
// };

// var ls =  getLandsat(startDate1,endDate2,startJulian,endJulian).map(addIndices);
// print(ls.first())
// var ls1 = ls.filterDate(startDate1,endDate1);//,startJulian,endJulian);
// var ls2 = ls.filterDate(startDate2,endDate2);//,startJulian,endJulian);


// // Map.addLayer(ls1.median(),vizParams,'Raw First Composite');
// // Map.addLayer(ls2.median(),vizParams,'Raw Second Composite');

// ls1 = ls1.map(addqa);
// ls2 = ls2.map(addqa);

// var lsAlt = simpleTDOM2(ls.map(function(img){return img.updateMask(img.mask().and(img.select(['cloud'])))}))
// var l1Alt = lsAlt.filterDate(startDate1,endDate1).median();
// var l2Alt = lsAlt.filterDate(startDate2,endDate2).median();

// // var y1 = ee.Date(startDate1).get('year');
// // var y2 = ee.Date(endDate2).get('year');
// // var compositingPeriod = 1;
// // var years = ee.List.sequence(y1,y2,compositingPeriod);
// // var dummyImage = ee.Image(lsAlt.first());

// // lsAlt = years.map(function(yr){
// //     var startDateT = ee.Date.fromYMD(yr,1,1);
// //     var endDateT = ee.Date.fromYMD(yr+compositingPeriod-1,12,31);
    
// //     var lsAltT = lsAlt.filterDate(startDateT,endDateT);
// //     // lsAltT = fillEmptyCollections(lsAltT,dummyImage);
// //     return lsAltT.median().set('system:time_start',ee.Date.fromYMD(yr,6,1).millis());
// // });
// lsAlt = lsAlt.map(addDateBand);
// // lsAlt =  despikeCollection(lsAlt);
// var linearModel = lsAlt.select(['year',whichIndex]).reduce(ee.Reducer.linearFit());
// chartCollection  = lsAlt.map(function(img){
//               var predicted = img.select(['year']).multiply(linearModel.select([0])).add(linearModel.select([1])).rename('Predicted '+whichIndex)
//               return img.select(['ndvi','nbr'  ]).addBands(predicted);
//             }).sort('system:time_start')
// var trend = linearModel.select([0]);
// var l1 =   ls1.qualityMosaic('B5_3');
// var l2 =   ls2.qualityMosaic('B5_3');
// Map.addExport(l1.float(),'FREL Composite 1',30,30,450,15,false,vizParams);
// Map.addExport(l2.float(),'FREL Composite 2',30,30,450,15,false,vizParams);
// Map.addExport(l1Alt.float(),'Alt Composite 1',30,30,450,15,false,vizParams);
// Map.addExport(l2Alt.float(),'Alt Composite 2',30,30,450,15,false,vizParams);


// Map.addREST(function(coord, zoom) {
//                         var quadCode = tileXYZToQuadKey(coord.x, coord.y, zoom)
//                         // var url = 'http://t0.tiles.virtualearth.net/tiles/a'+quadCode +'.jpeg?g=854&mkt=en-US&token=AuuPPqhBI42-Eo2fJ3uRFFrvhzM93YmEIZFbOtqjFPNW2V641ZvU4g53ivnlVT8b'
//                         var url = 'http://ecn.t1.tiles.virtualearth.net/tiles/h'+quadCode +'.jpeg?g=5256&mkt=en-US&token=AuuPPqhBI42-Eo2fJ3uRFFrvhzM93YmEIZFbOtqjFPNW2V641ZvU4g53ivnlVT8b'
                    
//                     return url
//                 },'Bing Hybrid',false,19,'Bing high resolution satellite data with vector data')

// Map.addREST(standardTileURLFunction('http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/'),'ESRI World Imagery',false,19, 'ESRI worldwide high resolution imagery')
// // Map.addREST(standardTileURLFunction('https://nowcoast.noaa.gov/arcgis/rest/services/nowcoast/forecast_meteoceanhydro_sfc_ndfd_windgust_offsets/MapServer/tile/'),'Temp',false,19, 'ESRI worldwide high resolution imagery')

// Map.addREST(function(coord, zoom) {
//                     var tilesPerGlobe = 1 << zoom;
//                     var x = coord.x % tilesPerGlobe;
//                     if (x < 0) {x = tilesPerGlobe+x;}
//                     return "http://api.tiles.mapbox.com/v4/digitalglobe.nal0mpda/" + zoom + "/" + x + "/" + coord.y + ".png?access_token=pk.eyJ1IjoiZGlnaXRhbGdsb2JlIiwiYSI6ImNpczBqYml6bzAzb2Yyb3A2YWI1aWVvc2sifQ.ZWTVvUKerHzA6AX41LS9LQ";
//                 },'Digital Globe Hybrid',false,19,'Digital Globe recent high resolution satellite data with vector data')

// // Map.addExport(trend.float(),'Trend '+whichIndex,30,30,450,15,false,trendViz );
// // Map.addExport(trend.float(),'Trend',30,30,450,15,false,trendViz);
// Map.addLayer(l1,vizParams  ,'FREL Composite 1',false)
// Map.addLayer(l2,vizParams,'FREL Composite 2',false)

// Map.addLayer(l1Alt,vizParams,'Alt Composite 1',false)
// Map.addLayer(l2Alt,vizParams,'Alt Composite 2',false)
// // print(ls1.first());

// var dFREL = l2.subtract(l1).select([whichIndex]);
// var dAlt = l2Alt.subtract(l1Alt).select([whichIndex]);

// Map.addExport(trend.float(),'Trend '+whichIndex,30,30,450,15,false,trendVizParams );
// Map.addExport(dFREL.float(),'FREL Difference '+whichIndex,30,30,450,15,false,diffVizParams );
// Map.addExport(dAlt.float(),'Alt Difference '+whichIndex,30,30,450,15,false,diffVizParams );

// Map.addLayer(dFREL,diffVizParams ,'FREL d'+whichIndex,false,'d' + whichIndex,'FFF','Difference in '+whichIndex+' between specified compositing periods using FAO FREL compositing method');
// Map.addLayer(dAlt,diffVizParams,'Alt d'+whichIndex,false,'d' + whichIndex,'FFF','Difference in '+whichIndex+' between specified compositing periods using Google/RSAC compositing method');

// Map.addLayer(trend,trendVizParams,'Trend '+whichIndex,false, whichIndex + '/yr','FFF','Linear trend in '+whichIndex+' across entire time period using Google/RSAC compositing method-prepared data');



// //Functions to perform basic clump and elim
// function sieve(image,mmu){
//   var connected = image.connectedPixelCount(mmu+2)
//   var elim = connected.gt(mmu);
//   var mode = image.focal_mode(mmu/2,'circle');
//   mode = mode.mask(image.mask());
//   var filled = image.where(elim.not(),mode)
//   return filled;
// }

// var thematic = ee.Image(0);
// // thematic = thematic.where(trend.gte(-0.01),2);
// thematic = thematic.where(trend.gte(0.002),2);
// thematic = thematic.where(trend.lt(-0.008),3);
// thematic = thematic.where(trend.lt(-0.015),4);

// // thematic = thematic.focal_mode(4,'circle');
// // Map.addLayer(thematic.updateMask(thematic.gte(4)),{'min':1,'max':4,'palette':'0F0,888,FF0,F00'},'Thematic Decline Only')
// Map.addLayer(thematic.updateMask(thematic.gt(1)),{'min':2,'max':4,'palette':'0F0,FF0,F00'},'Thematic Change',false)







// // var linearModel = ls.select(['year','nbr']).reduce(ee.Reducer.linearFit());

// // chartCollection  = ls.map(function(img){
// //               var predicted = img.select(['year']).multiply(linearModel.select([0])).add(linearModel.select([1])).rename('Predicted NBR')
// //               return img.select(['ndvi','nbr']).addBands(predicted);
// //             }).sort('system:time_start')
// // ls1 = ls.filterDate(startDate1,endDate1).median();
// // ls2 = ls.filterDate(startDate2,endDate2).median();

// // var water1 = waterScore(ls1);
// // var water2 = waterScore(ls2);
// // var dWater = water2.subtract(water1);
// // Map.addLayer(water1.mask(water1.gt(0.5)),{'min':0.2,'max':0.8,'palette':'888,00F'},'Water Pre',false);
// // Map.addLayer(water2.mask(water2.gt(0.5)),{'min':0.2,'max':0.8,'palette':'888,00F'},'Water Post',false);
// // Map.addLayer(dWater.mask(dWater.lt(-0.6).or(dWater.gt(0.6))),{'min':-1.5,'max':1.5,'palette':'F00,888,00F'},'Water Change',false);

// // var var1 = ls1.select(['nbr']).reduce(ee.Reducer.stdDev());
// // Map.addLayer(var1,{'min':0,'max':0.0001},'var1')
// // var forest1 = l5s1.select(['nbr','ndvi']).reduce(ee.Reducer.min()).gt(0.4);
// // Map.addLayer(forest1.mask(forest1),{'min':1,'max':1,'palette':'008000'},'Forest Mask 1')


// // d = ls2.select(['nbr']).subtract(ls1.select(['nbr']));
// // trend = linearModel.select([0]);
// // var vizParams = eval('('+compositeViz+')');//{'min':0.1,'max':0.4,'bands':'swir1,nir,red'}
// // var trendVizParams = eval('('+trendViz+')');
// console.log('v')
// console.log(vizParams)
// var exportVizParams = {'min':-0.2,'max':0}

// // Map.addLayer(ls1,vizParams, 'Pre Composite',false,null,null,'Composite of Landsat data for speficied pre date range')
// // Map.addLayer(ls2,vizParams, 'Post Composite',false,null,null,'Composite of Landsat data for speficied post date range')



// // Map.addLayer(d,{'min':-0.2,'max':0.1,'palette':'F00,888,0F0'}, 'Difference ' +whichIndex,false,'d','FFF','Difference in NBR between specified compositing periods')
// // Map.addLayer(trend,trendVizParams,'Trend',false,'nbr/year','FFF','Trend of NBR for specified date range');


// // var image1 = ee.Image(1);
// // Map.addExport(ls1.float(),'Pre Composite',30,30,450,15,true,vizParams);
// // Map.addExport(ls2.float(),'Post Composite',30,30,450,15,true,vizParams);
// // Map.addExport(d,'Difference Image',30,30,450,15,false,{'min':-0.2,'max':0.1,'palette':'F00,888,0F0'});
// // Map.addExport(trend,'Trend Image',30,30,450,15,false,trendVizParams);
// // Map.addExport(image1,'Testt',30,20,50,10,false)



// var eProvChange = ee.Image('users/ianhousman/Zambia/Zambia_Eastern_Province_Change_1990_2010_Unclipped');
// eProvChange = eProvChange.mask(eProvChange.gt(0));

// var RCMRD2010 = ee.Image('users/ianhousman/Zambia/zambia_landcover_2010_scheme_ii');
// var RCMRD2000 = ee.Image('users/ianhousman/Zambia/Zambia_Landcover2_2000_Scheme_II');

// RCMRD2010 = RCMRD2010.mask(RCMRD2010.mask().and(RCMRD2010.neq(255)));
// RCMRD2000 = RCMRD2000.mask(RCMRD2000.mask().and(RCMRD2000.neq(255)));

// var hansen = ee.Image('UMD/hansen/global_forest_change_2015');
// Map.addLayer(hansen.select(['treecover2000']).mask(hansen.select(['treecover2000']).gt(0)),{'min':1,'max':90,'palette':'500000,00FF00'},'Hansen Tree Cover 2000',false,'% Forest Cover','FFF','Tree cover for year 2000 from Hansen Global Forest Change suite of products');


// Map.addLayer(RCMRD2000,{'min': 0, 'max': 12, 'addToLegend':false, 'palette': '89FFFF,005800,009300,78FF58,808000,FFFF00,805000,FF5000,FF7878,009393,0000FF,575757,FFFFFF'},'RCMRD 2000',false,'Cover class','FFF','Land-use/land-cover classification from RCMRD for year 2000');
// Map.addLayer(RCMRD2010,{'min': 0, 'max': 12, 'addToLegend':false, 'palette': '89FFFF,005800,009300,78FF58,808000,FFFF00,805000,FF5000,FF7878,009393,0000FF,575757,FFFFFF'},'RCMRD 2010',false,'Cover class','FFF','Land-use/land-cover classification from RCMRD for year 2010');

// Map.addLayer(eProvChange,{'min': 1, 'max': 8, 'addToLegend':false,'palette': 'FF8060,006029,003100,FF0000,600029,FFFF00,FF5000,606000'},'USFS Eastern Prov Change',false,'Cover/change class','FFF','Cover and change classes from USFS RSAC forest trend product');
// // Map.addLayer(hansen.select(['loss']).mask(hansen.select(['loss']).gt(0)),{'min':0,'max':1},'loss')
// // Map.addLayer(hansen.select(['gain']).mask(hansen.select(['gain']).gt(0)),{'min':0,'max':1},'gain')
// Map.addLayer(hansen.select(['lossyear']).mask(hansen.select(['lossyear']).gt(0)).add(2000),{'min':2000,'max':2015,'palette':'FFFF00,FF0000'},'Hansen Loss Year',false,'forest loss','FFF');

};
