# EU_soilproperties
Mapping EU soil properties using point data [LUCAS topsoil](http://esdac.jrc.ec.europa.eu/content/lucas-2009-topsoil-data) and [GEMAS](http://gemas.geolba.ac.at/) data sets.

Targeted soil properties (0-30 cm depths):
  - soil organic carbon,
  - soil pH in H2O and CaCl2,
  - soil texture fractions sand silt and clay,
  - selection of heavy metals / micro-nutrients As, Cd, Cu, Pb, Zn,

![Distribution of GEMAS points over Europe](https://github.com/ISRICWorldSoil/EU_soilproperties/blob/master/points/GEMASlogZn.png)

Conversion of predictions and units:

  - `log_As_topsoil_1km.tif`: log of As concentration; to convert to ppm use `exp(x/100)-1`;
  - `log_Cd_topsoil_1km.tif`: log of Cd concentration; to convert to ppm use `exp(x/100)-1`;
  - `log_Cu_topsoil_1km.tif`: log of Cu concentration; to convert to ppm use `exp(x/100)-1`;
  - `log_Zn_topsoil_1km.tif`: log of Zn concentration; to convert to ppm use `exp(x/100)-1`;
  - `ORCDRC_topsoil_1km.tif`: soil organic carbon content in g/kg;
  - `PHICAL_topsoil_1km.tif`: soil pH in H2O; to convert to pH units divide by 10;
  - `PHIHOX_topsoil_1km.tif`: soil pH in H2O; to convert to pH units divide by 10;
  - `SLTPPT_topsoil_1km.tif`: silt content in w%;
  - `SNDPPT_topsoil_1km.tif`: sand content in w%;
  - `CLYPPT_topsoil_1km.tif`: clay content in w%;

Contacts: Tom.Hengl@isric.org and Gerard.Heuvelink@wur.nl  
