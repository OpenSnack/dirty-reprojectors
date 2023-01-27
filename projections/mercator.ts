import { geoMercator } from 'd3-geo';
const R = 6378137.0; // radius of Earth in meters
module.exports = geoMercator().translate([0, 0]).scale(R);
