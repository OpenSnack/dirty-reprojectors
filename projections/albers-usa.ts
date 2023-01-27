import { geoAlbersUsa } from 'd3-geo';
const R = 6378137.0 // radius of Earth in meters
export const albersUsa = geoAlbersUsa().translate([0, 0]).scale(R);
