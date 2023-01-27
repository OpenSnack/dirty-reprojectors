import { GeoProjection } from 'd3-geo';
import { Geometry } from 'geojson';
import projections from '../projections';
type ReprojectOptions = {
    forward?: string | GeoProjection;
    reverse?: string | GeoProjection;
    projections: {
        [key: string]: GeoProjection;
    };
};
/**
 * Reprojects the given geometry coordinate array _in place_, with
 * unprojectable points or degenerate geometries removed. If both
 * `options.forward` and `options.reverse` are supplied, then `forward` is
 * performed first.
 *
 * @param {Object} options
 * @param {Function|string} [options.forward] The forward projection to use.
 * @param {Function|string} [options.reverse] The reverse projection to use.
 * @param {Object} [options.projections] A map of named projections to use.  If provided, then string values of `options.forward` or `options.reverse` will be used as keys to look up the projection function in `options.projections`.  For an extensive list provided by d3-geo-projection, use `require('dirty-reprojectors/projections')`.
 * @param {Object} geometry A GeoJSON geometry object
 */
declare function reproject(options: ReprojectOptions, geometry: Geometry): Geometry;
export { reproject, projections };
