import { geoAlbers, geoConicEqualArea, GeoStream } from 'd3-geo';
import { multiplex, isNil } from './helpers';

function geoAlbersUsaTerritories() {
    const epsilon = 0.000001;
    let cache: GeoStream | null,
        cacheStream: GeoStream | null,
        lower48Point: GeoStream,
        alaskaPoint: GeoStream,
        hawaiiPoint: GeoStream,
        puertoRicoPoint: GeoStream,
        guamMarianaPoint: GeoStream,
        americanSamoaPoint: GeoStream,
        point;
        
    const lower48 = geoAlbers(),
        alaska = geoConicEqualArea().rotate([154, 0]).center([-2, 58.5]).parallels([55, 65]),
        hawaii = geoConicEqualArea().rotate([157, 0]).center([-3, 19.9]).parallels([8, 18]),
        puertoRico = geoConicEqualArea().rotate([66, 0]).center([0, 18]).parallels([8, 18]),
        guamMariana = geoConicEqualArea().rotate([-145, 0]).center([0, 16]).parallels([10, 20]),
        americanSamoa = geoConicEqualArea().rotate([170, 0]).center([0, -14]).parallels([-14, 0]),
        pointStream = { point: function(x: number, y: number) { point = [x, y]; } } as GeoStream;

    function albersUsaTerritories(coordinates: [number, number]): [number, number] {
        const x = coordinates[0];
        const y = coordinates[1];
        point = null;
        return (lower48Point.point(x, y), point)
            || (alaskaPoint.point(x, y), point)
            || (hawaiiPoint.point(x, y), point)
            || (puertoRicoPoint.point(x, y), point)
            || (guamMarianaPoint.point(x, y), point)
            || (americanSamoaPoint.point(x, y), point);
    }

    function reset() {
        cache = cacheStream = null;
        return albersUsaTerritories;
    }

    albersUsaTerritories.invert = function(coordinates: [number, number]) {
        const k = lower48.scale(),
            t = lower48.translate(),
            x = (coordinates[0] - t[0]) / k,
            y = (coordinates[1] - t[1]) / k;
        return (y >= 0.120 && y < 0.234 && x >= -0.390 && x < -0.185 ? alaska
            : y >= 0.166 && y < 0.234 && x >= -0.185 && x < -0.080 ? hawaii
                : y >= 0.204 && y < 0.234 && x >= 0.300 && x < 0.380 ? puertoRico
                    : y >= 0.050 && y < 0.210 && x >= -0.450 && x < - 0.390 ? guamMariana
                        : y >= 0.210 && y < 0.234 && x >= -0.450 && x < -0.390 ? americanSamoa
                            : lower48).invert?.(coordinates) ?? null;
    };

    albersUsaTerritories.stream = function(stream: GeoStream) {
        return cache && cacheStream === stream
            ? cache
            : cache = multiplex([
                lower48.stream(cacheStream = stream),
                alaska.stream(stream),
                hawaii.stream(stream),
                puertoRico.stream(stream),
                guamMariana.stream(stream),
                americanSamoa.stream(stream)
            ]);
    };

    function precision(precision: number): ReturnType<typeof reset>;
    function precision(precision?: null): number;
    function precision(precision?: number | null) {
        if (isNil(precision)) return lower48.precision();

        lower48.precision(precision);
        alaska.precision(precision);
        hawaii.precision(precision);
        puertoRico.precision(precision);
        guamMariana.precision(precision);
        americanSamoa.precision(precision);

        return reset();
    }

    function scale(scale?: number): ReturnType<typeof reset>;
    function scale(scale?: null): number;
    function scale(scale?: number | null) {
        if (isNil(scale)) return lower48.scale();

        lower48.scale(scale);
        alaska.scale(scale * 0.35);
        hawaii.scale(scale);
        puertoRico.scale(scale);
        guamMariana.scale(scale);
        americanSamoa.scale(scale);

        return albersUsaTerritories.translate(lower48.translate());
    }

    function translate(translate: [number, number]): ReturnType<typeof reset>;
    function translate(translate?: null): [number, number];
    function translate(translate?: [number, number] | null) {
        if (isNil(translate)) return lower48.translate();

        const k = lower48.scale()
        const x = +translate[0], y = +translate[1];

        lower48Point = lower48
            .translate(translate)
            .clipExtent([[x - 0.455 * k, y - 0.238 * k], [x + 0.455 * k, y + 0.238 * k]])
            .stream(pointStream);

        alaskaPoint = alaska
            .translate([x - 0.275 * k, y + 0.201 * k])
            .clipExtent([[x - 0.390 * k + epsilon, y + 0.120 * k + epsilon], [x - 0.185 * k - epsilon, y + 0.234 * k - epsilon]])
            .stream(pointStream);

        hawaiiPoint = hawaii
            .translate([x - 0.180 * k, y + 0.212 * k])
            .clipExtent([[x - 0.185 * k + epsilon, y + 0.166 * k + epsilon], [x - 0.080 * k - epsilon, y + 0.234 * k - epsilon]])
            .stream(pointStream);

        puertoRicoPoint = puertoRico
            .translate([x + 0.335 * k, y + 0.224 * k])
            .clipExtent([[x + 0.300 * k, y + 0.204 * k], [x + 0.380 * k, y + 0.234 * k]])
            .stream(pointStream);

        guamMarianaPoint = guamMariana
            .translate([x - 0.415 * k, y + 0.140 * k])
            .clipExtent([[x - 0.450 * k, y + 0.050 * k], [x - 0.390 * k, y + 0.210 * k]])
            .stream(pointStream);

        americanSamoaPoint = americanSamoa
            .translate([x - 0.415 * k, y + 0.215 * k])
            .clipExtent([[x - 0.450 * k, y + 0.210 * k], [x - 0.390 * k, y + 0.234 * k]])
            .stream(pointStream);

        return reset();
    }

    albersUsaTerritories.precision = precision;
    albersUsaTerritories.scale = scale;
    albersUsaTerritories.translate = translate;

    return albersUsaTerritories.scale(1070);
}

export default geoAlbersUsaTerritories;
