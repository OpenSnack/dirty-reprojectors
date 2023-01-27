import { geoAlbers, geoConicEqualArea, type GeoStream } from 'd3-geo';
import { multiplex, isNil } from './helpers';

function geoAlbersUsaPr() {
    const epsilon = 0.000001;
    let cache: GeoStream | null,
        cacheStream: GeoStream | null,
        lower48Point: GeoStream,
        alaskaPoint: GeoStream,
        hawaiiPoint: GeoStream,
        puertoRicoPoint: GeoStream,
        point;

    const lower48 = geoAlbers(),
        alaska = geoConicEqualArea().rotate([154, 0]).center([-2, 58.5]).parallels([55, 65]),
        hawaii = geoConicEqualArea().rotate([157, 0]).center([-3, 19.9]).parallels([8, 18]),
        puertoRico = geoConicEqualArea().rotate([66, 0]).center([0, 18]).parallels([8, 18]),
        pointStream = { point: function(x: number, y: number) { point = [x, y]; } } as GeoStream;

    function albersUsaPr(coordinates: [number, number]): [number, number] {
        const x = coordinates[0]
        const y = coordinates[1];
        return point = null,
        (lower48Point.point(x, y), point)
            || (alaskaPoint.point(x, y), point)
            || (hawaiiPoint.point(x, y), point)
            || (puertoRicoPoint.point(x, y), point);
    }

    function reset() {
        cache = cacheStream = null;
        return albersUsaPr;
    }

    albersUsaPr.invert = function(coordinates: [number, number]) {
        const k = lower48.scale(),
            t = lower48.translate(),
            x = (coordinates[0] - t[0]) / k,
            y = (coordinates[1] - t[1]) / k;
        return (y >= 0.120 && y < 0.234 && x >= -0.425 && x < -0.214 ? alaska
            : y >= 0.166 && y < 0.234 && x >= -0.214 && x < -0.115 ? hawaii
                : y >= 0.204 && y < 0.234 && x >= 0.320 && x < 0.380 ? puertoRico
                    : lower48).invert?.(coordinates) ?? null;
    };

    albersUsaPr.stream = function(stream: GeoStream) {
        return cache && cacheStream === stream
            ? cache
            : cache = multiplex([
                lower48.stream(cacheStream = stream),
                alaska.stream(stream),
                hawaii.stream(stream),
                puertoRico.stream(stream)
            ]);
    };

    function precision(precision: number): ReturnType<typeof reset>;
    function precision(precision?: null): number;
    function precision(precision?: number | null) {
        if (isNil(precision)) return lower48.precision();

        lower48.precision(precision), alaska.precision(precision), hawaii.precision(precision), puertoRico.precision(precision);
        return reset();
    }

    function scale(scale: number): ReturnType<typeof reset>;
    function scale(scale?: null): number;
    function scale(scale?: number | null) {
        if (isNil(scale)) return lower48.scale();

        lower48.scale(scale), alaska.scale(scale * 0.35), hawaii.scale(scale), puertoRico.scale(scale);
        return albersUsaPr.translate(lower48.translate());
    }

    function translate(translate: [number, number]): ReturnType<typeof reset>;
    function translate(translate?: null): [number, number];
    function translate(translate?: [number, number] | null) {
        if (isNil(translate)) return lower48.translate();

        const k = lower48.scale(), x = +translate[0], y = +translate[1];

        lower48Point = lower48
            .translate(translate)
            .clipExtent([[x - 0.455 * k, y - 0.238 * k], [x + 0.455 * k, y + 0.238 * k]])
            .stream(pointStream);

        alaskaPoint = alaska
            .translate([x - 0.307 * k, y + 0.201 * k])
            .clipExtent([[x - 0.425 * k + epsilon, y + 0.120 * k + epsilon], [x - 0.214 * k - epsilon, y + 0.234 * k - epsilon]])
            .stream(pointStream);

        hawaiiPoint = hawaii
            .translate([x - 0.205 * k, y + 0.212 * k])
            .clipExtent([[x - 0.214 * k + epsilon, y + 0.166 * k + epsilon], [x - 0.115 * k - epsilon, y + 0.234 * k - epsilon]])
            .stream(pointStream);

        puertoRicoPoint = puertoRico
            .translate([x + 0.350 * k, y + 0.224 * k])
            .clipExtent([[x + 0.320 * k, y + 0.204 * k], [x + 0.380 * k, y + 0.234 * k]])
            .stream(pointStream);

        return reset();
    }

    albersUsaPr.precision = precision;
    albersUsaPr.scale = scale;
    albersUsaPr.translate = translate;

    return albersUsaPr.scale(1070);
}

export default geoAlbersUsaPr;
