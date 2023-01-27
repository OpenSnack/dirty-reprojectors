import { GeoStream } from 'd3-geo';
declare function geoAlbersUsaTerritories(): {
    (coordinates: [number, number]): [number, number];
    invert(coordinates: [number, number]): [number, number] | null;
    stream(stream: GeoStream): GeoStream;
    precision: {
        (precision: number): any;
        (precision?: null): number;
    };
    scale: {
        (scale?: number): any;
        (scale?: null): number;
    };
    translate: {
        (translate: [number, number]): any;
        (translate?: null): [number, number];
    };
};
export default geoAlbersUsaTerritories;
