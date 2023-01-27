import type { GeoStream } from 'd3-geo';
export declare function multiplex(streams: GeoStream[]): {
    point(x: number, y: number): void;
    sphere(): void;
    lineStart(): void;
    lineEnd(): void;
    polygonStart(): void;
    polygonEnd(): void;
};
export declare function isNil(n: unknown): n is null | undefined;
