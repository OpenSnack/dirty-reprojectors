#!/usr/bin/env node

import bbox from '@turf/bbox';
import geojsonStream from 'geojson-stream';
// import through from 'through2';
import { Transform } from 'stream';
import minimist from 'minimist';
import projections from './projections';
import { reproject } from '.';
  
function usage (message: string) {
    if (message) { console.error(message) }
    console.error('Usage:')
    console.error('cat normal.geojson | dirty-reproject [--forward PROJECTION] [--reverse PROJECTION=mercator] > weird.geojson')
    console.error('cat normal.geojson | dirty-reproject --forward PROJECTION --no-reverse > projected-space.geojson')
    console.error('dirty-reproject --list (for a list of supported projections)')
    process.exit(1)
}

const argv = minimist(process.argv.slice(2), {
    alias: {
        f: 'forward',
        r: 'reverse',
        l: 'list'
    },
    default: {
        reverse: 'mercator'
    }
});

if (argv.list) {
    console.log(Object.keys(projections).join('\n'));
    process.exit(0);
}

if (!argv.forward && !argv.reverse) {
    usage('Please supply either a --forward or --reverse projection (or both).');
}

process.stdin.pipe(geojsonStream.parse())
    .pipe(new Transform({
        transform(feature, _, next) {
            feature.geometry = reproject({
                forward: argv.forward,
                reverse: argv.reverse,
                projections
            }, feature.geometry);
            feature.bbox = bbox(feature);
            (this as any).push(feature)
            next();
        },
        objectMode: true
    }))
    .pipe(geojsonStream.stringify())
    .pipe(process.stdout);
