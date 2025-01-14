# Dirty Reprojectors

Quick and dirty re-projections to trick your web maps out of web mercator.

## Install

    For command line use:

    npm install -g @opensnack/dirty-reprojectors

    For browser/node API use:

    npm install @opensnack/dirty-reprojectors

## Usage

### CLI

    cat input.geojson | dirty-reproject --forward PROJECTION [--reverse PROJECTION=mercator] > output.geojson

Example: to reproject some geojson so that web mapping libraries will render it
looking like 'albersUsa':

    cat input.geojson | dirty-reproject --forward albersUsa > output.geojson

For a list of supported projections, `dirty-reproject --list`


### API

`import { reproject, projections } from @opensnack/dirty-reprojectors;`

#### reproject

Reprojects the given geometry coordinate array _in place_, with
unprojectable points or degenerate geometries removed. If both
`options.forward` and `options.reverse` are supplied, then `forward` is
performed first.

**Parameters**

-   `options` **[Object](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Object)** 
    -   `options.forward` **([Function](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Statements/function) \| [string](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/String))?** The forward projection to use.
    -   `options.reverse` **([Function](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Statements/function) \| [string](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/String))?** The reverse projection to use.
    -   `options.projections` **[Object](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Object)?** A map of named projections to use.  If provided, then string values of `options.forward` or `options.reverse` will be used as keys to look up the projection function in `options.projections`.  For an extensive list provided by d3-geo-projection, use `import { projections }`.
-   `coordinates` **[Array](https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Array)** 

## How it works

Take, for example:

    cat input.geojson | dirty-reproject --forward albersUsa > output.geojson

What this actually does is:

1. Project `input.geojson` from WGS 84 (longitude/latitude) into `albersUsa`, with the target coordinates scaled to match the dimensions of Web Mercator.
2. Reverse-project the result back to WGS84 _as if_ it had been projected with Web Mercator.  So now, when your favorite web mapping library tries to project it into mercator, the geometries end up looking like they were projected using Albers.

The main catch is that if you actually look at the longitude/latitude
coordinates in `output.geojson`, they are totally wrong.  (There are other,
subtler catches, too, having to do with Web Mercator's limited latitude range,
varying loss of precision, and probably many other nuances I am not aware of.)

## Credits

 - Inspired by [this trick](https://www.mapbox.com/blog/mapping-arctic-ice-polar-projection/)
 - All the heavy lifting here is thanks to Mike Bostock's excellent [d3-geo](https://github.com/d3/d3-geo) and [de-geo-projection](https://github.com/d3/d3-geo-projection)

