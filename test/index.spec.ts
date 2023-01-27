import { describe, it, expect } from 'vitest';
import * as fs from 'node:fs';
import path from 'node:path';
import { reproject } from '../';
import projections from '../projections';

describe('dirty-reprojectors', () => {
    // disabled for now because numbers are very very very slightly different.....
    it.skip('each named projection', () => {
        for (const name in projections) {
            console.log(name);
            const inputFile = path.join(__dirname, './fixtures/utah.input.geojson')
            const outputFile = path.join(__dirname, `./fixtures/utah.output.${name}.geojson`)
    
            const input = JSON.parse(fs.readFileSync(inputFile).toString());
            const result = reproject({
                forward: name,
                reverse: 'mercator',
                projections: projections
            }, input.geometry)
    
            const output = Object.assign({}, input, { geometry: result });
            const expected = JSON.parse(fs.readFileSync(outputFile).toString());
    
            expect(output).toEqual(expected);
        }
    });

    it('dateline wrapping', () => {
        const input = JSON.parse(fs.readFileSync(path.join(__dirname, './fixtures/wrap.input.geojson')).toString());

        const result = reproject({
            forward: 'mercator',
            reverse: 'mercator',
            projections: projections
        }, input.geometry);

        const coordinates = result.coordinates[0];

        for (let i = 0; i < coordinates.length - 1; i++) {
            const delta = Math.abs(coordinates[i][0] - coordinates[i + 1][0])
            expect(delta).toBeLessThanOrEqual(90);
        }
    });
});

// test('dateline wrapping', function (t) {
//     const input = JSON.parse(fs.readFileSync(path.join(__dirname, './fixtures/wrap.input.geojson')))
//     const result = reproject({
//         forward: 'mercator',
//         reverse: 'mercator',
//         projections: projections
//     }, input.geometry)

//     const coordinates = result.coordinates[0]

//     for (let i = 0; i < coordinates.length - 1; i++) {
//         const delta = Math.abs(coordinates[i][0] - coordinates[i + 1][0])
//         t.ok(delta <= 90, 'coordinates are not weirdly wrapped after fwd+reverse projection')
//     }

//     t.end()
// })

// test('each named projection', {skip: true}, function (t) {
//   for (const name in projections) {
//     const inputFile = path.join(__dirname, './fixtures/utah.input.geojson')
//     const outputFile = path.join(__dirname, `./fixtures/utah.output.${name}.geojson`)

//     t.test(name, function (t) {
//       const input = JSON.parse(fs.readFileSync(inputFile))
//       const result = reproject({
//         forward: name,
//         reverse: 'mercator',
//         projections: projections
//       }, input.geometry)

//       const output = Object.assign({}, input, { geometry: result })

//       if (process.env.UPDATE) {
//         fs.writeFileSync(outputFile, JSON.stringify(output, null, 2))
//       }
//       const expected = JSON.parse(fs.readFileSync(outputFile))

//       t.same(output, expected)
//       t.end()
//     })
//   }

//   t.end()
// })
