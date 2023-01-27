const { commonjs } = require("@hyrious/esbuild-plugin-commonjs");
const esbuild = require('esbuild');

// browser
esbuild.build({
    entryPoints: ['index.ts'],
    bundle: true,
    format: 'esm',
    platform: 'browser',
    outfile: 'dist/reproject.mjs',
    plugins: [commonjs()],
});

// node
esbuild.build({
    entryPoints: ['index.ts'],
    bundle: true,
    format: 'esm',
    platform: 'node',
    outfile: 'dist/reproject.node.mjs',
    plugins: [commonjs()],
});

// bin
esbuild.build({
    entryPoints: ['bin.ts'],
    bundle: true,
    format: 'esm',
    platform: 'node',
    outfile: 'bin/reproject.bin.mjs',
    plugins: [commonjs()],
});
