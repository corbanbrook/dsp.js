## DSP.js [![npm](https://img.shields.io/npm/v/dsp.js.svg?style=flat-square)](https://www.npmjs.com/package/dsp.js)

DSP.js is a comprehensive digital signal processing library for javascript.
It includes many functions for signal analysis and generation, including
Oscillators(sine, saw, square, triangle), Window functions (Hann, Hamming, etc),
Envelopes(ADSR), IIR Filters(lowpass, highpass, bandpass, notch), FFT and DFT
transforms, Delays, Reverb.

## Install

__Node__

Via npm: `npm install --save dsp.js` and require the modules:

```js
const { DSP, WindowFunction, Oscillator } = require('dsp.js')
```

You can require individual modules:

```js
const DSP = require('dsp.js/lib/dsp')
const WindowFunction = require('dsp.js/lib/window-function')
```

__Browser__

Or grab the [dsp.js](https://github.com/corbanbrook/dsp.js/blob/master/dist/dsp.min.js) and include it in your html:

```html
<script src="dsp.min.js"></script>
```

## Usage

[Read the API documentation](https://github.com/corbanbrook/dsp.js/blob/master/docs/API.md)

## Tests and scripts

To setup you local machine to work with the code, first you have to clone this repository and install dependencies (`node` and `npm` are assumed to be installed): `npm install`

- Run the tests: `npm test`
- Generate API reference and distribution files: `npm run dist`
- Run benchmarks: `npm run bench`


## License

MIT License

Copyright (c) 2010 Corban Brook @corban                                                    
