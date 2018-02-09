"use strict";

/*
 *  DSP.js - a comprehensive digital signal processing  library for javascript
 *
 *  Created by Corban Brook <corbanbrook@gmail.com> on 2010-01-01.
 *  Copyright 2010 Corban Brook. All rights reserved.
 *
 */
var DSPJS = {
  DSP: require("./dsp"),
  DFT: require("./dft"),
  FFT: require("./fft"),
  RFFT: require("./rfft"),
  Sampler: require("./sampler"),
  Oscillator: require("./oscillator"),
  ADSR: require("./adsr"),
  IIRFilter: require("./iir-filter"),
  IIRFilter2: require("./iir-filter2"),
  WindowFunction: require("./window-function"),
  sinh: require("./sinh"),
  Biquad: require("./biquad"),
  GraphicalEq: require("./graphical-eq"),
  MultiDelay: require("./multi-delay"),
  SingleDelay: require("./single-delay"),
  Reverb: require("./reverb")
};

if (typeof module === "object" && module.exports) module.exports = DSPJS;
if (typeof window !== "undefined") {
  Object.keys(DSPJS).forEach(function (k) {
    window[k] = DSPJS[k];
  });
}
