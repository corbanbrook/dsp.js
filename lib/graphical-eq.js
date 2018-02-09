/* global Float64Array */
var DSP = require("./dsp");
var Biquad = require("./biquad");

/**
 *  Create a Graphical Equalizer
 *
 *  Implementation of a graphic equalizer with a configurable bands-per-octave
 *  and minimum and maximum frequencies
 *
 *  Created by Ricard Marxer <email@ricardmarxer.com> on 2010-05-23.
 *  Copyright 2010 Ricard Marxer. All rights reserved.
 *
 * @constructor
 * @param {SampleRate}
 * @example
 * var eq = new GraphicalEq(44100)
 */
function GraphicalEq(sampleRate) {
  this.FS = sampleRate;
  this.minFreq = 40.0;
  this.maxFreq = 16000.0;

  this.bandsPerOctave = 1.0;

  this.filters = [];
  this.freqzs = [];

  this.calculateFreqzs = true;

  this.recalculateFilters = function() {
    var bandCount = Math.round(Math.log(this.maxFreq/this.minFreq) * this.bandsPerOctave/ Math.LN2);

    this.filters = [];
    for (var i=0; i<bandCount; i++) {
      var freq = this.minFreq*(Math.pow(2, i/this.bandsPerOctave));
      var newFilter = new Biquad(DSP.PEAKING_EQ, this.FS);
      newFilter.setDbGain(0);
      newFilter.setBW(1/this.bandsPerOctave);
      newFilter.setF0(freq);
      this.filters[i] = newFilter;
      this.recalculateFreqz(i);
    }
  };

  this.setMinimumFrequency = function(freq) {
    this.minFreq = freq;
    this.recalculateFilters();
  };

  this.setMaximumFrequency = function(freq) {
    this.maxFreq = freq;
    this.recalculateFilters();
  };

  this.setBandsPerOctave = function(bands) {
    this.bandsPerOctave = bands;
    this.recalculateFilters();
  };

  this.setBandGain = function(bandIndex, gain) {
    if (bandIndex < 0 || bandIndex > (this.filters.length-1)) {
      throw "The band index of the graphical equalizer is out of bounds.";
    }

    if (!gain) {
      throw "A gain must be passed.";
    }

    this.filters[bandIndex].setDbGain(gain);
    this.recalculateFreqz(bandIndex);
  };

  this.recalculateFreqz = function(bandIndex) {
    if (!this.calculateFreqzs) {
      return;
    }

    if (bandIndex < 0 || bandIndex > (this.filters.length-1)) {
      throw "The band index of the graphical equalizer is out of bounds. " + bandIndex + " is out of [" + 0 + ", " + this.filters.length-1 + "]";
    }

    if (!this.w) {
      this.w = new Float64Array(400);
      for (var i=0; i<this.w.length; i++) {
         this.w[i] = Math.PI/this.w.length * i;
      }
    }

    var b = [this.filters[bandIndex].b0, this.filters[bandIndex].b1, this.filters[bandIndex].b2];
    var a = [this.filters[bandIndex].a0, this.filters[bandIndex].a1, this.filters[bandIndex].a2];

    this.freqzs[bandIndex] = DSP.mag2db(DSP.freqz(b, a, this.w));
  };

  this.process = function(buffer) {
    var output = buffer;

    for (var i = 0; i < this.filters.length; i++) {
      output = this.filters[i].process(output);
    }

    return output;
  };

  this.processStereo = function(buffer) {
    var output = buffer;

    for (var i = 0; i < this.filters.length; i++) {
      output = this.filters[i].processStereo(output);
    }

    return output;
  };
}

module.exports = GraphicalEq;
