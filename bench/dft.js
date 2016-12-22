/* global FREQUENCY */
/* eslint-disable no-console */
var benchmark = require("./bench");
var dsp = require("..");
var bufferSize = 2048;
var sampleRate = 44100;
var frequency = FREQUENCY || 440;

var dft = new dsp.DFT(bufferSize, sampleRate);
var osc = new dsp.Oscillator(dsp.DSP.SAW, frequency, 1.0, bufferSize, sampleRate);
var signal = osc.generate();

var duration = benchmark(function() { dft.forward(signal); }, 20);

var peakBand = 0;

for (var i = 0; i < dft.spectrum.length; i++) {
  peakBand = (dft.spectrum[i] > dft.spectrum[peakBand]) ? i : peakBand;
}

var peakFreq = dft.getBandFrequency(dft.peakBand);

console.log("Detected peak: " + peakFreq + " Hz (error " + Math.abs(peakFreq - frequency) + " Hz)");
console.log("20 DFTs: " + (duration) + " ms (" + ((duration) / 20) + "ms per DFT)\n");
