/* global Float64Array */
/**
 * MultiDelay effect by Almer Thie (http://code.almeros.com).
 * Copyright 2010 Almer Thie. All rights reserved.
 * Example: http://code.almeros.com/code-examples/delay-firefox-audio-api/
 *
 * This is a delay that feeds it's own delayed signal back into its circular
 * buffer. Also known as a CombFilter.
 *
 * Compatible with interleaved stereo (or more channel) buffers and
 * non-interleaved mono buffers.
 *
 * @param {Number} maxDelayInSamplesSize Maximum possible delay in samples (size of circular buffer)
 * @param {Number} delayInSamples Initial delay in samples
 * @param {Number} masterVolume Initial master volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 * @param {Number} delayVolume Initial feedback delay volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 *
 * @constructor
 */
function MultiDelay(maxDelayInSamplesSize, delayInSamples, masterVolume, delayVolume) {
  this.delayBufferSamples   = new Float64Array(maxDelayInSamplesSize); // The maximum size of delay
  this.delayInputPointer     = delayInSamples;
  this.delayOutputPointer   = 0;

  this.delayInSamples   = delayInSamples;
  this.masterVolume     = masterVolume;
  this.delayVolume     = delayVolume;
}

/**
 * Change the delay time in samples.
 *
 * @param {Number} delayInSamples Delay in samples
 */
MultiDelay.prototype.setDelayInSamples = function (delayInSamples) {
  this.delayInSamples = delayInSamples;

  this.delayInputPointer = this.delayOutputPointer + delayInSamples;

  if (this.delayInputPointer >= this.delayBufferSamples.length-1) {
    this.delayInputPointer = this.delayInputPointer - this.delayBufferSamples.length;
  }
};

/**
 * Change the master volume.
 *
 * @param {Number} masterVolume Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 */
MultiDelay.prototype.setMasterVolume = function(masterVolume) {
  this.masterVolume = masterVolume;
};

/**
 * Change the delay feedback volume.
 *
 * @param {Number} delayVolume Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify)
 */
MultiDelay.prototype.setDelayVolume = function(delayVolume) {
  this.delayVolume = delayVolume;
};

/**
 * Process a given interleaved or mono non-interleaved float value Array and adds the delayed audio.
 *
 * @param {Array} samples Array containing Float values or a Float64Array
 *
 * @returns A new Float64Array interleaved or mono non-interleaved as was fed to this function.
 */
MultiDelay.prototype.process = function(samples) {
  // NB. Make a copy to put in the output samples to return.
  var outputSamples = new Float64Array(samples.length);

  for (var i=0; i<samples.length; i++) {
    // delayBufferSamples could contain initial NULL's, return silence in that case
    var delaySample = (this.delayBufferSamples[this.delayOutputPointer] === null ? 0.0 : this.delayBufferSamples[this.delayOutputPointer]);

    // Mix normal audio data with delayed audio
    var sample = (delaySample * this.delayVolume) + samples[i];

    // Add audio data with the delay in the delay buffer
    this.delayBufferSamples[this.delayInputPointer] = sample;

    // Return the audio with delay mix
    outputSamples[i] = sample * this.masterVolume;

    // Manage circulair delay buffer pointers
    this.delayInputPointer++;
    if (this.delayInputPointer >= this.delayBufferSamples.length-1) {
      this.delayInputPointer = 0;
    }

    this.delayOutputPointer++;
    if (this.delayOutputPointer >= this.delayBufferSamples.length-1) {
      this.delayOutputPointer = 0;
    }
  }

  return outputSamples;
};

module.exports = MultiDelay;
