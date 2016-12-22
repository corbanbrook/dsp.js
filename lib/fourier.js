/* global Float64Array */

/**
 * A Mixin for every fourier module
 * Fourier Transform Module used by DFT, FFT, RFFT
 * @private
 */
module.exports = function FourierTransform(bufferSize, sampleRate) {
  this.bufferSize = bufferSize;
  this.sampleRate = sampleRate;
  this.bandwidth  = 2 / bufferSize * sampleRate / 2;

  this.spectrum   = new Float64Array(bufferSize/2);
  this.real       = new Float64Array(bufferSize);
  this.imag       = new Float64Array(bufferSize);

  this.peakBand   = 0;
  this.peak       = 0;

  /**
   * Calculates the *middle* frequency of an FFT band.
   *
   * @param {Number} index The index of the FFT band.
   *
   * @returns The middle frequency in Hz.
   */
  this.getBandFrequency = function(index) {
    return this.bandwidth * index + this.bandwidth / 2;
  };

  /**
   * Calculate the spectrum (amplitude magnitudes)
   */
  this.calculateSpectrum = function() {
    var spectrum  = this.spectrum,
      real      = this.real,
      imag      = this.imag,
      bSi       = 2 / this.bufferSize,
      sqrt      = Math.sqrt,
      rval,
      ival,
      mag;

    for (var i = 0, N = bufferSize/2; i < N; i++) {
      rval = real[i];
      ival = imag[i];
      mag = bSi * sqrt(rval * rval + ival * ival);

      if (mag > this.peak) {
        this.peakBand = i;
        this.peak = mag;
      }

      spectrum[i] = mag;
    }
  };
};
