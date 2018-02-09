var DSP = require("./dsp");

/**
 * WindowFunction
 *
 * @constructor
 * @param {Number} type
 * @param {Number} alpha
 */
function WindowFunction(type, alpha) {
  this.alpha = alpha;

  switch(type) {
    case DSP.BARTLETT:
      this.func = WindowFunction.Bartlett;
      break;

    case DSP.BARTLETTHANN:
      this.func = WindowFunction.BartlettHann;
      break;

    case DSP.BLACKMAN:
      this.func = WindowFunction.Blackman;
      this.alpha = this.alpha || 0.16;
      break;

    case DSP.COSINE:
      this.func = WindowFunction.Cosine;
      break;

    case DSP.GAUSS:
      this.func = WindowFunction.Gauss;
      this.alpha = this.alpha || 0.25;
      break;

    case DSP.HAMMING:
      this.func = WindowFunction.Hamming;
      break;

    case DSP.HANN:
      this.func = WindowFunction.Hann;
      break;

    case DSP.LANCZOS:
      this.func = WindowFunction.Lanczoz;
      break;

    case DSP.RECTANGULAR:
      this.func = WindowFunction.Rectangular;
      break;

    case DSP.TRIANGULAR:
      this.func = WindowFunction.Triangular;
      break;
  }
}

/**
 * Process a buffer
 * @param {Array} buffer
 */
WindowFunction.prototype.process = function(buffer) {
  var length = buffer.length;
  for ( var i = 0; i < length; i++ ) {
    buffer[i] *= this.func(length, i, this.alpha);
  }
  return buffer;
};

WindowFunction.Bartlett = function(length, index) {
  return 2 / (length - 1) * ((length - 1) / 2 - Math.abs(index - (length - 1) / 2));
};

WindowFunction.BartlettHann = function(length, index) {
  return 0.62 - 0.48 * Math.abs(index / (length - 1) - 0.5) - 0.38 * Math.cos(DSP.TWO_PI * index / (length - 1));
};

WindowFunction.Blackman = function(length, index, alpha) {
  var a0 = (1 - alpha) / 2;
  var a1 = 0.5;
  var a2 = alpha / 2;

  return a0 - a1 * Math.cos(DSP.TWO_PI * index / (length - 1)) + a2 * Math.cos(4 * Math.PI * index / (length - 1));
};

WindowFunction.Cosine = function(length, index) {
  return Math.cos(Math.PI * index / (length - 1) - Math.PI / 2);
};

WindowFunction.Gauss = function(length, index, alpha) {
  return Math.pow(Math.E, -0.5 * Math.pow((index - (length - 1) / 2) / (alpha * (length - 1) / 2), 2));
};

WindowFunction.Hamming = function(length, index) {
  return 0.54 - 0.46 * Math.cos(DSP.TWO_PI * index / (length - 1));
};

WindowFunction.Hann = function(length, index) {
  return 0.5 * (1 - Math.cos(DSP.TWO_PI * index / (length - 1)));
};

WindowFunction.Lanczos = function(length, index) {
  var x = 2 * index / (length - 1) - 1;
  return Math.sin(Math.PI * x) / (Math.PI * x);
};

WindowFunction.Rectangular = function() {
  return 1;
};

WindowFunction.Triangular = function(length, index) {
  return 2 / length * (length / 2 - Math.abs(index - (length - 1) / 2));
};

module.exports = WindowFunction;
