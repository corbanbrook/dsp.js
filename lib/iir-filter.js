var DSP = require("./dsp");
var ADSR = require("./adsr");

/**
 * IIRFilter
 * @constructor
 */
function IIRFilter(type, cutoff, resonance, sampleRate) {
  this.sampleRate = sampleRate;

  switch(type) {
    case DSP.LOWPASS:
    case DSP.LP12:
      this.func = new IIRFilter.LP12(cutoff, resonance, sampleRate);
      break;
  }
}

IIRFilter.prototype.__defineGetter__("cutoff",
  function() {
    return this.func.cutoff;
  }
);

IIRFilter.prototype.__defineGetter__("resonance",
  function() {
    return this.func.resonance;
  }
);

/**
 * Set filter parameters
 * @param {Number} cutoff
 * @param {Number} resonance
 */
IIRFilter.prototype.set = function(cutoff, resonance) {
  this.func.calcCoeff(cutoff, resonance);
};

/**
 * Process a buffer
 * @param {Array} buffer
 */
IIRFilter.prototype.process = function(buffer) {
  this.func.process(buffer);
};

/**
 * Add an envelope to the filter
 * @param {ADSR} envelope
 */
IIRFilter.prototype.addEnvelope = function(envelope) {
  if ( envelope instanceof ADSR ) {
    this.func.addEnvelope(envelope);
  } else {
    throw "Not an envelope.";
  }
};

/**
 * LP12 filter
 * @constructor
 */
IIRFilter.LP12 = function(cutoff, resonance, sampleRate) {
  this.sampleRate = sampleRate;
  this.vibraPos   = 0;
  this.vibraSpeed = 0;
  this.envelope = false;

  this.calcCoeff = function(cutoff, resonance) {
    this.w = 2.0 * Math.PI * cutoff / this.sampleRate;
    this.q = 1.0 - this.w / (2.0 * (resonance + 0.5 / (1.0 + this.w)) + this.w - 2.0);
    this.r = this.q * this.q;
    this.c = this.r + 1.0 - 2.0 * Math.cos(this.w) * this.q;

    this.cutoff = cutoff;
    this.resonance = resonance;
  };

  this.calcCoeff(cutoff, resonance);

  this.process = function(buffer) {
    for ( var i = 0; i < buffer.length; i++ ) {
      this.vibraSpeed += (buffer[i] - this.vibraPos) * this.c;
      this.vibraPos   += this.vibraSpeed;
      this.vibraSpeed *= this.r;

      /*
      var temp = this.vibraPos;

      if ( temp > 1.0 ) {
        temp = 1.0;
      } else if ( temp < -1.0 ) {
        temp = -1.0;
      } else if ( temp != temp ) {
        temp = 1;
      }

      buffer[i] = temp;
      */

      if (this.envelope) {
        buffer[i] = (buffer[i] * (1 - this.envelope.value())) + (this.vibraPos * this.envelope.value());
        this.envelope.samplesProcessed++;
      } else {
        buffer[i] = this.vibraPos;
      }
    }
  };
};

/**
 * Add an envelope to the filter
 * @param {ADSR} envelope
 */
IIRFilter.LP12.prototype.addEnvelope = function(envelope) {
  this.envelope = envelope;
};

module.exports = IIRFilter;
