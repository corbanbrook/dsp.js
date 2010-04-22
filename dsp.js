/*  
 *  DSP.js - a comprehensive digital signal processing  library for javascript 
 *  
 *  Created by Corban Brook <corbanbrook@gmail.com> on 2010-01-01.
 *  Copyright 2010 Corban Brook. All rights reserved.
 *
 */

DSP = function() {};

// Waveforms
DSP.SINEWAVE     = 1;
DSP.SQUAREWAVE   = 2;
DSP.SAWWAVE      = 3;
DSP.TRIANGLEWAVE = 4;

// IIR Filters
DSP.LOWPASS      = 5;
DSP.HIGHPASS     = 6;

DSP.TWO_PI = 2*Math.PI;

DFT = function(bufferSize, sampleRate) {
  this.bufferSize = bufferSize;
  this.sampleRate = sampleRate;

  var N = bufferSize/2 * bufferSize;
      
  this.sinTable = new Array(N);
  this.cosTable = new Array(N);
  
  for ( var i = 0; i < N; i++ ) {
    this.sinTable[i] = Math.sin(i * DSP.TWO_PI / bufferSize);
    this.cosTable[i] = Math.cos(i * DSP.TWO_PI / bufferSize);
  }
  
  this.spectrum = new Array(bufferSize/2);
  this.complexValues = new Array(bufferSize/2);
};

DFT.prototype.forward = function(buffer) {
  var real, imag;

  for ( var k = 0; k < this.bufferSize/2; k++ ) {
    real = 0.0;
    imag = 0.0;

    for ( var n = 0; n < buffer.length; n++ ) {
      real += this.cosTable[k*n] * signal[n];
      imag += this.sinTable[k*n] * signal[n];
    }

    this.complexValues[k] = {real: real, imag: imag};
  }
  
  for ( var i = 0; i < this.bufferSize/2; i++ ) {
    this.spectrum[i] = 2 * Math.sqrt(Math.pow(this.complexValues[i].real, 2) + Math.pow(this.complexValues[i].imag, 2)) / this.bufferSize;
  }

  return this.spectrum;
};

FFT = function(bufferSize, sampleRate) {
  this.bufferSize = bufferSize;
  this.sampleRate = sampleRate;

  this.spectrum      = new Array(bufferSize/2);
  this.complexValues = new Array(bufferSize);
    
  this.reverseTable  = new Array(bufferSize);
  this.reverseTable[0] = 0;

  var limit = 1;
  var bit = bufferSize >> 1;

  while ( limit < bufferSize ) {
    for ( var i = 0; i < limit; i++ ) {
      this.reverseTable[i + limit] = this.reverseTable[i] + bit;
    }

    limit = limit << 1;
    bit = bit >> 1;
  }
    
  this.sinTable = new Array(bufferSize);
  this.cosTable = new Array(bufferSize);
  
  for ( var i = 0; i < bufferSize; i++ ) {
    this.sinTable[i] = Math.sin(-Math.PI/i);
    this.cosTable[i] = Math.cos(-Math.PI/i);
  }
};

FFT.prototype.forward = function(buffer) {
  if ( this.bufferSize % 2 != 0 )         { throw "Invalid buffer size, must be a power of 2."; }
  if ( this.bufferSize != buffer.length ) { throw "Supplied buffer is not the same size as defined FFT. FFT Size: " + this.bufferSize + " Buffer Size: " + buffer.length; }

  for ( var i = 0; i < buffer.length; i++ ) {
    this.complexValues[i] = {real: buffer[this.reverseTable[i]], imag: 0.0};
  }

  var halfSize = 1;

  while ( halfSize < buffer.length ) {
    var phaseShiftStepReal = this.cosTable[halfSize];
    var phaseShiftStepImag = this.sinTable[halfSize];
    var currentPhaseShiftReal = 1.0;
    var currentPhaseShiftImag = 0.0;

    for ( var fftStep = 0; fftStep < halfSize; fftStep++ ) {
      var i = fftStep;

      while ( i < buffer.length ) {
        var off = i + halfSize;
        var tr = (currentPhaseShiftReal * this.complexValues[off].real) - (currentPhaseShiftImag * this.complexValues[off].imag);
        var ti = (currentPhaseShiftReal * this.complexValues[off].imag) + (currentPhaseShiftImag * this.complexValues[off].real);

        this.complexValues[off].real = this.complexValues[i].real - tr;
        this.complexValues[off].imag = this.complexValues[i].imag - ti;
        this.complexValues[i].real += tr;
        this.complexValues[i].imag += ti;

        i += halfSize << 1;
      }

      var tmpReal = currentPhaseShiftReal;
      currentPhaseShiftReal = (tmpReal * phaseShiftStepReal) - (currentPhaseShiftImag * phaseShiftStepImag);
      currentPhaseShiftImag = (tmpReal * phaseShiftStepImag) + (currentPhaseShiftImag * phaseShiftStepReal);
    }

    halfSize = halfSize << 1;
  }

  for ( var i = 0; i < this.bufferSize/2; i++ ) {
    this.spectrum[i] = 2 * Math.sqrt(Math.pow(this.complexValues[i].real, 2) + Math.pow(this.complexValues[i].imag, 2)) / this.bufferSize;
  }
  
  return this.spectrum;
};
  
/*  Oscillator Signal Generator
 *    
 *  Usage: var sine = Oscillator(SINEWAVE, 440.0, 1, 2048, 44100);
 *         var signal = sine.generate();
 *
 */

Oscillator = function(waveFunc, frequency, amplitude, bufferSize, sampleRate) {
    this.waveFunc   = waveFunc;
    this.frequency  = frequency;
    this.amplitude  = amplitude;
    this.bufferSize = bufferSize;
    this.sampleRate = sampleRate;
    this.frameCount = 0;
    
    this.cyclesPerSample = frequency / sampleRate;

    this.signal = new Array(bufferSize);
    this.envelope = null;
    this.envelopedSignal = new Array(bufferSize);

    this.generate();
}; 

Oscillator.prototype.setAmp = function(amplitude) {
  if (amplitude >= 0 && amplitude <= 1) {
    this.amplitude = amplitude;
    //this.generate();
  } else {
    throw "Amplitude out of range (0..1).";
  }
};
      
Oscillator.prototype.setFreq = function(frequency) {
  this.frequency = frequency;
  this.cyclesPerSample = frequency / this.sampleRate;
  //this.generate();
};
      
// Add an oscillator
Oscillator.prototype.add = function(oscillator) {
  for ( var i = 0; i < this.bufferSize; i++ ) {
    //self.signal[i] += oscillator.valueAt(i);
    this.signal[i] += oscillator.signal[i];
  }
  
  return this.signal;
};
      
// Add a signal to the current generated osc signal
Oscillator.prototype.addSignal = function(signal) {
  for ( var i = 0; i < signal.length; i++ ) {
    if ( i >= this.bufferSize ) {
      break;
    }
    this.signal[i] += signal[i];
  }
  return self.signal;
};
      
// Add an envelope to the oscillator
Oscillator.prototype.addEnvelope = function(envelope) {
  this.envelope = envelope;
};
      
Oscillator.prototype.valueAt = function(offset) {
  return this.waveLength[offset % this.waveLength.length];
};
      
Oscillator.prototype.generate = function() {
  var frameOffset = this.frameCount * this.bufferSize, step;

  for ( var i = 0; i < this.bufferSize; i++ ) {
    step = (frameOffset + i) * this.cyclesPerSample % 1;

    this.signal[i] = this.waveFunc(step) * this.amplitude;
  }

  this.frameCount++;

  return this.signal;
};

Oscillator.Sine = function(step) {
  return Math.sin(DSP.TWO_PI * step);
};

Oscillator.Square = function(step) {
  return step < 0.5 ? 1 : -1;
};

Oscillator.Saw = function(step) {
  return 2 * (step - Math.round(step));
};

Oscillator.Triangle = function(step) {
  return 1 - 4 * Math.abs(Math.round(step) - step);
};
  
ADSR = function(attackLength, decayLength, sustainLevel, sustainLength, releaseLength, sampleRate) {
  this.attackLength  = attackLength;
  this.decayLength   = decayLength;
  this.sustainLevel  = sustainLevel;
  this.sustainLength = sustainLength;
  this.releaseLength = releaseLength;
  this.sampleRate    = sampleRate;
  
  this.attackSamples  = attackLength  * sampleRate;
  this.decaySamples   = decayLength   * sampleRate;
  this.sustainSamples = sustainLength * sampleRate;
  this.releaseSamples = releaseLength * sampleRate;
  
  this.attack         =                this.attackSamples;
  this.decay          = this.attack  + this.decaySamples;
  this.sustain        = this.decay   + this.sustainSamples;
  this.release        = this.sustain + this.releaseSamples;
    
  this.samplesProcessed = 0;
};


ADSR.prototype.trigger = function() {
  this.samplesProcessed = 0;
};

ADSR.prototype.processSample = function(sample) {
  var amplitude = 0;

  if ( this.samplesProcessed <= this.attack ) {
    amplitude = 0 + (1 - 0) * ((this.samplesProcessed - 0) / (this.attack - 0));
  } else if ( this.samplesProcessed > this.attack && this.samplesProcessed <= this.decay ) {
    amplitude = 1 + (this.sustainLevel - 1) * ((this.samplesProcessed - this.attack) / (this.decay - this.attack));
  } else if ( this.samplesProcessed > this.decay && this.samplesProcessed <= this.sustain ) {
    amplitude = this.sustainLevel;
  } else if ( this.samplesProcessed > this.sustain && this.samplesProcessed <= this.release ) {
    amplitude = this.sustainLevel + (0 - this.sustainLevel) * ((this.samplesProcessed - this.sustain) / (this.release - this.sustain));
  }
  
  return sample * amplitude;
};
      
ADSR.prototype.process = function(buffer) {
  for ( var i = 0; i < buffer.length; i++ ) {
    /*
    var amplitude = 0;
    
    if ( this.samplesProcessed <= this.attack ) {
      amplitude = 0 + (1 - 0) * ((this.samplesProcessed - 0) / (this.attack - 0));
    } else if ( this.samplesProcessed > this.attack && this.samplesProcessed <= this.decay ) {
      amplitude = 1 + (this.sustainLevel - 1) * ((this.samplesProcessed - this.attack) / (this.decay - this.attack));
    } else if ( this.samplesProcessed > this.decay && this.samplesProcessed <= this.sustain ) {
      amplitude = this.sustainLevel;
    } else if ( this.samplesProcessed > this.sustain && this.samplesProcessed <= this.release ) {
      amplitude = this.sustainLevel + (0 - this.sustainLevel) * ((this.samplesProcessed - this.sustain) / (this.release - this.sustain));
    }
    
    buffer[i] *= amplitude;

    this.samplesProcessed++;
    */

    buffer[i] = this.processSample(buffer[i]);

    this.samplesProcessed++;
  }
  
  return buffer;
};
      
      
ADSR.prototype.isActive = function() {
  if ( this.samplesProcessed > this.release ) {
    return false;
  } else {
    return true;
  }
};
  
LP12 = function(cutoff, resonance, sampleRate) {
  this.sampleRate = sampleRate;
  this.vibraPos   = 0; 
  this.vibraSpeed = 0;
  
  this.calcCoeff = function(cutoff, resonance) {
    this.w = 2.0 * Math.PI * cutoff / this.sampleRate;
    this.q = 1.0 - this.w / (2.0 * (resonance + 0.5 / (1.0 + this.w)) + this.w - 2.0);
    this.r = this.q * this.q;
    this.c = this.r + 1.0 - 2.0 * Math.cos(this.w) * this.q;
    
    this.cutoff = cutoff;
    this.resonance = resonance;
  };

  this.calcCoeff(cutoff, resonance);
};  

LP12.prototype.set = function(cutoff, resonance) {
  this.calcCoeff(cutoff, resonance);
};

LP12.prototype.process = function(buffer) {
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
    buffer[i] = this.vibraPos;
  }
};

WindowFunction = function(windowFunc) {
  this.windowFunc = windowFunc;
};

WindowFunction.prototype.process = function(buffer) {
  var length = buffer.length;
  for ( var i = 0; i < length; i++ ) {
    buffer[i] *= this.windowFunc(i, length);
  }
};

WindowFunction.Hann = function(index, length) {
  return 0.5 * (1 - Math.cos(DSP.TWO_PI * index / (length - 1)));
};

