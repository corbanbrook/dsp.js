load('audio-harness.js');
load('dsp.js');

var iterations = 1000;

var envelope = new ADSR(0.01, 0.1, 0.5, 0.1, 0.2, 44100);

var calcADSR = function() {
  var fb     = getFramebuffer(),
      signal = DSP.getChannel(DSP.MIX, fb);

  envelope.process(signal);
};

runTest(calcADSR, iterations);
