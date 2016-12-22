load('audio-harness.js');
load('dsp.js');

var iterations = 1000;

var filter = new IIRFilter(DSP.LOWPASS, 200, 44100);

var calcIIRFilter = function() {
  var fb     = getFramebuffer(),
      signal = DSP.getChannel(DSP.MIX, fb);

  filter.process(signal);
};

runTest(calcIIRFilter, iterations);
