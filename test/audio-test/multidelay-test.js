load('audio-harness.js');
load('dsp.js');

var iterations = 1000;

var delay = new MultiDelay(44100*5, 44100*1, 1.0, 0.6);

var calcDelay = function() {
  var fb     = getFramebuffer(),
      signal = DSP.getChannel(DSP.MIX, fb);

  delay.process(signal);
};

runTest(calcDelay, iterations);
