load('audio-harness.js');
load('dsp.js');

var iterations = 1000;

var reverb = new Reverb(20000, 6500, 0.8, 0.5, 0.9, 4500);

var calcReverb = function() {
  var fb     = getFramebuffer(),
      signal = DSP.getChannel(DSP.MIX, fb);

  reverb.process(signal);
};

runTest(calcReverb, iterations);
