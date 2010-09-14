load('audio-harness.js');
load('dsp.js');

var iterations = 100;
var dft = new DFT(frameBufferLength / channels, rate);

var calcDFT = function() {
  var fb     = getFramebuffer(),
      signal = DSP.getChannel(DSP.MIX, fb);

  dft.forward(signal);
};

runTest(calcDFT, iterations);
