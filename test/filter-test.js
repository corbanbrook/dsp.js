load('audio-harness.js');
load('dsp.js');

var iterations = 1000;

var sampleRate = 44100;
var lp12 = new IIRFilter(DSP.LOWPASS, 22050, 0, sampleRate);
lp12.set(2500, 0.3);

//var frameSize = 2048;

var calcFilter = function() {
  var fb = getFramebuffer();
  var signal = DSP.getChannel(DSP.MIX, fb);
  lp12.process(signal);
};

runTest(calcFilter, iterations);
