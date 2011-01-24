var bufferSize = 2048;
var sampleRate = 44100;
var frequency = FREQUENCY || 440;

var dft = new DFT(bufferSize, sampleRate);
var osc = new Oscillator(DSP.SAW, frequency, 1.0, bufferSize, sampleRate);
var signal = osc.generate();

var duration = benchmark(function() { dft.forward(signal); }, 20);

var peakBand = 0;

for (var i = 0; i < dft.spectrum.length; i++) {
  peakBand = (dft.spectrum[i] > dft.spectrum[peakBand]) ? i : peakBand;
}

var peakFreq = dft.getBandFrequency(dft.peakBand);

print("Detected peak: " + peakFreq + " Hz (error " + Math.abs(peakFreq - frequency) + " Hz)");
print("20 DFTs: " + (duration) + " ms (" + ((duration) / 20) + "ms per DFT)\n");
