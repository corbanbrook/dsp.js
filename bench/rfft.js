var bufferSize = 2048;
var sampleRate = 44100;
var frequency = FREQUENCY || 440;

var fft = new RFFT(bufferSize, sampleRate);
var osc = new Oscillator(DSP.SAW, frequency, 1.0, bufferSize, sampleRate);
var signal = osc.generate();

var duration = benchmark(function() { fft.forward(signal); });

var peakBand = 0;

for (var i = 0; i < fft.spectrum.length; i++) {
  peakBand = (fft.spectrum[i] > fft.spectrum[peakBand]) ? i : peakBand;
}

var peakFreq = fft.getBandFrequency(fft.peakBand);

print("Detected peak: " + peakFreq + " Hz (error " + Math.abs(peakFreq - frequency) + " Hz)");
print("10000 FFTs: " + (duration) + " ms (" + ((duration) / 10000) + "ms per FFT)\n");
