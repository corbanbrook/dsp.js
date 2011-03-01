var bufferSize = 2048;
var sampleRate = 44100;

var buffer1 = new Float32Array(bufferSize);
var buffer2 = new Float32Array(bufferSize);
var buffer3 = new Float32Array(bufferSize);
var buffer4 = new Float32Array(bufferSize);

for (var i = 0; i < bufferSize; i++) {
  buffer1[i] = (i % 2 === 0) ? -Math.random() : Math.random();
}

for (var i = 0; i < bufferSize; i++) {
  buffer2[i] = (i % 2 === 0) ? -Math.random() : Math.random();
}

for (var i = 0; i < bufferSize; i++) {
  buffer3[i] = (i % 2 === 0) ? -Math.random() : Math.random();
}

for (var i = 0; i < bufferSize; i++) {
  buffer4[i] = (i % 2 === 0) ? -Math.random() : Math.random();
}

var channel;
var temp;

var duration = benchmark(function() { 
  channel = DSP.deinterleave(DSP.MIX, buffer1);

  // cycle buffers
  temp = buffer1;
  buffer1 = buffer2;
  buffer2 = buffer3;
  buffer3 = buffer4;
  buffer4 = temp;
}, 100000);

print("Channel length: " + channel.length);
print("100000 iterations: " + (duration) + " ms (" + ((duration) / 100000) + "ms per iter)\n");
