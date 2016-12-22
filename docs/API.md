## Classes

<dl>
<dt><a href="#ADSR">ADSR</a></dt>
<dd></dd>
<dt><a href="#Biquad">Biquad</a></dt>
<dd></dd>
<dt><a href="#DFT">DFT</a></dt>
<dd></dd>
<dt><a href="#FFT">FFT</a></dt>
<dd></dd>
<dt><a href="#GraphicalEq">GraphicalEq</a></dt>
<dd></dd>
<dt><a href="#IIRFilter">IIRFilter</a></dt>
<dd></dd>
<dt><a href="#IIRFilter2">IIRFilter2</a></dt>
<dd></dd>
<dt><a href="#MultiDelay">MultiDelay</a></dt>
<dd></dd>
<dt><a href="#Oscillator">Oscillator</a></dt>
<dd></dd>
<dt><a href="#Reverb">Reverb</a></dt>
<dd></dd>
<dt><a href="#RFFT">RFFT</a></dt>
<dd></dd>
<dt><a href="#Sampler">Sampler</a></dt>
<dd></dd>
<dt><a href="#SingleDelay">SingleDelay</a></dt>
<dd></dd>
<dt><a href="#WindowFunction">WindowFunction</a></dt>
<dd></dd>
</dl>

## Members

<dl>
<dt><a href="#DSP">DSP</a></dt>
<dd><p>DSP is an object which contains general purpose utility functions and constants</p>
</dd>
</dl>

## Functions

<dl>
<dt><a href="#sinh">sinh(num)</a></dt>
<dd><p>Returns the hyperbolic sine of the number</p>
</dd>
</dl>

<a name="ADSR"></a>

## ADSR
**Kind**: global class  

* [ADSR](#ADSR)
    * [new ADSR(attack, decay, sustain, release, sampleRate)](#new_ADSR_new)
    * [.noteOn()](#ADSR+noteOn)
    * [.noteOff()](#ADSR+noteOff)
    * [.processSample()](#ADSR+processSample)
    * [.value()](#ADSR+value) ⇒ <code>Number</code>
    * [.process(buffer)](#ADSR+process)
    * [.isActive()](#ADSR+isActive) ⇒ <code>Boolean</code>
    * [.disable()](#ADSR+disable)

<a name="new_ADSR_new"></a>

### new ADSR(attack, decay, sustain, release, sampleRate)
ADSR Envelope


| Param | Type | Description |
| --- | --- | --- |
| attack | <code>Number</code> | The attack length in seconds |
| decay | <code>Number</code> | The decay length in seconds |
| sustain | <code>Number</code> | The sustain level |
| release | <code>Number</code> | The release length in seconds |
| sampleRate | <code>Number</code> | The the sample rate |

<a name="ADSR+noteOn"></a>

### adsR.noteOn()
Start the envelope

**Kind**: instance method of <code>[ADSR](#ADSR)</code>  
<a name="ADSR+noteOff"></a>

### adsR.noteOff()
Stop the envelope

Send a note off when using a sustain of infinity to let the envelope enter
the release phase

**Kind**: instance method of <code>[ADSR](#ADSR)</code>  
<a name="ADSR+processSample"></a>

### adsR.processSample()
Process sample

**Kind**: instance method of <code>[ADSR](#ADSR)</code>  
<a name="ADSR+value"></a>

### adsR.value() ⇒ <code>Number</code>
Get current value

**Kind**: instance method of <code>[ADSR](#ADSR)</code>  
**Returns**: <code>Number</code> - amplitude  
<a name="ADSR+process"></a>

### adsR.process(buffer)
Process a buffer

**Kind**: instance method of <code>[ADSR](#ADSR)</code>  

| Param | Type |
| --- | --- |
| buffer | <code>Array</code> | 

<a name="ADSR+isActive"></a>

### adsR.isActive() ⇒ <code>Boolean</code>
Test if the envelope is active

**Kind**: instance method of <code>[ADSR](#ADSR)</code>  
<a name="ADSR+disable"></a>

### adsR.disable()
Disable the envelope

**Kind**: instance method of <code>[ADSR](#ADSR)</code>  
<a name="Biquad"></a>

## Biquad
**Kind**: global class  
<a name="new_Biquad_new"></a>

### new Biquad(type, sampleRate)
Biquad filter

Created by Ricard Marxer <email@ricardmarxer.com> on 2010-05-23.
 Copyright 2010 Ricard Marxer. All rights reserved.

Implementation based on:
http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt


| Param | Type |
| --- | --- |
| type | <code>Number</code> | 
| sampleRate | <code>Number</code> | 

<a name="DFT"></a>

## DFT
**Kind**: global class  

* [DFT](#DFT)
    * [new DFT(bufferSize, sampleRate)](#new_DFT_new)
    * [.forward(buffer)](#DFT+forward) ⇒

<a name="new_DFT_new"></a>

### new DFT(bufferSize, sampleRate)
DFT is a class for calculating the Discrete Fourier Transform of a signal.


| Param | Type | Description |
| --- | --- | --- |
| bufferSize | <code>Number</code> | The size of the sample buffer to be computed |
| sampleRate | <code>Number</code> | The sampleRate of the buffer (eg. 44100) |

<a name="DFT+forward"></a>

### dfT.forward(buffer) ⇒
Performs a forward transform on the sample buffer.
Converts a time domain signal to frequency domain spectra.

**Kind**: instance method of <code>[DFT](#DFT)</code>  
**Returns**: The frequency spectrum array  

| Param | Type | Description |
| --- | --- | --- |
| buffer | <code>Array</code> | The sample buffer |

<a name="FFT"></a>

## FFT
**Kind**: global class  

* [FFT](#FFT)
    * [new FFT(bufferSize, sampleRate)](#new_FFT_new)
    * [.forward(buffer)](#FFT+forward) ⇒
    * [.inverse(real, imag)](#FFT+inverse) ⇒

<a name="new_FFT_new"></a>

### new FFT(bufferSize, sampleRate)
FFT is a class for calculating the Discrete Fourier Transform of a signal
with the Fast Fourier Transform algorithm.


| Param | Type | Description |
| --- | --- | --- |
| bufferSize | <code>Number</code> | The size of the sample buffer to be computed. Must be power of 2 |
| sampleRate | <code>Number</code> | The sampleRate of the buffer (eg. 44100) |

<a name="FFT+forward"></a>

### ffT.forward(buffer) ⇒
Performs a forward transform on the sample buffer.
Converts a time domain signal to frequency domain spectra.

**Kind**: instance method of <code>[FFT](#FFT)</code>  
**Returns**: The frequency spectrum array  

| Param | Type | Description |
| --- | --- | --- |
| buffer | <code>Array</code> | The sample buffer. Buffer Length must be power of 2 |

<a name="FFT+inverse"></a>

### ffT.inverse(real, imag) ⇒
Performs a inverse FFT transformation
Converts a frequency domain spectra to a time domain signal

**Kind**: instance method of <code>[FFT](#FFT)</code>  
**Returns**: The time domain signal  

| Param | Type |
| --- | --- |
| real | <code>Array</code> | 
| imag | <code>Array</code> | 

<a name="GraphicalEq"></a>

## GraphicalEq
**Kind**: global class  
<a name="new_GraphicalEq_new"></a>

### new GraphicalEq(sampleRate)
Create a Graphical Equalizer

 Implementation of a graphic equalizer with a configurable bands-per-octave
 and minimum and maximum frequencies

 Created by Ricard Marxer <email@ricardmarxer.com> on 2010-05-23.
 Copyright 2010 Ricard Marxer. All rights reserved.


| Param | Type |
| --- | --- |
| sampleRate | <code>SampleRate</code> | 

**Example**  
```js
var eq = new GraphicalEq(44100)
```
<a name="IIRFilter"></a>

## IIRFilter
**Kind**: global class  

* [IIRFilter](#IIRFilter)
    * [new IIRFilter()](#new_IIRFilter_new)
    * _instance_
        * [.set(cutoff, resonance)](#IIRFilter+set)
        * [.process(buffer)](#IIRFilter+process)
        * [.addEnvelope(envelope)](#IIRFilter+addEnvelope)
    * _static_
        * [.LP12](#IIRFilter.LP12)
            * [new IIRFilter.LP12()](#new_IIRFilter.LP12_new)
            * [.addEnvelope(envelope)](#IIRFilter.LP12+addEnvelope)

<a name="new_IIRFilter_new"></a>

### new IIRFilter()
IIRFilter

<a name="IIRFilter+set"></a>

### iirFilter.set(cutoff, resonance)
Set filter parameters

**Kind**: instance method of <code>[IIRFilter](#IIRFilter)</code>  

| Param | Type |
| --- | --- |
| cutoff | <code>Number</code> | 
| resonance | <code>Number</code> | 

<a name="IIRFilter+process"></a>

### iirFilter.process(buffer)
Process a buffer

**Kind**: instance method of <code>[IIRFilter](#IIRFilter)</code>  

| Param | Type |
| --- | --- |
| buffer | <code>Array</code> | 

<a name="IIRFilter+addEnvelope"></a>

### iirFilter.addEnvelope(envelope)
Add an envelope to the filter

**Kind**: instance method of <code>[IIRFilter](#IIRFilter)</code>  

| Param | Type |
| --- | --- |
| envelope | <code>[ADSR](#ADSR)</code> | 

<a name="IIRFilter.LP12"></a>

### IIRFilter.LP12
**Kind**: static class of <code>[IIRFilter](#IIRFilter)</code>  

* [.LP12](#IIRFilter.LP12)
    * [new IIRFilter.LP12()](#new_IIRFilter.LP12_new)
    * [.addEnvelope(envelope)](#IIRFilter.LP12+addEnvelope)

<a name="new_IIRFilter.LP12_new"></a>

#### new IIRFilter.LP12()
LP12 filter

<a name="IIRFilter.LP12+addEnvelope"></a>

#### lP12.addEnvelope(envelope)
Add an envelope to the filter

**Kind**: instance method of <code>[LP12](#IIRFilter.LP12)</code>  

| Param | Type |
| --- | --- |
| envelope | <code>[ADSR](#ADSR)</code> | 

<a name="IIRFilter2"></a>

## IIRFilter2
**Kind**: global class  

* [IIRFilter2](#IIRFilter2)
    * [new IIRFilter2(type, cutoff, resonance, sampleRate)](#new_IIRFilter2_new)
    * [.process(buffer)](#IIRFilter2+process)
    * [.addEnvelope(envelope)](#IIRFilter2+addEnvelope)
    * [.set(cutoff, resonance)](#IIRFilter2+set)

<a name="new_IIRFilter2_new"></a>

### new IIRFilter2(type, cutoff, resonance, sampleRate)
IIRFilter2


| Param | Type |
| --- | --- |
| type | <code>Number</code> | 
| cutoff | <code>Number</code> | 
| resonance | <code>Number</code> | 
| sampleRate | <code>Number</code> | 

<a name="IIRFilter2+process"></a>

### iirFilter2.process(buffer)
Process a buffer

**Kind**: instance method of <code>[IIRFilter2](#IIRFilter2)</code>  

| Param | Type |
| --- | --- |
| buffer | <code>Array</code> | 

<a name="IIRFilter2+addEnvelope"></a>

### iirFilter2.addEnvelope(envelope)
Add an envelope to the filter

**Kind**: instance method of <code>[IIRFilter2](#IIRFilter2)</code>  

| Param | Type |
| --- | --- |
| envelope | <code>[ADSR](#ADSR)</code> | 

<a name="IIRFilter2+set"></a>

### iirFilter2.set(cutoff, resonance)
Set filter parameters

**Kind**: instance method of <code>[IIRFilter2](#IIRFilter2)</code>  

| Param | Type |
| --- | --- |
| cutoff | <code>Number</code> | 
| resonance | <code>Number</code> | 

<a name="MultiDelay"></a>

## MultiDelay
**Kind**: global class  

* [MultiDelay](#MultiDelay)
    * [new MultiDelay(maxDelayInSamplesSize, delayInSamples, masterVolume, delayVolume)](#new_MultiDelay_new)
    * [.setDelayInSamples(delayInSamples)](#MultiDelay+setDelayInSamples)
    * [.setMasterVolume(masterVolume)](#MultiDelay+setMasterVolume)
    * [.setDelayVolume(delayVolume)](#MultiDelay+setDelayVolume)
    * [.process(samples)](#MultiDelay+process) ⇒

<a name="new_MultiDelay_new"></a>

### new MultiDelay(maxDelayInSamplesSize, delayInSamples, masterVolume, delayVolume)
MultiDelay effect by Almer Thie (http://code.almeros.com).
Copyright 2010 Almer Thie. All rights reserved.
Example: http://code.almeros.com/code-examples/delay-firefox-audio-api/

This is a delay that feeds it's own delayed signal back into its circular
buffer. Also known as a CombFilter.

Compatible with interleaved stereo (or more channel) buffers and
non-interleaved mono buffers.


| Param | Type | Description |
| --- | --- | --- |
| maxDelayInSamplesSize | <code>Number</code> | Maximum possible delay in samples (size of circular buffer) |
| delayInSamples | <code>Number</code> | Initial delay in samples |
| masterVolume | <code>Number</code> | Initial master volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |
| delayVolume | <code>Number</code> | Initial feedback delay volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |

<a name="MultiDelay+setDelayInSamples"></a>

### multiDelay.setDelayInSamples(delayInSamples)
Change the delay time in samples.

**Kind**: instance method of <code>[MultiDelay](#MultiDelay)</code>  

| Param | Type | Description |
| --- | --- | --- |
| delayInSamples | <code>Number</code> | Delay in samples |

<a name="MultiDelay+setMasterVolume"></a>

### multiDelay.setMasterVolume(masterVolume)
Change the master volume.

**Kind**: instance method of <code>[MultiDelay](#MultiDelay)</code>  

| Param | Type | Description |
| --- | --- | --- |
| masterVolume | <code>Number</code> | Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |

<a name="MultiDelay+setDelayVolume"></a>

### multiDelay.setDelayVolume(delayVolume)
Change the delay feedback volume.

**Kind**: instance method of <code>[MultiDelay](#MultiDelay)</code>  

| Param | Type | Description |
| --- | --- | --- |
| delayVolume | <code>Number</code> | Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |

<a name="MultiDelay+process"></a>

### multiDelay.process(samples) ⇒
Process a given interleaved or mono non-interleaved float value Array and adds the delayed audio.

**Kind**: instance method of <code>[MultiDelay](#MultiDelay)</code>  
**Returns**: A new Float64Array interleaved or mono non-interleaved as was fed to this function.  

| Param | Type | Description |
| --- | --- | --- |
| samples | <code>Array</code> | Array containing Float values or a Float64Array |

<a name="Oscillator"></a>

## Oscillator
**Kind**: global class  

* [Oscillator](#Oscillator)
    * [new Oscillator(type, frequency, amplitude, bufferSize, sampleRate)](#new_Oscillator_new)
    * [.setAmp(amplitude)](#Oscillator+setAmp)
    * [.setFreq(frequency)](#Oscillator+setFreq)
    * [.add(oscillator)](#Oscillator+add) ⇒ <code>Array</code>
    * [.addSignal(signal)](#Oscillator+addSignal)
    * [.addEnvelope(envelope)](#Oscillator+addEnvelope)
    * [.applyEnvelope()](#Oscillator+applyEnvelope)
    * [.valueAt(offset)](#Oscillator+valueAt)
    * [.generate()](#Oscillator+generate) ⇒ <code>Array</code>

<a name="new_Oscillator_new"></a>

### new Oscillator(type, frequency, amplitude, bufferSize, sampleRate)
Oscillator class for generating and modifying signals


| Param | Type | Description |
| --- | --- | --- |
| type | <code>Number</code> | A waveform constant (eg. DSP.SINE) |
| frequency | <code>Number</code> | Initial frequency of the signal |
| amplitude | <code>Number</code> | Initial amplitude of the signal |
| bufferSize | <code>Number</code> | Size of the sample buffer to generate |
| sampleRate | <code>Number</code> | The sample rate of the signal |

<a name="Oscillator+setAmp"></a>

### oscillator.setAmp(amplitude)
Set the amplitude of the signal

**Kind**: instance method of <code>[Oscillator](#Oscillator)</code>  

| Param | Type | Description |
| --- | --- | --- |
| amplitude | <code>Number</code> | The amplitude of the signal (between 0 and 1) |

<a name="Oscillator+setFreq"></a>

### oscillator.setFreq(frequency)
Set the frequency of the signal

**Kind**: instance method of <code>[Oscillator](#Oscillator)</code>  

| Param | Type | Description |
| --- | --- | --- |
| frequency | <code>Number</code> | The frequency of the signal |

<a name="Oscillator+add"></a>

### oscillator.add(oscillator) ⇒ <code>Array</code>
Add an oscillator

**Kind**: instance method of <code>[Oscillator](#Oscillator)</code>  
**Returns**: <code>Array</code> - the current oscillator signal  

| Param | Type | Description |
| --- | --- | --- |
| oscillator | <code>[Oscillator](#Oscillator)</code> | The oscillator to be added to |

<a name="Oscillator+addSignal"></a>

### oscillator.addSignal(signal)
Add a signal to the current generated osc signal

**Kind**: instance method of <code>[Oscillator](#Oscillator)</code>  

| Param | Type |
| --- | --- |
| signal | <code>Array</code> | 

<a name="Oscillator+addEnvelope"></a>

### oscillator.addEnvelope(envelope)
Add an envelope to the oscillator

**Kind**: instance method of <code>[Oscillator](#Oscillator)</code>  

| Param | Type |
| --- | --- |
| envelope | <code>[ADSR](#ADSR)</code> | 

<a name="Oscillator+applyEnvelope"></a>

### oscillator.applyEnvelope()
Apply the oscillator envelope to its signal

**Kind**: instance method of <code>[Oscillator](#Oscillator)</code>  
<a name="Oscillator+valueAt"></a>

### oscillator.valueAt(offset)
Get value

**Kind**: instance method of <code>[Oscillator](#Oscillator)</code>  

| Param | Type |
| --- | --- |
| offset | <code>Number</code> | 

<a name="Oscillator+generate"></a>

### oscillator.generate() ⇒ <code>Array</code>
Generate the oscillator signal

**Kind**: instance method of <code>[Oscillator](#Oscillator)</code>  
**Returns**: <code>Array</code> - the signal  
<a name="Reverb"></a>

## Reverb
**Kind**: global class  

* [Reverb](#Reverb)
    * [new Reverb(maxDelayInSamplesSize, delayInSamples, masterVolume, mixVolume, delayVolume, dampFrequency)](#new_Reverb_new)
    * [.setDelayInSamples(delayInSamples)](#Reverb+setDelayInSamples)
    * [.setMasterVolume(masterVolume)](#Reverb+setMasterVolume)
    * [.setMixVolume(mixVolume)](#Reverb+setMixVolume)
    * [.setDelayVolume(delayVolume)](#Reverb+setDelayVolume)
    * [.setDampFrequency(dampFrequency)](#Reverb+setDampFrequency)
    * [.process(samples)](#Reverb+process) ⇒

<a name="new_Reverb_new"></a>

### new Reverb(maxDelayInSamplesSize, delayInSamples, masterVolume, mixVolume, delayVolume, dampFrequency)
Reverb effect by Almer Thie (http://code.almeros.com).
Copyright 2010 Almer Thie. All rights reserved.
Example: http://code.almeros.com/code-examples/reverb-firefox-audio-api/

This reverb consists of 6 SingleDelays, 6 MultiDelays and an IIRFilter2
for each of the two stereo channels.

Compatible with interleaved stereo buffers only!


| Param | Type | Description |
| --- | --- | --- |
| maxDelayInSamplesSize | <code>Number</code> | Maximum possible delay in samples (size of circular buffers) |
| delayInSamples | <code>Number</code> | Initial delay in samples for internal (Single/Multi)delays |
| masterVolume | <code>Number</code> | Initial master volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |
| mixVolume | <code>Number</code> | Initial reverb signal mix volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |
| delayVolume | <code>Number</code> | Initial feedback delay volume for internal (Single/Multi)delays. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |
| dampFrequency | <code>Number</code> | Initial low pass filter frequency. 0 to 44100 (depending on your maximum sampling frequency) |

<a name="Reverb+setDelayInSamples"></a>

### reverb.setDelayInSamples(delayInSamples)
Change the delay time in samples as a base for all delays.

**Kind**: instance method of <code>[Reverb](#Reverb)</code>  

| Param | Type | Description |
| --- | --- | --- |
| delayInSamples | <code>Number</code> | Delay in samples |

<a name="Reverb+setMasterVolume"></a>

### reverb.setMasterVolume(masterVolume)
Change the master volume.

**Kind**: instance method of <code>[Reverb](#Reverb)</code>  

| Param | Type | Description |
| --- | --- | --- |
| masterVolume | <code>Number</code> | Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |

<a name="Reverb+setMixVolume"></a>

### reverb.setMixVolume(mixVolume)
Change the reverb signal mix level.

**Kind**: instance method of <code>[Reverb](#Reverb)</code>  

| Param | Type | Description |
| --- | --- | --- |
| mixVolume | <code>Number</code> | Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |

<a name="Reverb+setDelayVolume"></a>

### reverb.setDelayVolume(delayVolume)
Change all delays feedback volume.

**Kind**: instance method of <code>[Reverb](#Reverb)</code>  

| Param | Type | Description |
| --- | --- | --- |
| delayVolume | <code>Number</code> | Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |

<a name="Reverb+setDampFrequency"></a>

### reverb.setDampFrequency(dampFrequency)
Change the Low Pass filter frequency.

**Kind**: instance method of <code>[Reverb](#Reverb)</code>  

| Param | Type | Description |
| --- | --- | --- |
| dampFrequency | <code>Number</code> | low pass filter frequency. 0 to 44100 (depending on your maximum sampling frequency) |

<a name="Reverb+process"></a>

### reverb.process(samples) ⇒
Process a given interleaved float value Array and copies and adds the reverb signal.

**Kind**: instance method of <code>[Reverb](#Reverb)</code>  
**Returns**: A new Float64Array interleaved buffer.  

| Param | Type | Description |
| --- | --- | --- |
| samples | <code>Array</code> | Array containing Float values or a Float64Array |

<a name="RFFT"></a>

## RFFT
**Kind**: global class  

* [RFFT](#RFFT)
    * [new RFFT(bufferSize, sampleRate)](#new_RFFT_new)
    * [.forward(buffer)](#RFFT+forward) ⇒

<a name="new_RFFT_new"></a>

### new RFFT(bufferSize, sampleRate)
RFFT is a class for calculating the Discrete Fourier Transform of a signal
with the Fast Fourier Transform algorithm.

This method currently only contains a forward transform but is highly optimized.


| Param | Type | Description |
| --- | --- | --- |
| bufferSize | <code>Number</code> | The size of the sample buffer to be computed. Must be power of 2 |
| sampleRate | <code>Number</code> | The sampleRate of the buffer (eg. 44100) |

<a name="RFFT+forward"></a>

### rffT.forward(buffer) ⇒
Performs a forward transform on the sample buffer.
Converts a time domain signal to frequency domain spectra.

**Kind**: instance method of <code>[RFFT](#RFFT)</code>  
**Returns**: The frequency spectrum array  

| Param | Type | Description |
| --- | --- | --- |
| buffer | <code>Array</code> | The sample buffer. Buffer Length must be power of 2 |

<a name="Sampler"></a>

## Sampler
**Kind**: global class  
<a name="new_Sampler_new"></a>

### new Sampler(file, bufferSize, sampleRate, playStart, playEnd, loopStart, loopEnd, loopMode)
Sampler


| Param | Type |
| --- | --- |
| file |  | 
| bufferSize | <code>Number</code> | 
| sampleRate | <code>Number</code> | 
| playStart | <code>Number</code> | 
| playEnd | <code>Number</code> | 
| loopStart | <code>Number</code> | 
| loopEnd | <code>Number</code> | 
| loopMode | <code>Number</code> | 

<a name="SingleDelay"></a>

## SingleDelay
**Kind**: global class  

* [SingleDelay](#SingleDelay)
    * [new SingleDelay(maxDelayInSamplesSize, delayInSamples, delayVolume)](#new_SingleDelay_new)
    * [.setDelayInSamples(delayInSamples)](#SingleDelay+setDelayInSamples)
    * [.setDelayVolume(delayVolume)](#SingleDelay+setDelayVolume)
    * [.process(samples)](#SingleDelay+process) ⇒

<a name="new_SingleDelay_new"></a>

### new SingleDelay(maxDelayInSamplesSize, delayInSamples, delayVolume)
SingleDelay effect by Almer Thie (http://code.almeros.com).
Copyright 2010 Almer Thie. All rights reserved.
Example: See usage in Reverb class

This is a delay that does NOT feeds it's own delayed signal back into its
circular buffer, neither does it return the original signal. Also known as
an AllPassFilter(?).

Compatible with interleaved stereo (or more channel) buffers and
non-interleaved mono buffers.


| Param | Type | Description |
| --- | --- | --- |
| maxDelayInSamplesSize | <code>Number</code> | Maximum possible delay in samples (size of circular buffer) |
| delayInSamples | <code>Number</code> | Initial delay in samples |
| delayVolume | <code>Number</code> | Initial feedback delay volume. Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |

<a name="SingleDelay+setDelayInSamples"></a>

### singleDelay.setDelayInSamples(delayInSamples)
Change the delay time in samples.

**Kind**: instance method of <code>[SingleDelay](#SingleDelay)</code>  

| Param | Type | Description |
| --- | --- | --- |
| delayInSamples | <code>Number</code> | Delay in samples |

<a name="SingleDelay+setDelayVolume"></a>

### singleDelay.setDelayVolume(delayVolume)
Change the return signal volume.

**Kind**: instance method of <code>[SingleDelay](#SingleDelay)</code>  

| Param | Type | Description |
| --- | --- | --- |
| delayVolume | <code>Number</code> | Float value: 0.0 (silence), 1.0 (normal), >1.0 (amplify) |

<a name="SingleDelay+process"></a>

### singleDelay.process(samples) ⇒
Process a given interleaved or mono non-interleaved float value Array and
returns the delayed audio.

**Kind**: instance method of <code>[SingleDelay](#SingleDelay)</code>  
**Returns**: A new Float64Array interleaved or mono non-interleaved as was fed to this function.  

| Param | Type | Description |
| --- | --- | --- |
| samples | <code>Array</code> | Array containing Float values or a Float64Array |

<a name="WindowFunction"></a>

## WindowFunction
**Kind**: global class  

* [WindowFunction](#WindowFunction)
    * [new WindowFunction(type, alpha)](#new_WindowFunction_new)
    * [.process(buffer)](#WindowFunction+process)

<a name="new_WindowFunction_new"></a>

### new WindowFunction(type, alpha)
WindowFunction


| Param | Type |
| --- | --- |
| type | <code>Number</code> | 
| alpha | <code>Number</code> | 

<a name="WindowFunction+process"></a>

### windowFunction.process(buffer)
Process a buffer

**Kind**: instance method of <code>[WindowFunction](#WindowFunction)</code>  

| Param | Type |
| --- | --- |
| buffer | <code>Array</code> | 

<a name="DSP"></a>

## DSP
DSP is an object which contains general purpose utility functions and constants

**Kind**: global variable  

* [DSP](#DSP)
    * [.deinterleave](#DSP.deinterleave) ⇒
    * [.getChannel](#DSP.getChannel) ⇒
    * [.invert(buffer)](#DSP.invert) ⇒
    * [.interleave(left, right)](#DSP.interleave) ⇒
    * [.mixSampleBuffers(sampleBuffer1, sampleBuffer2, negate, volumeCorrection)](#DSP.mixSampleBuffers) ⇒
    * [.RMS(buffer)](#DSP.RMS)
    * [.Peak(buffer)](#DSP.Peak)
    * [.mag2db(@buffer)](#DSP.mag2db) ⇒
    * [.freqz(b, a, w)](#DSP.freqz) ⇒

<a name="DSP.deinterleave"></a>

### DSP.deinterleave ⇒
Converts a stereo-interleaved sample buffer into split-stereo (dual mono) sample buffers

**Kind**: static property of <code>[DSP](#DSP)</code>  
**Returns**: an Array containing left and right channels  

| Param | Type | Description |
| --- | --- | --- |
| buffer | <code>Array</code> | A stereo-interleaved sample buffer |

<a name="DSP.getChannel"></a>

### DSP.getChannel ⇒
Separates a channel from a stereo-interleaved sample buffer

**Kind**: static property of <code>[DSP](#DSP)</code>  
**Returns**: an Array containing a signal mono sample buffer  

| Param | Type | Description |
| --- | --- | --- |
| buffer | <code>Array</code> | A stereo-interleaved sample buffer |
| channel | <code>Number</code> | A channel constant (LEFT, RIGHT, MIX) |

<a name="DSP.invert"></a>

### DSP.invert(buffer) ⇒
Inverts the phase of a signal

**Kind**: static method of <code>[DSP](#DSP)</code>  
**Returns**: The inverted sample buffer  

| Param | Type | Description |
| --- | --- | --- |
| buffer | <code>Array</code> | A sample buffer |

<a name="DSP.interleave"></a>

### DSP.interleave(left, right) ⇒
Converts split-stereo (dual mono) sample buffers into a stereo interleaved sample buffer

**Kind**: static method of <code>[DSP](#DSP)</code>  
**Returns**: The stereo interleaved buffer  

| Param | Type | Description |
| --- | --- | --- |
| left | <code>Array</code> | A sample buffer |
| right | <code>Array</code> | A sample buffer |

<a name="DSP.mixSampleBuffers"></a>

### DSP.mixSampleBuffers(sampleBuffer1, sampleBuffer2, negate, volumeCorrection) ⇒
Helper method (for Reverb) to mix two (interleaved) samplebuffers. It's possible
to negate the second buffer while mixing and to perform a volume correction
on the final signal.

**Kind**: static method of <code>[DSP](#DSP)</code>  
**Returns**: A new Float64Array interleaved buffer.  

| Param | Type | Description |
| --- | --- | --- |
| sampleBuffer1 | <code>Array</code> | Array containing Float values or a Float64Array |
| sampleBuffer2 | <code>Array</code> | Array containing Float values or a Float64Array |
| negate | <code>Boolean</code> | When true inverts/flips the audio signal |
| volumeCorrection | <code>Number</code> | When you add multiple sample buffers, use this to tame your signal ;) |

<a name="DSP.RMS"></a>

### DSP.RMS(buffer)
Find RMS of signal

**Kind**: static method of <code>[DSP](#DSP)</code>  

| Param | Type |
| --- | --- |
| buffer | <code>Array</code> | 

<a name="DSP.Peak"></a>

### DSP.Peak(buffer)
Find Peak of signal

**Kind**: static method of <code>[DSP](#DSP)</code>  

| Param | Type |
| --- | --- |
| buffer | <code>Array</code> | 

<a name="DSP.mag2db"></a>

### DSP.mag2db(@buffer) ⇒
Magnitude to decibels

Created by Ricard Marxer <email@ricardmarxer.com> on 2010-05-23.
Copyright 2010 Ricard Marxer. All rights reserved.

**Kind**: static method of <code>[DSP](#DSP)</code>  
**Returns**: the array in decibels  

| Param | Type | Description |
| --- | --- | --- |
| @buffer | <code>Array</code> | The array of magnitudes to convert to decibels |

<a name="DSP.freqz"></a>

### DSP.freqz(b, a, w) ⇒
Frequency response

 Created by Ricard Marxer <email@ricardmarxer.com> on 2010-05-23.
 Copyright 2010 Ricard Marxer. All rights reserved.

 Calculates the frequency response at the given points.

**Kind**: static method of <code>[DSP](#DSP)</code>  
**Returns**: the frequency response in magnitude  

| Param | Type | Description |
| --- | --- | --- |
| b | <code>Number</code> | The coefficients of the filter |
| a | <code>Number</code> | The coefficients of the filter |
| w | <code>Number</code> | The points (normally between -PI and PI) where to calculate the frequency response |

<a name="sinh"></a>

## sinh(num)
Returns the hyperbolic sine of the number

**Kind**: global function  
**Meta**: version: 1004.2314  
**Meta**: discuss at: http://phpjs.org/functions/sinh  
**Meta**: original by: Onno Marsman  

| Param | Type |
| --- | --- |
| num | <code>Number</code> | 

**Example**  
```js
sinh(-0.9834330348825909); // => -1.1497971402636502
```
