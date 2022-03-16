# Digital Sampling DIY Edition

## Objective

The primary objective of these experiments is to familiarize you with digital data acquisition of time-varying signals. This lab covers concepts in frequency analysis of time varying signals and sampling theory. It also provides an introduction to digital data acquisition (DAQ) systems. You wil use a personal DAQ to sample signals produced by waveforms stored in mp3 formats and converted to analog electrical signals. You will explore issues in: sampling, including the Nyquist limit and aliasing; filtering and its use for noise reduction; and digitization errors. You will also use your DAQ and a microphone to record and analyze sounds in an experiment chosen by you and a partner.

## Background

Most experimental measurements involve the dimension of **time**. Experimental data is acquired over the course of some time, and during this time the actual physical parameter of interest (the _**measurand**_) may change. This could be due to a transient, such as the stress induced in a material by a sudden impact, a periodic phenomena, like the bending and twist of a helicopter blade due to flutter, or random or chaotic fluctuations, like the turbulent velocity in a wind tunnel. Even when the measurand is nominally constant in time, other parameters that influence the measurement may vary, for example drifts in the measurement device. Thus, the experimenter is often interested in measuring a variable that can be described by the general function (or _**waveform**_),

$$
\large v=v(t)
$$

(1)

### Waveforms, Frequency Content and Discrete Sampling

#### Fourier Series

One of the simplest time-dependent functions we encounter is the sine (or cosine\*),

{% hint style="info" %}
\*Either function is acceptable, since $$\sin(wt)=\cos(wt-\pi/2)$$, i.e., the two functions are identical except for a phase difference of π/2 or 90°, meaning that shifted by one-fourth of a cycle, cosine looks just like sine.
{% endhint %}

$$
\large v(t)=A\sin(\omega t+\phi)=A\sin(2\pi ft+\phi)
$$

(2)

where $$A$$ is the amplitude, $$\omega$$ is the circular frequency (e.g., rad/s), $$f$$is the cyclic frequency (e.g., cycles/s, Hertz or $$s^{-1}$$), and $$\phi$$ is the phase, which represents the time-shift of the sine-wave from some reference time that defines _t_=0. Such a function is often denoted as a _simple harmonic_ waveform.

More general periodic waveforms, which repeat themselves with a period _T_ and thus have a frequency $$f$$= 1_/T_, can be **written as a linear combination of simple harmonic** _**modes**_. There is the basic, _**fundamental**_ mode (with frequency$$f$$), and _**harmonics**_ of the fundamental mode, with integer multiples of its frequency (2$$f$$, 3$$f$$, …). For example, we could describe the vibrations of a tuning fork or the acoustic oscillations in a pipe this way. Mathematically, this linear combination of modes is expressed as a **Fourier series expansion**,

$$
\large v(t)=a_0+\sum_{n=1}^\infty\left[a_n\cos(2\pi nft)+b_n\sin(2\pi nft)\right]
$$

(3)

where $$nf$$ represents the frequency of the $$n^{th}$$ mode ($$n$$=1 for the fundamental, $$n$$=2 for the first harmonic, etc.), $$a_0$$ represents the steady component of the waveform, and the $$a_n$$, $$b_n$$ are the harmonic coefficients (or amplitudes) of each mode. The steady amplitude, $$a_0$$, is often called the DC component of the waveform, in reference to classical electrical power systems, which are either **D**irect **C**urrent (steady) or **A**lternating **C**urrent (sinusoidal with a zero average).

![](https://lh4.googleusercontent.com/o3qZr2lYi6rcjH7sJeKmBTw\_2UWN9lMQPjXd0QjSEXS3SLgvAuAg4XScFAJLa9WeW\_SLi4NW6dDWGzjbesBhe1usVjkIKje9VN2jrXL4bmrJNbKhFGaSHHACqVBDqn8jDfR5yB7q)

> **Figure 1.** A waveform composed of a fundamental mode (at 50 Hz) and its 9th harmonic (at 10 times the fundamental frequency, or 500 Hz). The waveform also has a DC, or time-averaged, component of 4 mV. Specifically, the signal (in millivolts) is 4+sin(100πt)+2sin(1000πt), or equivalently, based on cosines, 4+cos(100πt-π/2)+2cos(1000πt-π/2), which simply represents a phase shift of -π/2.

For example, Figure 1 shows a simple waveform composed of two frequencies, a fundamental mode at 50 Hz and its 9th harmonic. Thus the complete waveform is repeated every 20 ms (period=1/fundamental frequency =1/50 s). The waveform shown in the figure also has a DC component. In other words, the signal has a nonzero value when averaged over its period. In general, we can write the DC amplitude as

$$
\large a_0=\frac{1}{T}\int_{-T/2}^{T/2}v(t)dt=f\int_{-T/2}^{T/2}v(t)dt
$$

(4)

The other coefficients of the Fourier expansion are given by

$$
\large a_n=2f\int_{-T/2}^{T/2}v(t)\cos(2\pi nft)dt \\ \space \\ \space \\ b_n=2f\int_{-T/2}^{T/2}v(t)\sin(2\pi nft)dt
$$

(5)

and they can be combined into a complex number (since, $$e^{-ix}=\cos x-i\sin x$$),

$$
\large a_n-ib_n=2f\int_{-T/2}^{T/2}v(t)e^{-i2\pi nft}dt
$$

(6)

The _**power**_, $$P$$, contained in single mode is given by the square of the amplitude

$$
\large P(n)=a_{n}^{2}+b_{n}^{2}
$$

(7)

and the _**phase**_ $$\phi$$ (or phase angle) of a mode is given by

$$
\large \phi (n)=\tan^{-1}(b_n/a_n)
$$

(8)

A second example that shows the ability of a combination of sine waves to create an arbitrary periodic function is shown in Fig. 2. Five sine waves and a DC component (see Fig. 3) were combined to create a function approaching a square wave. While the constructed function resembles a square wave, it is clear that more sine waves would be needed to produce a sharp square wave.

![](<.gitbook/assets/fig 2.png>)

> **Figure 2.** Partial reconstruction of a square wave using five sine waves, each with a different amplitude, frequency and phase, and a separate DC component. The individual waves are shown in Fig. 3.

![](.gitbook/assets/fig3.png)

> **Figure 3.** The five sine waves and constant function used to construct the square wave shown in Fig. 2.

#### Fourier Transforms

The procedure outlined above for periodic functions can be extended to general functions, which are not necessarily periodic, by considering any arbitrary function to be periodic with an infinitely long period. This approach leads to the **Fourier Transform**. Given a function $$v(t)$$, its Fourier Transform $$V(f)$$ is a complex function defined by

$$
\large V(f)=\int_{-\infty}^{\infty}v(t)e^{-i2\pi ft}dt
$$

(9)

in parallel to the complex Fourier function of equation (6). The function $$V(f)$$ represents the information given by $$v(t)$$ _**transformed from the time domain to the frequency domain**_. The transformation is nearly identical in the reverse direction, with simply a change in the phase (note the sign of the exponent), i.e.,

$$
\large v(t)=\int_{-\infty}^{\infty}V(f)e^{+i2\pi ft}df
$$

(10)

For example, Figure 4 graphically shows the Fourier transforms of various functions, including sine and cosine waves, a rectangle function (Π), a triangle function (Λ) and a constant, or DC, function. The sine, cosine and DC waveforms result in Fourier transforms that are nonzero at a single frequency\*\*;

{% hint style="info" %}
\*\*The negative frequencies relate to phase information for the sine and cosine and do not actually represent different frequencies, i.e., for real functions $$v(t)$$ , it can be shown that $$|V(f)|=|V(-f)|$$. That means that if you take the absolute value of V, the part of V below 0 frequency looks like a reflection of the part for f>0.
{% endhint %}

in other words, they contain information at only one frequency (the DC function, which does not change in time, is associated with a frequency of zero). The Fourier transforms of the rectangle and triangle functions result in $$sinc$$ and $$sinc^2$$ functions, where $$sinc(f)=\sin(\pi f)/\pi f$$ , which contains information at many frequencies, but with multiple frequency “peaks”.

![](.gitbook/assets/fig4.png)

> **Figure 4.** Fourier transforms of various functions (left and right pairs). The arrows represent impulse functions (i.e., delta functions), which extend infinitesimally along the x-axis, but have a integrated area corresponding to the height indicated by the arrow. The dashed regions indicate imaginary values.

Instead of looking at the Fourier transform, we often are interested in the _**power spectrum**_ (or _**power spectral density,**_ **PSD**) of a waveform. This represents the amount of power or energy in a region between $$f$$ and $$f+df$$. For real (noncomplex) functions $$v(t)$$, this is given by

$$
\large PSD(f)=|V(f)|^2
$$

(11)

where it is sufficient to consider only 0<$$f$$<∞ since the PSD of a real function is symmetric about $$f$$=0.\*\*

Thus the PSD of the rectangle function, $$\Pi(x)$$ as shown in Figure 4, is the square of its Fourier transform, or $$sin c^2(f)$$ (also shown in Figure 4).

Extensions of the Fourier Transform method have been developed for non-continuous functions, specifically for signals that have been discretely sampled by a computer, data acquisition system, or produced by digital means. These are generally known as Discrete Fourier Transforms. In addition, methods to quickly compute the Fourier Transform have also been developed, e.g., the Fast Fourier Transform. These concepts are described in detail in references 2 and 4. The digital data acquisition system you will use employs these techniques to compute the power and phase spectra of the signals that are sampled in this lab.

#### Discrete Sampling

In most situations, especially for digital data acquisition, the continuous function $$v(t)$$ is sampled (i.e., the data is acquired) at evenly spaced, discrete intervals in time, separated by an amount $$\Delta t$$. The sampling frequency (or data acquisition rate) is thus $$f_s=1/\Delta t$$.

For a given sampling rate, we might ask how accurately the discretely acquired data can reproduce the actual waveform being sampled. The answer depends on the frequency content of the waveform and a special frequency, called the _**Nyquist frequency**_ $$(f_N)$$, which is half the sampling frequency, i.e., $$f_N=f_s/2$$. If the waveform contains no components above the Nyquist frequency, then the waveform can be completely determined by the sampled data (assuming no errors in the measurement).\*\*\* This is known as the _**Nyquist/Nyquist-Shannon Sampling Theorem**_.

{% hint style="info" %}
\*\*\*A waveform that has information in only a limited range of frequencies is called **bandwidth limited**. Due to phase ambiguity, the sampling frequency should actually be more than twice the maximum frequency in the waveform. For example, a sine wave sampled at 0, $$\pi$$, $$2\pi$$, etc. would always have a 0 result and could be confused with a null function.
{% endhint %}

As a simple example, consider a single sine wave. If we know we are dealing with a single frequency sine wave, it takes at least two measurements per period to determine its frequency, which means we must sample at twice the sine wave’s frequency. If we sample any slower, we actually infer a lower frequency than the actual frequency of the sine wave (you will see this in the lab). This process, by which information at a higher frequency shows up at a lower frequency is known as _**aliasing**_**.**

Aliasing occurs for any sampled waveform having components with frequencies above the sampling system’s Nyquist frequency, i.e., $$(f>f_N)$$. One way to remove this problem is to filter the data before it is sampled. This can be accomplished by a low pass filter, a filter that only passes frequencies below some cut-off frequency. One would set the cut-off at or below the Nyquist frequency. The high frequency information is thus removed before it can be aliased. In essence, the filter produces a bandwidth limited waveform.

### Digital Data Acquisition Systems

Data will be acquired with a standalone digital data acquisition system (**DAQ**) that communicates with your computer through a USB connection and using a LabView™ software interface. Most DAQs can be connected to more than one input source; each signal (e.g., a voltage) is connected to one channel of the DAQ. A typical DAQ consists of a multiplexer, a sample-and-hold device, an amplifier, an analog-to-digital converter, a memory buffer, a microcontroller, and an interface to a computer (see Figure 5).

![](<.gitbook/assets/DAQ schematic 1.png>)

**Figure 5.** Schematic of multiplexed, sequential sampling, digital data acquisition system and its connection to a computer.

The _**multiplexer**_** ** (MUX) is a switch that connects one of a number of input channels (usually numbered starting at 0) to the _**sample-and-hold**_** ** (S/H). The input voltage on the channel switched by the MUX “charges up” the sample-and-hold during some time interval, which is a fraction of the sampling period (the time between samples). This circuit is then disconnected from the input voltage, and some of the stored charge is drained from it. The amount of charge leaving during this time is proportional to the original input voltage. The output of the S/H is amplified and then converted to a digital value by the _**analog-to-digital converter**_** ** (ADC). The digital result is then moved to the buffer memory, and communicated to the computer.&#x20;

The digital value produced by the ADC (sometimes referred to as a “word” of data) depends not only on the input voltage, but also on the _**voltage range**_ and _**number of bits**_ of the ADC/amplifer system. The range is given by the _minimum_ and _maximum_ voltages that the ADC/amplifier can read (e.g., 0 and 5 V). The number of bits  ($$N$$) in the ADC determines its digital dynamic range (= $$2^N-1$$). Thus the relation between the digitizer output and the voltage input is given by

$$
\large output=\frac{input-minimum}{maximum-minimum}\times\left(2^N-1\right)
$$

(12)

where output has to be an _integer value_. As an example, for a 2.05 V input into a DAQ with a 0-10 V range, and an 8-bit digitizer (possible digital values of 0-255), the output value would be 52 (not 52.275). Any signal amplitude variations below the difference between two adjacent quantized levels are lost; this is known as the _**quantization error**_ =$$(maximum-minimum)/2^N$$. In the example above, we can only say the input value was 2.039V$$\pm$$0.0196 V (assuming the example ADC rounds rather than truncates). One would normally choose an ADC with a number of bits sufficiently high that the quantization error is less than the dominant sources of error in the measurement. Other factors, though, may influence the choice of ADC bits, including cost and data storage requirements, both of which increase with the added number of bits.&#x20;

Multiple signal inputs are recorded by using the MUX to cycle through each of the input channels at a rate that must be faster than the overall sampling rate (how often a given channel is read) times the number of input channels being read. In the sequential sampling system illustrated in Fig. 5 (and which is representative of the system you will be using), note that the channels are _not read at exactly the same time_. There is a time delay (**skew**) between when one channel and the next is read. The skew determined by the maximum switching and reading rates of the MUX, S/H and ADC. This is illustrated in Fig. 6. Simultaneous data acquisition systems, which have negligible skew, typically employ multiple, synchronized S/H systems just upstream of the MUX (see Fig. 7).

![](https://lh6.googleusercontent.com/C\_ONWpEQ-MIyPf-a\_pjkBNn0w03MAVeAAQuU9N19xQx7Uhz6DZORXVICw6ZwGgtYg384BdTbTlqj-fKJFZ-oI8uZkfNqopfebYQG3pRhhqEP6UvHu4j7Va1BwOR2YK9yh\_PF9wvRhnYQG6YjRQ)

> **Figure 6.** Time delay (skew) between successive channels in sequential sampling system.
>
>

![](<.gitbook/assets/DAQ schematic 2 (1).png>)

> **Figure 7.** Schematic of simultaneous sampling, digital data acquisition system.

In this lab, you control the data acquisition process through a software interface called a LabView _virtual instrument_ (VI). The VI creates a display on the computer screen that lets you think of the data acquisition system as a box with “knobs”, “dials”, and other displays. For this experiment, the VI allows you to control parameters such as the minimum and maximum voltages read by the DAQ, the sampling rate$$(f_s)$$, and the number of samples recorded.



### Sampling/Digital Data Acquisition Terminology

* **Sample** = a single measurement (i.e., at an "instant" in time) captured by the DAQ from one channel
* **Sampling period** = the time between two successive samples &#x20;
* **Sampling rate** = 1/sampling period&#x20;
* **Record** = a group of successive samples acquired by the DAQ&#x20;
* **Record length** = the number of samples in a record
* **Record duration** = the time between the first and last sample in a record&#x20;



## Procedure

**Week 1** - During your lab time

1. **Pickup DIY kit:**
   * _Labjack DAQ, USB cord, 3.5 mm audio jack cord_
   * _Microphone and cord_
   * _Watch video on components and connections_
2. **Download and software:**
   * Watch video on installing software and testing system
   * Download and install executable Labview VI&#x20;
   * Download mp3 tracks to your computer or mobile device
3. **Test your devices:**
   * Connect DAQ to computer using USB cable
   * Connect device audio output to DAQ using 3.5mm audio jack cable
   * Test your system - if you have trouble, ask TA for help

**Week 2** - At a  time of your choosing (but working individually)

1. **Perform initial waveform identification:**
   * There are 12 audio tracks. Each audio track contains a different periodic signal. These signals include: **single sine waves** (at different frequencies), a **sum of three sine waves** (each at a different frequency)**,** a **product of two sine waves** (e.g., sin(At) sin(Bt), also known as amplitude modulation), **a sine wave of a sine wave** (e.g., sin(sin(At), also known as frequency modulation), and periodic waveforms that are not sine waves: **square waves**, **triangle waves**, and **ramps**. Some tracks also have "noisy" versions of some of these waveforms.
   * Play each track on your device, listen to it, and write down your guess at what that track is (it is okay to guess wrong)&#x20;
2. **Setup DAQ system**
   * Connect DAQ to computer using USB cable
   * Connect audio output of the device you are using to play the tracks to the DAQ using the 3.5mm audio jack cable
   * Open the Labview VI
3. **Use DAQ to perform waveform identification**
   * &#x20;Run VI
   * Use the VI to choose the following settings
     * Sampling rate = 25kHz
     * Record length = 10,000
     * ???
   * Use your device to play one track at a time, in repeat mode
   * Observe the time plot and power spectrum and use them to identify the waveform on each track (you should adjust the output/volume level to make sure most of the tracks have a peak voltage of a few volts).&#x20;
     * For each track, write down the frequencies for each of the peaks you see in the power spectrum&#x20;
4. **Examine complex waveforms and interpret power spectrum**
   1. Play the track you identified as product of sines (**amplitude modulation**)&#x20;
      * From the time plot, determine the period of the wave (time between peaks)
        * _Tip: With the waveform displayed as you like, you can toggle the Continuous/Hold switch to the Hold position so that the display just shows the last data captured (doesn't keep taking new samples)._&#x20;
      * Is the frequency of the waveform based on the period the same as any of the frequencies you wrote down from the power spectrum? If not, why?  &#x20;
        * _Hint: the power spectrum shows the frequencies needed to produce that signal from a **sum** of sines!!!_&#x20;
   2. Play the track you identified as the triangle wave
      * &#x20;__ From the time plot, again determine the period of the wave (time between peaks); how is it related to the frequencies you wrote down from the power spectrum?
   3. Alternate playing the triangle wave track and the square wave track
      * Compare the heights (power) of each peak in the power spectrum; which waveform has more power at high frequencies?&#x20;
      * Can you think of a reason why that waveform contains more high frequency content?
5. **Examine quantization error**
   * Play the track with the square wave, and make sure the VI is set to _Continuous_ and both y-axes are set to _Autoscale_
   * Observe both the time plot and power spectrum as you reduce the output volume on the device playing the tracks; what changes do you observe when the volume is set very low &#x20;
6.  **Gather data to explore the effects of varying record length and sampling rate:**

    * Play the track containing the 1 kHz sine wave
    * At the sampling rate and record length combinations shown in the table below, determine **record time**, **power spectrum frequency range**, **power spectrum frequency resolution**, and **number of points in the power spectrum**. Do this by adjusting the x-axis limits on both the time history and the power spectrum as needed, directly observing and noting down each of the required variables.
    *
    * Frequency resolution = the frequency spacing between two points in the power spectrum
    *

    ****
7.  ****

    **Gather data to understand Nyquist sampling theory and aliasing:**

    * In this step we will determine how **aliasing,** brought on by sampling rates below the Nyquist frequency, affects our ability to accurately reconstruct/analyze a signal.
    * Disconnect the T-connector and connect the 3.5mm-to-BNC cable directly to the DAQ's AI0 port. Turn off the oscilloscope.
    * Find and play the track containing the 1 kHz sine wave.
    * With a sampling rate of 2500 S/s and record length of 2500 S, acquire a power spectrum. Record at what frequency in the spectrum the peak occurs (i.e. the frequency with the maximum power).
    * Repeat the above step for the following seven sampling rates: 2000, 1500, 1200, 1000, 800, 675 and 665 S/s, in each case setting the record length value to the sampling rate value (i.e. capture X samples at X S/s).
    * From the 8 observed frequencies, identify at which sampling rate(s) aliasing is occurring.
    * For at least two additional sampling rates of your choice (below 650 S/s, with a matching record length as before), first predict whether aliasing will occur. If you believe aliasing will occur, predict the specific aliasing frequency, then acquire data to verify experimentally.
8. **Gather data to explore the effects of varying record length and sampling rate:**
   * Some important terms:
     * Sample = a single measurement captured by the DAQ
     * Record = a batch of samples collected by the DAQ before downloading to the VI
     * Record length = the number of samples in a record
     * Record time = the period over which the record was captured
     * Sampling time = the period between successive samples
     * Sampling rate = the number of samples acquired in a given period of time
     * Frequency resolution = the frequency spacing between two points in the power spectrum
   * Again play the track containing the 1 kHz sine wave.
   * At the sampling rate and record length combinations shown in the table below, determine **record time**, **power spectrum frequency range**, **power spectrum frequency resolution**, and **number of points in the power spectrum**. Do this by adjusting the x-axis limits on both the time history and the power spectrum as needed, directly observing and noting down each of the required variables.\
     ![](.gitbook/assets/SamplingPoints.PNG)
9. **Explore the implementation of a low pass filter to remove unwanted noise:**
   * Reconfigure the cables/connectors:
     * Remove the BNC cable from **AI0** and replace the T-connector back into AI0.
     * Connect the 3.5mm-to-BNC cable from the MP3 player to the T-connector.
     * Connect the remaining open port on the T-connector to either input of the Krohn-Hite filter using a coaxial jumper cable.
     * Connect the corresponding output from the Krohn-Hite filter to **AI1** on the DAQ, ensuring the switch is set to FS (floating source), with a coaxial jumper cable.
   * Power on the filter and set it to LOW PASS x100 mode
   * Locate and play the track containing the sum of three sine waves at three frequencies on the MP3 player.
   * Set the following VI parameters:
     * Sampling rate = 22 kS/s
     * Record length = 1000 S
     * Number of averages = 1
     * Display Settings
       * Window = None (Uniform) | Vrms | Linear
       * Plot = Amplitude | Radians
     * In both time plots, turn off x-axis autoscale and set the limits from 0 to 0.005 s (will round up to 0.1 after hitting enter)
   * In this step, the sum of sines represents a fictional scenario whereby a signal with two low frequency components of interest (e.g. vibration data from a structures experiment) are subject to a high frequency noise component. Your goal is to remove the high frequency noise without altering the two low frequency components of the signal of interest. This is a very common scenario for signal processing in engineering and science.
   * Paying attention to the time history and power spectrum of both filtered and unfiltered signals, adjust the cutoff frequency dial of the low pass filter until you obtain a cleaned up filtered signal. When you are happy with your results, take a screenshot of the VI for your report and note down the cutoff frequency on the filter.
10. **Explore the ramification of low pass filtering on signal phase**
    * Continue playing the previous track (sum of three sine waves).
    * Set the Krohn-Hite filter to a cutoff frequency of 10 kHz
    * Set the following VI settings:
      * Sampling rate = 22 kS/s
      * Record length = 1000 S
      * Number of averages = 1
      * Display Settings
        * Window = None (Uniform) | Vrms | dB
        * Plot = Phase | Radians
    * Set the VI to CONTINUOUS and press HOLD to freeze the plots after a few seconds, once new records have downloaded.
    * Successively zoom in to the 3 frequencies of interest in each phase plot by adjusting the x-axis limits.
    * Record the phase of the unfiltered and filtered signals at each frequency. Have the TAs check your data.
11. **Gather data to determine the transfer function of the low pass filter at 800 Hz**
    * Locate and play the repetitive sweeping track on the MP3 player
    * Set the following VI settings:
      * Sampling rate = 6000 S/s
      * Record length = 3000 S
      * Number of averages = 1
      * Display Settings
        * Window = None (Uniform) | Vrms | Linear
        * Plot = Amplitude | Radians
    * Set the Krohn-Hite cutoff frequency dial to 800 Hz.
    * With data acquiring, click HOLD and then TAKE NEXT as many times as needed to display a relatively flat and noise-free power spectrum. The TAs will help you achieve this.
    * When you are happy, click SAVE and choose a useful filename, being sure to add .xls as an extension (this can be added in Windows Explorer afterwards if this step is forgotten).
    * Change Number of averages to 10, repeating the previous 2 steps to acquire a new data set.
12. **Explore the implementation of a band pass filter to remove all frequency content except for one frequency of interest:**
    * Locate and play the excessively noisy single sine wave track on the MP3 player. Ensure that the MP3 player volume is set to 32.
    * Set the following VI settings:
      * Sampling rate = 22000 S/s
      * Record length = 5000 S
      * Number of averages = 1
      * Display Settings
        * Window = None (Uniform) | Vrms | dB
        * Plot = Amplitude | Radians
    * Remove the output of the lowpass filter from DAQ AI1 and instead connect it to input of the currently unused channel of Krohn-Hite filter.
    * Connect the output of this channel to DAQ AI1 using a further coaxial jumper cable.
    * Set the second filter channel to HIGH PASS x1 mode and set the cutoff frequency dial such that the frequency is as low as possible (around 20 Hz).
    * Change the original low pass filter to x100 and set the cutoff frequency dial such that the frequency is much higher than the frequency we are trying to preserve (around 20 kHz).
    * On **both** of the filtered and unfiltered time history plots:
      * Turn off x-axis auto-scale and set the x-axis limits between 0 and 0.01 s.
      * Turn off y-axis auto-scale and set the y-axis limits between -0.1 and 0.1 V.
    * On both power spectrum plots:
      * Turn off x-axis auto-scale and set the x-axis limits between 0 and 11 kHz.
      * Turn off y-axis auto-scale and set the y-axis limits between -130 and -30 dBVrms (-30 on the top-most limit).
    * Take a screenshot of the VI for future reference.
    * Gain a high-level perspective of how both filters affect the power spectrum:
      * Whilst watching the filtered power spectrum and the time history plots, slowly reduce the low-pass filter cutoff frequency to 2000 Hz. Observe how both plots change.
      * Zoom closer into the frequency of interest by adjusting the x-axis limits of both power spectra to between 0 and 2000 Hz.
      * Whilst watching the filtered power spectrum, slowly increase the high-pass filter cutoff frequency to 200 Hz. Observe how the power spectrum and the time history plots change.
      * On each power spectrum, drag marker "m1" until it snaps onto the peak of the frequency of interest. The y-axis readout (in orange below the x-axis) for each marker tells you the dBVrms of each marker, filtered and unfiltered, and thus the magnitude of the signal at that frequency. If the above steps have been followed correctly, you should see almost identical dB values at these points; have the TA check your setup if this is not the case.
      * Note down the dBVrms of each peak.
    * Dial in the cut-off frequencies to complete your band-pass filter design:
      * Keeping the low pass cutoff frequency at 2000 Hz, adjust LOW PASS mode to x10.
      * Keeping the high pass cutoff frequency at 200 Hz, adjust HIGH PASS mode to x10.
      * Successively adjusting both dials until the filtered power spectrum has as much unwanted frequency content removed as possible, without exceeding a 3 dBVrms drop at the frequency of interest. The final filter cutoff frequencies should be evenly spaced around the frequency of interest (i.e. 1000 +/- X Hz).
      * Note down your final band-pass filter frequencies, the final dBVrms values of each peak (unfiltered and filtered).
      * Take a screenshot of the LabView VI at its current zoom level.
      * Set both power spectrum x-axis limits between 0 and 11 kHz before taking another screenshot of the LabView VI.
13. **Lab shutdown procedure**
    1. Plug the MP3 player in to charge
    2. Unplug all coaxial cables and arrange neatly on the desk
    3. Turn off filter, oscilloscope, and DAQ
    4. Have the TA upload your data/files to Canvas

## Data to be Taken

1. The waveform shape/description and peak frequency(s) for each of the 11 tracks.
2. For the 1 kHz waveform:
   * [ ] The peak frequency from the power spectrum for the 10 different sampling rate/record length combinations.
   * [ ] The predicted and measured peak frequency for the 2 selected sampling rate/record length combinations.
   * [ ] A note saying which, if any, of the frequencies are aliased frequencies.
3. For the 1 kHz waveform at the 9 different sampling rate/record length combinations:
   * [ ] Record time.
   * [ ] Power spectrum x-axis range.
   * [ ] Power spectrum frequency resolution.
   * [ ] The number of points in the power spectrum.
4. For the sum of 3 sines waveform:
   * [ ] The cutoff frequency of the low-pass filter design and a screenshot of the LabView VI with the lowpass filter set at the cutoff frequency.
5. For the sum of 3 sines waveform run through a lowpass filter set at 10 kHz:
   * [ ] The phase of the filtered and unfiltered signal at all 3 frequencies of interest (6 values total).
6. For the sweeping sine waveform run through a lowpass filter set at 800 Hz:
   * [ ] Two excel documents, one containing unfiltered and filtered power spectrum data from the single record run, and one for the average of 10 records run.
7. For the excessively noise sine waveform:
   * [ ] Designed bandpass filter range.
   * [ ] dBVrms values of both peaks before and after filtering.
   * [ ] All three screenshots taken of the LabView VI.

## **​Data Reduction**

1. From your data taken in step 3 of the Data to be Taken section, come up with rules for:
   1. The number of points in the power spectrum as a function of record length
   2. The power spectrum frequency resolution as a function of record time
   3. The power spectrum frequency resolution as a function of sampling rate and record length
   4. The power spectrum x-axis limit as a function of the sampling rate.
2. From your data taken in step 5 of the Data to be Taken section, calculate the phase difference between the unfiltered and filtered signals at all 3 frequencies of interest.
3. Calculate the ratio of the power spectra for the filtered data (i.e., the output from the filter) and unfiltered data (i.e., the input to the filter) recorded in step 6 of the Data to be Taken section. Find the Transfer Function of the filter as a function of frequency. Do this for both the single record run and the average of 10 records run.

## **​Results Needed for Data Report**

**Note: This is a Data Report, so you must follow instructions on how to prepare the Data Report (not the FORMAL report) on the Canvas course page. Also, be sure to answer any supplemental questions listed on Canvas for this lab.**

1. A table describing of each of the 11 tracks including all relevant characteristics (Step 1 of Data to be Taken).
2. A table containing the data from Step 2 of the Data to be Taken section.
3. A table containing the data from Step 3 of the Data to be Taken section.
4. All 4 equations for Step 1 of the Data Reduction section.
5. A table containing the design cutoff frequency of the lowpass filter (step 4 of Data to be Taken section).
6. A figure of the lowpass filtered and unfiltered time history and power spectra from step 4 of the Data to be Taken section.
7. A table containing the data from step 5 of the Data to be Taken Section and Step 2 of the Data Reduction section.
8. Plots of the single and average power spectra, filtered and unfiltered, for the waveform with the rapidly sweeping frequency (step 6 of the Data to be Taken section)
9. Plots of the filter transfer function described in step 3 of the Data Reduction, one based on 1 record and another based on the 10-record average results.
10. For the bandpass filter design (Step 7 of Data to be Taken):
    1. A table containing:
       * [ ] The final bandpass filter design cutoff frequencies.
       * [ ] Unfiltered/filtered peak dbVrms values prior to finalizing the filter design.
       * [ ] Unfiltered/filtered peak dbVrms values after finalizing the filter design.
    2. Three screenshots of the LabView VI during filter development

## **Further Reading**

{% embed url="https://www.youtube.com/watch?ab_channel=UncleDoug&v=ueOup-XBexU" %}
**How to use an analog oscilloscope**
{% endembed %}

1.  R. V. Churchill and J. W. Brown, _Fourier Series and Boundary Value Problems_, 3rd ed., McGraw-Hill, 1978.

    ​
2.  R. N. Bracewell, _The Fourier Transform and Its Applications_, 2nd ed., McGraw-Hill, 1978.

    ​
3.  T. G. Beckwith, R. D. Marangoni and J. H. Lienhard V, _Mechanical Measurements_, 5th ed., Addison-Wesley, 1995.

    ​
4. W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flanner\_, Numerical Recipes - The Art of Scientific Computing\_, 2nd ed., Cambridge University Press, 1992.

**​**
