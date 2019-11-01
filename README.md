# TrackerMCTune
Code for tuning the simulation and response of the tracker.

## Setup
To build, just run make
To run the included fcl, you will need to source env.sh so that Offline knows where to find geometry files etc.

## Procedure for tuning

### Determining waveform parameters

1. cd waveform, run FitCurrentPulse.py
2. run ExtractCurrentPulse.py to get values for
   - strawElectronics.currentSigmas
   - strawElectronics.currentT0s
   - strawElectronics.reflectionTimeShift
   - strawElectronics.reflectionVelocity
   - strawElectronics.reflectionALength
   - strawElectronics.reflectionFrac
   - strawElectronics.currentMeans
   - strawElectronics.currentNormalizations
3. run FitADCPulse.py to get values for the Fermilab setup for
   - strawElectronics.ADCPoles
4. run FitADCPulse2.py to get values for the LBL setup for
   - strawElectronics.ADCPoles
5. Set these parameters in fcl/epilog.fcl

### Getting gain and noise parameters for Fe55 tune

1. Analyze pedestals to get ADC noise
   - bin/pedestal_data
     * strawElectronics.adcAnalogNoise
2. Separately get gain of ADC circuit for Fermilab (updated) setup and LBL (older) setup (integrator has changed since prototype data was taken). For each setup:
   - bin/fe55adc_data
   - mu2e -c fcl/Fe55adc_gain_<lbl/fermi>.fcl
   - bin/fe55adc_sim
   - compare sim and data with bin/fe55adc_compare
   - in fcl/epilog_Fe55_<lbl/fermi>.fcl, adjust until you get agreement values of
     * strawElectronics.defaultAdcdVdI
     * strawElectronics.adcAnalogNoise/strawNoise
3. The preamps should be the same between both setups, so the threshold gain is the same
   - bin/fe55threshold_data
   - mu2e -c fcl/Fe55_thresholdgain_<12/16/20/24>.fcl
   - bin/fe55threshold_sim
   - compare with bin/fe55threshold_compare, in fcl/epilog_Fe55.fcl adjust values of
     * strawElectronics.defaultThresholddVdI
     * strawElectronics.thresholdAnalogNoise/strawNoise
4. Saturation effects also are the same between setups. Here we have data from Fermilab setup
   - bin/fe55saturation_data
   - mu2e -c fcl/fermi_Fe55_saturation_<1/2/4/6/8/10/12/14>.fcl
   - bin/fe55saturation_sim
   - bin/fe55saturation_compare, in fcl/epilog.fcl adjust value of
     * strawElectronics.saturationVoltage
