# xraySimulation
Eric Christie
Based on simulation created by Taylor Barton (maskSim)

Models photons from a laser reflecting off of the Spatial Light Modulator
and reconstruction.

maskSim files are original simulations from Taylor Barton. Some calculations were incorrect and hard to follow so I corrected them with the slmSim_Eric files and intensityCorrelationSim.

Various images are included for testing and verification. The larger the image, the longer it takes to process each pixel.

Files currently in use:
- slmSim_Eric_v1.py: basic reconstruction based on theoretical R(k) value and only low detector noise (+/-2%)
- slmSim_Eric_v1.5_const_phase_shift.py: Adds constant phase shift over whole incoming beam, stepping stone to v2
- slmSim_Eric_v2.py: adds in Static Phase Disorder to reconstruction and uses correction mask (from paper) to correct it on SLM and after values are calculated (When plotting, dots appear throughout image. If the field of view is zoomed into those regions, those dots disappear so I believe it's an error with plotting, not the calculation)
- intensityCorrelationSim_v1.py: Based on generated data, it processed raw frames to generate real R(k) values and reconstruct the image. This pulls mostly from Supplemental Section 1. Needs way to generate data.
