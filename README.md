# generalizedchip

### A generalized precession parameter to interpret gravitational-wave data

> Originally designed for waveform approximants, the effective precession parameter `chip` is the most commonly used quantity to characterize spin-precession effects in gravitational-wave observations of black-hole binary coalescences. We point out that the current definition of `chip` retains some, but not all, variations taking place on the precession timescale. We rectify this inconsistency and propose more general definitions that either fully consider or fully average those oscillations. Our generalized parameter 0<=`chip`<=2 presents an exclusive region `chip`>1 that can only be populated by binaries with two precessing spins. We apply our prescriptions to current LIGO/Virgo events and find that posterior distributions of `chip` tend to show longer tails at larger values. This appears to be a  generic feature, implying that (i) current `chip` measurement errors might be underestimated, but also that (ii) evidence for spin precession in current data might be stronger than previously inferred. Among the gravitational-wave events released to date, that which shows the most striking behavior is GW190521. 


This repository contains a short code supporting [arXiv:2011.XXXXX](https://arxiv.org/abs/2011.XXXXX). We are very happy if you find this useful for your research; please cite our paper. For a DOI pointing to this repository: [![DOI](https://zenodo.org/badge/XXXXXXX.svg)](https://zenodo.org/badge/latestdoi/XXXXXXX)

This code is developed and maintained by [Davide Gerosa](https://davidegerosa.com/). To report bugs, please open an issue on GitHub. If you want to contact me, it's `d.gerosa@bham.ac.uk`. 


### Usage

The main usage is as follows:

      q=0.7
      chi1=0.3
      chi2=1
      theta1=np.pi/3
      theta2=np.pi/4
      deltaphi=np.pi/5
      fref= 20 # Hz
      M_msun=60 # Msun
      print("Heuristic chip:", chip_heuristic(q,chi1,chi2,theta1,theta2))
      print("Asymptotic chip:", chip_asymptotic(q,chi1,chi2,theta1,theta2))
      print("Generalized chip:", chip_generalized(q,chi1,chi2,theta1,theta2,deltaphi))
      print("Averaged chip:", chip_averaged(q,chi1,chi2,theta1,theta2,deltaphi,fref=fref,M_msun=M_msun))

Alternatively, one can specify `r` instead of `fref`:

      r=10
      print("Averaged chip:", chip_averaged(q,chi1,chi2,theta1,theta2,deltaphi,r=r))

The code is numpy-vectorized:

      q=np.array([0.7,0.7])
      chi1=np.array([0.3,0.3])
      chi2=np.array([1,1])
      theta1=np.array([np.pi/3,np.pi/3])
      theta2=np.array([np.pi/4,np.pi/4])
      deltaphi=np.array([np.pi/5,np.pi/5])
      fref= np.array([20,20]) # Hz
      M_msun= np.array([60,60]) # Msun
      print("Heuristic chip:", chip_heuristic(q,chi1,chi2,theta1,theta2))
      print("Asymptotic chip:", chip_asymptotic(q,chi1,chi2,theta1,theta2))
      print("Generalized chip:", chip_generalized(q,chi1,chi2,theta1,theta2,deltaphi))
      print("Averaged chip:", chip_averaged(q,chi1,chi2,theta1,theta2,deltaphi,fref=fref,M_msun=M_msun))

