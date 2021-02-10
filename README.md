# RT_TDDFT


Here we are going to analyze data extracted from nwchem output in a **real-time tddft** tasl

If we add `print dipole ` in a rt-tddft it will print dipole vector in the given time-steps. Now, if we apply Fourier transformation over these dipole moments to find frequency-dependent linear polarizability $\alpha(\omega)$ and absorption spectra.


## Python script

To do this, you can use this python script <a id="raw-url" href="https://github.com/Yavar-Azar/RT_TDDFT/blob/main/rt-tddft.py">Download FILE</a> to generate spectra image
```bash
python rt-tddft.py output.nwo -wi 0.0 wf 25.0
```
where `wi` and `wf` are used to limit energy window into the desired interval.  
