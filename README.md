# CH4-CO-OH BoxModel
README file for the 2-box model from Nguyen et al.

Developed by Newton Nguyen, Alex Turner, and Christian Frankenberg

February 14, 2020

# Info
This code models the interactive methane, carbon monoxide, and
hydroxyl chemistry in Earth's atmosphere, which was used in Nguyen et
al, 2020. There are two tropospheric boxes and to optional stratospheric
boxes. The simplified tropospheric chemistry is from Prather 1994. 

# Info
- We have included all the tests seen in Nguyen et al, 2020 (Refer to below
  table)
- The MATLAB code can be run using the "DriverScript.m" file.
- The inversions use publicly available datasets.  As such, we have not included those datasets in our tarball.
- We have included scripts to download all of the datasets though.  This can be done by running the "download_data.csh" script.
- Users could also manually download the data by navigating the websites (see below), contacting the PIs for those datasets, contacting Alex Turner (aturner@fas.harvard.edu), or Christian Frankenberg (cfranken@caltech.edu).
- The code will still run without the original datasets and will just use the hemispheric averages that were bootstrapped from Turner et al.
- The code in the tarball is consistent with the code on the github (as of December 23, 2016).


# Public datasets:
- CH4 from NOAA/ESRL
 - PI: Ed Dlugokencky
 - URL: "ftp://aftp.cmdl.noaa.gov/data/trace_gases/ch4/flask/"
- MCF from NOAA/ESRL
 - PI: Steve Montzka
 - URL: "ftp://aftp.cmdl.noaa.gov/data/hats/solvents/CH3CCl3/"
- MCF from GAGE/AGAGE
 - PI: Ron Prinn
 - URL: "http://agage.mit.edu/"
- CH4C13 from NOAA/ESRL
 - PI: James White
 - URL: "ftp://aftp.cmdl.noaa.gov/data/trace_gases/ch4c13/flask/"
- CH4C13 from U. Heidelberg
 - PI: Ingeborg Levin
 - URL: "http://www.iup.uni-heidelberg.de/institut/forschung/groups/kk/Data_html"
- CH4C13 from U. Washington
 - Data is included in the tarball
- CH4C13 (MCF) from UC Irvine
 - Data is included in the tarball
- C2H6 data (wrote an observation operator but did not use in the manuscript)
 - URL: "http://ds.data.jma.go.jp/gmd/wdcgg/cgi-bin/wdcgg/download.cgi?para=VOCs"


# Table 1: Tests Scripts

|  Test Description | Corresponding figure in Nguyen et al | File Name | 
| --- | --- | --- |
| Test instantaneous 10-Tg-perturbation lifetime | Fig. 1A | ForwardModelTest.m
| Forward model runs to test equilibrium concentrations | Fig. 1B, C |
SteadyStateTest.m |
| El Nino Biomass Burning perturbation tests | Fig. 2 | SyntheticEnso.m |
| Idealized inversion with synthetic concentrations | Fig. 3 |
SyntheticInversion.m |
| Inversions constrained by Observations | Fig. 4 | DriverScriptCaseI.m |

# Table 2: Inversions from Fig. 4 of Nguyen et al

Please refer to Tables 1 and A1 in Nguyen et al.

| Case | File Name |
| --- | --- |
| -I +[OH] | DriverScript.m |
| -I | DriverScriptCase1.m |
| +I | DriverScriptCase2.m |
| +I +S<sub>CO</sub> | DriverScriptCase3.m|
| +I +S<sub>OH</sub> | DriverScriptCase4.m |
| +I +S<sub>OH</sub> +S<sub>CO</sub> | DriverScriptCase5.m |
